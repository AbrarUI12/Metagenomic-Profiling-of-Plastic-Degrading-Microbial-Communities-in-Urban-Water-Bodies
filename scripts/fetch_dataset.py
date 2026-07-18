"""Crash-safe FASTQ downloader for the project dataset (default: SRR23872596).

Downloads paired-end FASTQ files directly from EBI-ENA over HTTPS and (by default)
decompresses them to plain .fastq that ``run_pipeline.py`` expects.

Design for large, resumable-across-power-loss downloads
------------------------------------------------------
Each file is streamed to a temporary ``<name>.part`` while its md5 is computed.
Only after the full download matches ENA's reported size AND md5 is the ``.part``
renamed to its final name. Therefore:

* A file that already exists with the correct size + md5 is skipped (done).
* If the machine powers off mid-download, only a ``<name>.part`` is left behind.
  On the next run that partial file is DELETED and re-downloaded from scratch
  (no byte-range resume), then the remaining files continue. This matches the
  requirement: whatever was mid-download restarts cleanly; finished files are kept.

Usage
-----
    .venv/bin/python scripts/fetch_dataset.py                 # download + gunzip into repo root
    .venv/bin/python scripts/fetch_dataset.py --outdir data   # into data/
    .venv/bin/python scripts/fetch_dataset.py --no-gunzip     # keep only .fastq.gz
    .venv/bin/python scripts/fetch_dataset.py --keep-gz       # keep .gz after gunzip

Just re-run the exact same command after any interruption; it resumes safely.
"""

import argparse
import gzip
import hashlib
import shutil
import sys
from pathlib import Path
from urllib.request import Request, urlopen

ENA_FILEREPORT = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={srr}&result=read_run"
    "&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5&format=tsv"
)
CHUNK = 1024 * 1024  # 1 MiB
MAX_ATTEMPTS = 4  # per file; each attempt restarts from scratch (no resume)


def http_get_text(url):
    req = Request(url, headers={"User-Agent": "plastic-metagenome-pipeline/1.0"})
    with urlopen(req) as resp:
        return resp.read().decode("utf-8")


def query_ena(srr):
    """Return list of dicts: {url, name, size, md5} for each fastq file."""
    text = http_get_text(ENA_FILEREPORT.format(srr=srr))
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"ENA returned no file records for {srr}:\n{text}")
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    ftps = row.get("fastq_ftp", "").split(";")
    sizes = row.get("fastq_bytes", "").split(";")
    md5s = row.get("fastq_md5", "").split(";")
    files = []
    for ftp, size, md5 in zip(ftps, sizes, md5s):
        ftp = ftp.strip()
        if not ftp:
            continue
        # ENA lists an FTP host path; fetch over HTTPS from the same host.
        url = "https://" + ftp
        files.append(
            {
                "url": url,
                "name": ftp.rsplit("/", 1)[-1],
                "size": int(size) if size.strip() else None,
                "md5": md5.strip() or None,
            }
        )
    if not files:
        raise RuntimeError(f"No fastq_ftp entries parsed for {srr}")
    return files


def md5_of_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def human(n):
    if n is None:
        return "?"
    for unit in ["B", "KB", "MB", "GB"]:
        if n < 1024:
            return f"{n:.0f}{unit}" if unit == "B" else f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}TB"


def already_complete(final_path, expected_size, expected_md5):
    if not final_path.exists():
        return False
    if expected_size is not None and final_path.stat().st_size != expected_size:
        return False
    if expected_md5:
        print(f"  [verify] checking md5 of existing {final_path.name} ...", flush=True)
        if md5_of_file(final_path) != expected_md5:
            print("  [verify] md5 mismatch on existing file; will re-download.")
            return False
    return True


def download_one(file_info, outdir):
    name = file_info["name"]
    url = file_info["url"]
    expected_size = file_info["size"]
    expected_md5 = file_info["md5"]
    final_path = outdir / name
    part_path = outdir / (name + ".part")

    if already_complete(final_path, expected_size, expected_md5):
        print(f"[skip] {name} already present and verified ({human(expected_size)}).")
        return final_path

    for attempt in range(1, MAX_ATTEMPTS + 1):
        # Always start fresh: delete any leftover partial from a prior crash/attempt.
        if part_path.exists():
            print(f"[clean] removing stale partial {part_path.name}")
            part_path.unlink()

        print(
            f"[get] {name} ({human(expected_size)}) attempt {attempt}/{MAX_ATTEMPTS}",
            flush=True,
        )
        h = hashlib.md5()
        downloaded = 0
        try:
            req = Request(url, headers={"User-Agent": "plastic-metagenome-pipeline/1.0"})
            with urlopen(req) as resp, open(part_path, "wb") as out:
                while True:
                    chunk = resp.read(CHUNK)
                    if not chunk:
                        break
                    out.write(chunk)
                    h.update(chunk)
                    downloaded += len(chunk)
                    if expected_size:
                        pct = downloaded / expected_size * 100
                        print(
                            f"\r    {human(downloaded)}/{human(expected_size)} ({pct:5.1f}%)",
                            end="",
                            flush=True,
                        )
            print()
        except Exception as exc:  # network error, interrupted, etc.
            print(f"\n[warn] download failed: {exc}")
            if part_path.exists():
                part_path.unlink()
            if attempt == MAX_ATTEMPTS:
                raise
            continue

        # Verify size and md5 before accepting.
        if expected_size is not None and downloaded != expected_size:
            print(f"[warn] size mismatch ({downloaded} != {expected_size}); retrying.")
            part_path.unlink()
            continue
        if expected_md5 and h.hexdigest() != expected_md5:
            print(f"[warn] md5 mismatch; retrying.")
            part_path.unlink()
            continue

        part_path.rename(final_path)
        print(f"[done] {name} verified ({human(expected_size)}).")
        return final_path

    raise RuntimeError(f"Failed to download {name} after {MAX_ATTEMPTS} attempts.")


def gunzip_to_fastq(gz_path, keep_gz):
    """Decompress <name>.fastq.gz -> <name>.fastq, crash-safe via .part."""
    if gz_path.suffix != ".gz":
        return gz_path
    final = gz_path.with_suffix("")  # drop .gz
    part = final.with_name(final.name + ".part")
    if final.exists():
        print(f"[skip] {final.name} already decompressed.")
        if not keep_gz and gz_path.exists():
            gz_path.unlink()
        return final
    if part.exists():
        print(f"[clean] removing stale {part.name}")
        part.unlink()
    print(f"[gunzip] {gz_path.name} -> {final.name}", flush=True)
    with gzip.open(gz_path, "rb") as fin, open(part, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=CHUNK)
    part.rename(final)
    print(f"[done] {final.name} ({human(final.stat().st_size)}).")
    if not keep_gz:
        gz_path.unlink()
        print(f"[clean] removed {gz_path.name} (use --keep-gz to retain).")
    return final


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--srr", default="SRR23872596", help="SRA run accession")
    parser.add_argument("--outdir", default=".", help="destination directory")
    parser.add_argument(
        "--no-gunzip",
        dest="gunzip",
        action="store_false",
        help="keep .fastq.gz only; do not decompress",
    )
    parser.add_argument(
        "--keep-gz",
        action="store_true",
        help="keep the .fastq.gz after decompressing",
    )
    parser.set_defaults(gunzip=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Querying ENA for {args.srr} ...")
    files = query_ena(args.srr)
    total = sum(f["size"] or 0 for f in files)
    print(f"Found {len(files)} file(s), total {human(total)} compressed:")
    for f in files:
        print(f"  - {f['name']}  {human(f['size'])}  md5={f['md5']}")
    print()

    gz_paths = []
    for f in files:
        gz_paths.append(download_one(f, outdir))

    if args.gunzip:
        print()
        for gz in gz_paths:
            gunzip_to_fastq(gz, keep_gz=args.keep_gz)

    print("\nAll dataset files ready.")
    print("Next: run the pipeline, e.g.")
    r1 = args.srr + "_1.fastq"
    r2 = args.srr + "_2.fastq"
    prefix = "" if args.outdir == "." else args.outdir.rstrip("/") + "/"
    print(f"  .venv/bin/python run_pipeline.py --fastq1 {prefix}{r1} --fastq2 {prefix}{r2}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[interrupted] Re-run the same command to resume "
              "(the in-progress file restarts; finished files are kept).")
        sys.exit(130)
    except Exception as exc:
        print(f"[fetch_dataset] Error: {exc}", file=sys.stderr)
        sys.exit(1)
