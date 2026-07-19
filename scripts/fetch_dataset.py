"""Crash-safe FASTQ downloader for the project dataset (default: SRR23872596).

Downloads paired-end FASTQ files directly from EBI-ENA over HTTPS and (by default)
decompresses them to plain .fastq that ``run_pipeline.py`` expects.

Design for large downloads over an unreliable connection
--------------------------------------------------------
Each file is streamed to a temporary ``<name>.part`` while its md5 is computed.
Only after the full download matches ENA's reported size AND md5 is the ``.part``
renamed to its final name. Two failure modes are handled differently:

* A file that already exists with the correct size + md5 is skipped (done).
* **Transient blip while the run is live** (brief internet outage, stalled socket):
  a socket timeout breaks the stall, then the file RESUMES from where it stopped
  via an HTTP Range request (ENA supports 206) — no wasted re-download. Retries a
  few times with a short wait, verifying md5 over the whole file at the end.
* **Power-off / process killed:** on the NEXT run the leftover ``<name>.part`` is
  RESUMED via HTTP Range and the whole file is md5-verified at the end (a corrupt
  partial fails the check and is re-downloaded from scratch). Finished files are skipped.

The ``--acc`` argument accepts a single run (SRRxxxx) OR a whole BioProject
(PRJNAxxxxxx) — ENA returns every run in the project, and all of them are fetched.

Usage
-----
    # One run into the repo root (download + gunzip):
    .venv/bin/python scripts/fetch_dataset.py

    # The ENTIRE Brahmaputra BioProject as .gz onto the HDD (recommended for the full set):
    .venv/bin/python scripts/fetch_dataset.py --acc PRJNA944918 \
        --outdir /mnt/shinrinyoku/brahmaputra_metagenome --no-gunzip

    .venv/bin/python scripts/fetch_dataset.py --no-gunzip     # keep only .fastq.gz
    .venv/bin/python scripts/fetch_dataset.py --keep-gz       # keep .gz after gunzip

Just re-run the exact same command after any interruption — it resumes the interrupted
file from where it stopped (md5-verified) and skips already-finished files.
"""

import argparse
import gzip
import hashlib
import shutil
import sys
import time
from pathlib import Path
from urllib.request import Request, urlopen

ENA_FILEREPORT = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={acc}&result=read_run"
    "&fields=run_accession,sample_title,fastq_ftp,fastq_bytes,fastq_md5&format=tsv"
)
CHUNK = 1024 * 1024  # 1 MiB
MAX_ATTEMPTS = 6  # per file
# If no data arrives for this many seconds (e.g. the connection stalls during a
# brief internet outage), the read raises instead of hanging forever, so the
# retry loop can kick in. Without this the download can freeze indefinitely.
SOCKET_TIMEOUT = 60
RETRY_WAIT = 10  # seconds to wait before a retry, giving the network time to recover


class FileUnavailableError(RuntimeError):
    """The host served something other than the file — e.g. an HTML directory
    listing where the .fastq.gz should be. This happens when a run's object is
    missing/broken on the ENA mirror (occasionally a file gets replaced by an
    empty directory during a re-sync). Retrying within the same run won't help,
    so this is raised to fail this file fast instead of burning every attempt.
    """


def http_get_text(url):
    req = Request(url, headers={"User-Agent": "plastic-metagenome-pipeline/1.0"})
    with urlopen(req, timeout=SOCKET_TIMEOUT) as resp:
        return resp.read().decode("utf-8")


def query_ena(accession):
    """Return list of file dicts for a run (SRRxxx) or every run in a project (PRJNAxxx).

    Each dict: {run, sample, url, name, size, md5}.
    """
    text = http_get_text(ENA_FILEREPORT.format(acc=accession))
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"ENA returned no file records for {accession}:\n{text}")
    header = lines[0].split("\t")
    files = []
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t")))
        run = row.get("run_accession", "")
        sample = row.get("sample_title", "")
        ftps = row.get("fastq_ftp", "").split(";")
        sizes = row.get("fastq_bytes", "").split(";")
        md5s = row.get("fastq_md5", "").split(";")
        for ftp, size, md5 in zip(ftps, sizes, md5s):
            ftp = ftp.strip()
            if not ftp:
                continue
            # ENA lists an FTP host path; fetch over HTTPS from the same host.
            files.append(
                {
                    "run": run,
                    "sample": sample,
                    "url": "https://" + ftp,
                    "name": ftp.rsplit("/", 1)[-1],
                    "size": int(size) if size.strip() else None,
                    "md5": md5.strip() or None,
                }
            )
    if not files:
        raise RuntimeError(f"No fastq_ftp entries parsed for {accession}")
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

    # A leftover .part (from a power-off, a kill, or an earlier blip) is RESUMED, not
    # discarded: the loop continues it via HTTP Range, and the final md5 check over the
    # whole file guarantees integrity — a corrupt partial fails md5 and is re-downloaded.
    for attempt in range(1, MAX_ATTEMPTS + 1):
        resume_from = part_path.stat().st_size if part_path.exists() else 0

        # Already fully on disk but not yet verified/renamed (e.g. killed just before
        # the rename): verify and finish, don't re-request past EOF.
        if expected_size is not None and resume_from >= expected_size:
            if resume_from == expected_size and (
                not expected_md5 or md5_of_file(part_path) == expected_md5
            ):
                part_path.rename(final_path)
                print(f"[done] {name} verified ({human(expected_size)}).")
                return final_path
            part_path.unlink()
            resume_from = 0

        h = hashlib.md5()
        headers = {"User-Agent": "plastic-metagenome-pipeline/1.0"}
        if resume_from:
            # Re-hash bytes already on disk so the final md5 covers the whole file.
            with open(part_path, "rb") as f:
                for chunk in iter(lambda: f.read(CHUNK), b""):
                    h.update(chunk)
            headers["Range"] = f"bytes={resume_from}-"
            mode = "ab"
            print(f"[resume] {name} from {human(resume_from)}/{human(expected_size)} "
                  f"(attempt {attempt}/{MAX_ATTEMPTS})", flush=True)
        else:
            mode = "wb"
            print(f"[get] {name} ({human(expected_size)}) attempt {attempt}/{MAX_ATTEMPTS}",
                  flush=True)

        downloaded = resume_from
        try:
            req = Request(url, headers=headers)
            with urlopen(req, timeout=SOCKET_TIMEOUT) as resp:
                # If we asked to resume but the server ignored Range (200, not 206),
                # it resends the whole file; appending would corrupt it, so restart.
                status = getattr(resp, "status", None) or resp.getcode()
                # A real .fastq.gz is served as gzip/octet-stream. If the host hands
                # back an HTML page, urllib has followed a redirect to an error page or
                # a directory listing standing in for the missing file — downloading it
                # would just yield a few hundred bytes that fail the size check on every
                # attempt. Bail out immediately with a clear message instead.
                ctype = resp.headers.get_content_type()
                if ctype.startswith("text/"):
                    raise FileUnavailableError(
                        f"host returned {ctype} ({resp.geturl()}), not the file — "
                        "the object is missing/broken on the ENA mirror right now"
                    )
                if resume_from and status != 206:
                    print("[warn] server ignored resume; restarting this file from scratch")
                    part_path.unlink()
                    h = hashlib.md5()
                    downloaded = 0
                    mode = "wb"
                with open(part_path, mode) as out:
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
        except FileUnavailableError:
            # Not a network blip — the file isn't really there. Don't retry.
            raise
        except Exception as exc:  # stall/timeout/network drop
            print(f"\n[warn] connection problem: {exc}")
            # Keep the partial so the next attempt RESUMES (no wasted re-download).
            if attempt == MAX_ATTEMPTS:
                raise
            print(f"[retry] waiting {RETRY_WAIT}s, then resuming from "
                  f"{human(part_path.stat().st_size if part_path.exists() else 0)} ...")
            time.sleep(RETRY_WAIT)
            continue

        # Verify size and md5 before accepting.
        if expected_size is not None and downloaded != expected_size:
            print(f"[warn] size mismatch ({downloaded} != {expected_size}); restarting file.")
            if part_path.exists():
                part_path.unlink()
            continue
        if expected_md5 and h.hexdigest() != expected_md5:
            print("[warn] md5 mismatch; restarting file.")
            if part_path.exists():
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
    parser.add_argument(
        "--acc",
        "--srr",
        dest="acc",
        default="SRR23872596",
        help="SRA run (SRRxxxx) or whole BioProject (PRJNAxxxxxx)",
    )
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

    print(f"Querying ENA for {args.acc} ...")
    files = query_ena(args.acc)
    runs = sorted({f["run"] for f in files})
    total = sum(f["size"] or 0 for f in files)
    print(
        f"Found {len(runs)} run(s), {len(files)} file(s), "
        f"total {human(total)} compressed:"
    )
    for f in files:
        loc = (f["sample"] or "").split("_")[0]
        print(f"  - {f['name']:28} {human(f['size']):>9}  {loc}")
    print()

    gz_paths = []
    failed = []
    for idx, f in enumerate(files, 1):
        print(f"[{idx}/{len(files)}] {f['run']}", end="  ")
        try:
            gz_paths.append(download_one(f, outdir))
        except Exception as exc:
            # One bad file (e.g. missing on the ENA mirror) must not abort the whole
            # batch — record it, keep going, and report every failure at the end.
            print(f"[skip-file] {f['name']}: {exc}")
            failed.append((f["name"], str(exc)))

    if args.gunzip:
        print()
        for gz in gz_paths:
            gunzip_to_fastq(gz, keep_gz=args.keep_gz)

    if failed:
        print(f"\n[warning] {len(failed)} of {len(files)} file(s) could not be downloaded:")
        for name, why in failed:
            print(f"  - {name}: {why}")
        print("Re-run the same command later to retry only the missing files "
              "(finished files are skipped). If a file stays broken on ENA, fetch that "
              "run from NCBI instead (SRA .sra -> fasterq-dump).")
        sys.exit(1)

    print(f"\nAll dataset files ready in {outdir}/")
    if len(runs) == 1 and args.gunzip:
        r = runs[0]
        prefix = "" if args.outdir == "." else args.outdir.rstrip("/") + "/"
        print("Next: run the pipeline, e.g.")
        print(
            f"  .venv/bin/python run_pipeline.py "
            f"--fastq1 {prefix}{r}_1.fastq --fastq2 {prefix}{r}_2.fastq"
        )
    else:
        print(f"Downloaded {len(runs)} runs. Process each run through run_pipeline.py "
              "(one accession at a time), pointing --fastq1/--fastq2 at its files.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[interrupted] Re-run the same command; it resumes the interrupted "
              "file from where it stopped and skips finished files.")
        sys.exit(130)
    except Exception as exc:
        print(f"[fetch_dataset] Error: {exc}", file=sys.stderr)
        sys.exit(1)

# downloadaded some data manually later