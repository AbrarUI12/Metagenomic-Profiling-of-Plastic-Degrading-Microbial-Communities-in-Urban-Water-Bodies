import argparse
import csv
import platform
import sys
import tarfile
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.utils import ensure_dir, log_command, run_cmd, which


def read_family_map(path):
    mapping = {}
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mapping[row["protein_id"]] = row["family"]
    return mapping


def filter_hits(in_path, out_path, family_map, evalue, pident, min_aln_len):
    ensure_dir(Path(out_path).parent)
    with open(in_path, "r", encoding="utf-8") as f_in, open(
        out_path, "w", encoding="utf-8", newline=""
    ) as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(
            ["read_id", "protein_id", "family", "evalue", "bitscore", "pident", "aln_len"]
        )
        for line in f_in:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            read_id, protein_id, evalue_s, bitscore_s, pident_s, aln_len_s = parts[:6]
            protein_id = normalize_protein_id(protein_id)
            try:
                ev = float(evalue_s)
                pid = float(pident_s)
                aln_len = int(float(aln_len_s))
            except ValueError:
                continue
            if ev > evalue or pid < pident or aln_len < min_aln_len:
                continue
            family = family_map.get(protein_id, "unknown")
            writer.writerow(
                [read_id, protein_id, family, ev, bitscore_s, pid, aln_len]
            )


def prepare_fastq_subset(in_fastq, out_fastq, max_reads):
    ensure_dir(Path(out_fastq).parent)
    reads = 0
    with open(in_fastq, "r", encoding="utf-8", errors="replace") as fin, open(
        out_fastq, "w", encoding="utf-8"
    ) as fout:
        while reads < max_reads:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not qual:
                break
            fout.write(header)
            fout.write(seq)
            fout.write(plus)
            fout.write(qual)
            reads += 1
    return out_fastq


def conda_install_diamond(commands_log=None, skip=False):
    if skip:
        return False
    if not which("conda"):
        return False
    cmd = "conda install -y -c bioconda diamond"
    run_cmd(cmd, commands_log=commands_log, check=False)
    return which("diamond") is not None


DIAMOND_VERSION = "v2.2.4"


def diamond_binary_name():
    return "diamond.exe" if platform.system() == "Windows" else "diamond"


def download_diamond(commands_log=None):
    # Best-effort, platform-aware download. If it fails, caller falls back.
    system = platform.system()
    if system == "Windows":
        asset = "diamond-windows.zip"
    elif system == "Darwin":
        asset = "diamond-macos.tar.gz"
    else:
        asset = "diamond-linux64.tar.gz"
    url = (
        "https://github.com/bbuchfink/diamond/releases/download/"
        f"{DIAMOND_VERSION}/{asset}"
    )
    dest = Path("tools") / asset
    out_dir = Path("tools/diamond")
    ensure_dir(out_dir)
    try:
        if commands_log:
            log_command(commands_log, f"DOWNLOAD {url} -> {dest}")
        with urlopen(url) as resp:
            dest.write_bytes(resp.read())
        if asset.endswith(".zip"):
            with ZipFile(dest, "r") as zf:
                zf.extractall(out_dir)
        else:
            with tarfile.open(dest, "r:gz") as tf:
                for member in tf.getmembers():
                    # Flatten: extract only the diamond binary into out_dir
                    if Path(member.name).name == "diamond":
                        member.name = "diamond"
                        tf.extract(member, out_dir)
        binary = out_dir / diamond_binary_name()
        if binary.exists():
            binary.chmod(0o755)
        return binary.exists()
    except Exception:
        return False


def find_diamond():
    diamond = which("diamond")
    if diamond:
        return diamond
    local = Path("tools/diamond") / diamond_binary_name()
    if local.exists():
        return str(local)
    return None


def run_diamond(fastq, db_prefix, out_tsv, threads, evalue, pident, commands_log, skip_conda):
    diamond = find_diamond()
    if not diamond:
        if not conda_install_diamond(commands_log, skip=skip_conda):
            download_diamond(commands_log)
        diamond = find_diamond()
    if not diamond:
        return False
    db_cmd = f"\"{diamond}\" makedb --in \"{db_prefix}.fasta\" -d \"{db_prefix}\""
    run_cmd(db_cmd, commands_log=commands_log)
    cmd = (
        f"\"{diamond}\" blastx -d \"{db_prefix}\" -q \"{fastq}\" "
        f"-o \"{out_tsv}\" --evalue {evalue} --id {pident} "
        f"--max-target-seqs 1 --threads {threads} "
        "--outfmt 6 qseqid sseqid evalue bitscore pident length"
    )
    run_cmd(cmd, commands_log=commands_log)
    return True


def run_blastx(fastq, fasta_db, out_tsv, threads, evalue, commands_log):
    if not which("makeblastdb") or not which("blastx"):
        return False
    db_prefix = Path(fasta_db).with_suffix("").as_posix() + "_blastdb"
    cmd_db = f"makeblastdb -dbtype prot -in \"{fasta_db}\" -out \"{db_prefix}\""
    run_cmd(cmd_db, commands_log=commands_log)
    cmd = (
        f"blastx -query \"{fastq}\" -db \"{db_prefix}\" "
        f"-out \"{out_tsv}\" -evalue {evalue} "
        f"-max_target_seqs 1 -num_threads {threads} "
        "-outfmt \"6 qseqid sseqid evalue bitscore pident length\""
    )
    run_cmd(cmd, commands_log=commands_log)
    return True


CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def revcomp(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def translate(seq):
    seq = seq.upper()
    pep = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        pep.append(CODON_TABLE.get(codon, "X"))
    return "".join(pep)


def kmer_index_from_fasta(fasta_path, k):
    index = {}
    with open(fasta_path, "r", encoding="utf-8") as f:
        header = None
        seq_chunks = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    seq = "".join(seq_chunks)
                    add_kmers(index, header, seq, k)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header:
            seq = "".join(seq_chunks)
            add_kmers(index, header, seq, k)
    return index


def add_kmers(index, header, seq, k):
    acc = header.split("|")[1] if "|" in header else header.split()[0]
    seq = seq.replace("*", "")
    if len(seq) < k:
        return
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if kmer not in index:
            index[kmer] = set()
        index[kmer].add(acc)


def kmer_fallback_search(fastq, fasta_db, out_tsv, family_map, k=7, commands_log=None):
    if commands_log:
        log_command(commands_log, f"KMER_FALLBACK {fastq} vs {fasta_db} k={k}")
    index = kmer_index_from_fasta(fasta_db, k)
    ensure_dir(Path(out_tsv).parent)
    with open(fastq, "r", encoding="utf-8", errors="replace") as fin, open(
        out_tsv, "w", encoding="utf-8", newline=""
    ) as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(
            ["read_id", "protein_id", "family", "evalue", "bitscore", "pident", "aln_len"]
        )
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline().strip()
            fin.readline()
            fin.readline()
            read_id = header.strip().split()[0][1:]
            hits = find_kmer_hits(seq, index, k)
            if not hits:
                continue
            protein_id = sorted(hits)[0]
            family = family_map.get(protein_id, "unknown")
            writer.writerow([read_id, protein_id, family, 0.0, 0.0, 100.0, k])


def find_kmer_hits(seq, index, k):
    hits = set()
    frames = [seq, seq[1:], seq[2:], revcomp(seq), revcomp(seq)[1:], revcomp(seq)[2:]]
    for frame in frames:
        pep = translate(frame)
        for chunk in pep.split("*"):
            if len(chunk) < k:
                continue
            for i in range(len(chunk) - k + 1):
                kmer = chunk[i : i + k]
                if kmer in index:
                    hits.update(index[kmer])
                    if len(hits) > 0:
                        return hits
    return hits


def normalize_protein_id(protein_id):
    if "|" in protein_id:
        parts = protein_id.split("|")
        if len(parts) >= 2:
            return parts[1]
    return protein_id


def merge_outputs(outputs, merged):
    ensure_dir(Path(merged).parent)
    with open(merged, "w", encoding="utf-8", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(
            ["read_id", "protein_id", "family", "evalue", "bitscore", "pident", "aln_len"]
        )
        for path in outputs:
            with open(path, "r", encoding="utf-8") as f_in:
                next(f_in, None)
                for line in f_in:
                    f_out.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq1", required=True)
    parser.add_argument("--fastq2", default=None)
    parser.add_argument("--fasta-db", default="refs/plastic_enzymes.fasta")
    parser.add_argument("--family-map", default="refs/enzyme_family_map.tsv")
    parser.add_argument("--out", default="outputs/enzyme_hits/enzyme_hits.tsv")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--pident", type=float, default=30.0)
    parser.add_argument("--min-aln-len", type=int, default=50)
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--commands-log", default="outputs/logs/commands.log")
    parser.add_argument("--method-log", default="outputs/logs/search_method.txt")
    parser.add_argument("--skip-conda-install", action="store_true")
    args = parser.parse_args()

    family_map = read_family_map(args.family_map)

    fastq_inputs = []
    for label, fpath in [("read1", args.fastq1), ("read2", args.fastq2)]:
        if not fpath:
            continue
        fpath = Path(fpath)
        if args.max_reads and args.max_reads > 0:
            subset = Path("outputs/tmp") / f"{fpath.stem}_subset.fastq"
            prepare_fastq_subset(str(fpath), str(subset), args.max_reads)
            fastq_inputs.append((label, str(subset)))
        else:
            fastq_inputs.append((label, str(fpath)))

    outputs = []
    for label, fq in fastq_inputs:
        raw_out = Path("outputs/enzyme_hits") / f"raw_{label}.tsv"
        filtered_out = Path("outputs/enzyme_hits") / f"filtered_{label}.tsv"
        ensure_dir(raw_out.parent)
        db_prefix = Path(args.fasta_db).with_suffix("")
        method = ""
        diamond_ok = run_diamond(
            fq,
            str(db_prefix),
            str(raw_out),
            args.threads,
            args.evalue,
            args.pident,
            args.commands_log,
            args.skip_conda_install,
        )
        if not diamond_ok:
            blast_ok = run_blastx(
                fq, args.fasta_db, str(raw_out), args.threads, args.evalue, args.commands_log
            )
            if not blast_ok:
                kmer_fallback_search(
                    fq,
                    args.fasta_db,
                    str(filtered_out),
                    family_map,
                    k=max(7, args.min_aln_len),
                    commands_log=args.commands_log,
                )
                method = "kmer_fallback"
                outputs.append(str(filtered_out))
                log_command(args.method_log, f"{label}\t{method}\t{fq}")
                continue
            method = "blastx"
        else:
            method = "diamond"
        filter_hits(
            str(raw_out),
            str(filtered_out),
            family_map,
            args.evalue,
            args.pident,
            args.min_aln_len,
        )
        outputs.append(str(filtered_out))
        log_command(args.method_log, f"{label}\t{method}\t{fq}")

    merge_outputs(outputs, args.out)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[search] Error: {exc}", file=sys.stderr)
        sys.exit(1)
