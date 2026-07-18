import gzip
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def open_text(path, mode="r"):
    """Open a text file transparently whether or not it is gzip-compressed.

    Detects by the ``.gz`` suffix. Reads use errors='replace' so a stray byte in
    a large FASTQ never aborts a run. Used so the pipeline can consume .fastq.gz
    directly without a separate decompression step (DIAMOND already reads .gz).
    """
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_json(path, obj):
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True)


def append_text(path, text):
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("a", encoding="utf-8") as f:
        f.write(text)


def log_command(commands_log, cmd):
    timestamp = utc_now()
    append_text(commands_log, f"[{timestamp}] {cmd}\n")


def which(cmd):
    return shutil.which(cmd)


def run_cmd(cmd, commands_log=None, check=True):
    if commands_log:
        log_command(commands_log, cmd)
    return subprocess.run(cmd, shell=True, check=check)

