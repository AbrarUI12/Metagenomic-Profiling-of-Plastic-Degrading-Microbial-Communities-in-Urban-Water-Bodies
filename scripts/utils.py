import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


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

