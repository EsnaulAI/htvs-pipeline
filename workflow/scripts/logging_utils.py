import os
import sys


def log_info(message):
    print(f"[INFO] {message}")


def log_warn(message):
    print(f"[WARN] {message}")


def log_error(message):
    print(f"[ERROR] {message}")


def require_file(path, label="archivo de entrada"):
    if not path or not os.path.isfile(path):
        log_error(f"Falta {label}: {path}")
        sys.exit(1)


def ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def confirm_file(path, label="archivo de salida"):
    if not os.path.isfile(path):
        log_error(f"No se generó {label}: {path}")
        sys.exit(1)


def confirm_path(path, label="salida"):
    if not os.path.exists(path):
        log_error(f"No se generó {label}: {path}")
        sys.exit(1)
