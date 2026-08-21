#!/usr/bin/env python3
"""Download the public SPARC all-point RAR table used by the diagnostic.

The dataset is maintained by the SPARC collaboration. This helper does not
change the upstream terms; users should cite the SPARC/RAR papers listed in
README.md and verify current upstream data-use guidance.
"""
from __future__ import annotations
import argparse, hashlib, urllib.request
from pathlib import Path

URL = "https://astroweb.case.edu/SPARC/RAR.mrt"
EXPECTED_SHA256 = "24aa7059dab7fa44787f7c11191052489899819370f6508621674769f3b72833"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=root / "data" / "RAR.mrt")
    ap.add_argument("--allow-hash-change", action="store_true",
                    help="Keep the file if upstream content has changed; print the new hash.")
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(args.out.suffix + ".part")
    urllib.request.urlretrieve(URL, tmp)
    got = sha256(tmp)
    if got != EXPECTED_SHA256 and not args.allow_hash_change:
        tmp.unlink(missing_ok=True)
        raise SystemExit(
            "Downloaded SPARC table hash differs from the manuscript snapshot.\n"
            f"expected: {EXPECTED_SHA256}\nreceived: {got}\n"
            "Inspect upstream changes or rerun with --allow-hash-change."
        )
    tmp.replace(args.out)
    print(f"saved: {args.out}")
    print(f"sha256: {got}")

if __name__ == "__main__":
    main()
