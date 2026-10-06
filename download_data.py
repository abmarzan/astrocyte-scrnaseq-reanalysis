"""Download the two public GEO files this analysis needs (about 12.5 MB) and check them.

Source: NCBI GEO series GSE114000 (Batiuk, Martirosyan et al., Nature Communications 2020).
Usage:  python download_data.py
"""
import hashlib
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114000/suppl/"
FILES = {
    "GSE114000_Counts_Batiuk_Martirosyan_Supplementary_Data3.tsv.gz": "8283f391d38a9d5ae7a17395138dfa4a0851bdc40fc6bbf9024db12d6d458906",
    "GSE114000_Metadata_Batiuk_Martirosyan_Supplementary_Data1.xlsx": "fd6e5f5dfe2109d2667a6055a5d02387d58c2fab128bf1edb8d886c0d8be5624",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    data = Path("data")
    data.mkdir(exist_ok=True)
    for name, expected in FILES.items():
        target = data / name
        if not target.exists():
            print("downloading", name)
            urllib.request.urlretrieve(BASE + name, target)
        got = sha256(target)
        status = "ok" if got == expected else "CHECKSUM DIFFERS (the file on GEO may have been updated)"
        print(f"{name}: {target.stat().st_size:,} bytes, {status}")


if __name__ == "__main__":
    main()
