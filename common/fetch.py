"""Download helpers with a local cache, so every analysis re-runs offline once its
inputs have been fetched."""
from pathlib import Path

import requests

CACHE = Path(__file__).resolve().parent.parent / "data"


def download(url: str, name: str) -> Path:
    """Fetch `url` into data/<name> unless it is already there."""
    CACHE.mkdir(exist_ok=True)
    path = CACHE / name
    if not path.exists():
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        path.write_bytes(r.content)
    return path


def pdb_entry(code: str) -> Path:
    code = code.lower()
    return download(f"https://files.rcsb.org/download/{code}.cif", f"{code}.cif")


def alphafold_model(accession: str) -> Path:
    """Latest AlphaFold DB model (mmCIF) for a UniProt accession."""
    meta = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}", timeout=60)
    meta.raise_for_status()
    url = meta.json()[0]["cifUrl"]
    return download(url, f"AF-{accession}.cif")


def uniprot_entry(accession: str) -> dict:
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json", timeout=60)
    r.raise_for_status()
    return r.json()
