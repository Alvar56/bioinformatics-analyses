"""Small structural-bioinformatics helpers built on Biopython: CA extraction, sequence-based
residue matching and outlier-rejecting superposition."""
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one

_parser = MMCIFParser(QUIET=True)


def ca_trace(path: Path, chain: str, keep=None):
    """Return (resnums, sequence, CA coords, B-factors) for one chain.
    `keep(resnum)` can filter residues, e.g. to drop a fusion protein."""
    model = next(iter(_parser.get_structure(path.stem, str(path))))
    nums, seq, xyz, bfac = [], [], [], []
    for res in model[chain]:
        if res.id[0] != " " or "CA" not in res:
            continue
        if keep is not None and not keep(res.id[1]):
            continue
        try:
            aa = index_to_one(three_to_index(res.get_resname()))
        except (KeyError, ValueError):
            continue
        nums.append(res.id[1])
        seq.append(aa)
        xyz.append(res["CA"].coord)
        bfac.append(res["CA"].bfactor)
    return np.array(nums), "".join(seq), np.array(xyz, float), np.array(bfac, float)


def match_by_sequence(seq_a: str, seq_b: str):
    """Global alignment (BLOSUM62, affine gaps). Returns index pairs (i, j) of aligned
    residues and the percent identity over those pairs."""
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    al.mode = "global"
    aln = al.align(seq_a, seq_b)[0]
    pairs = []
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        pairs += list(zip(range(a0, a1), range(b0, b1)))
    pairs = np.array(pairs)
    ident = np.mean([seq_a[i] == seq_b[j] for i, j in pairs]) * 100
    return pairs, ident


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Rotation R and translation t minimising |(P @ R.T + t) - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - pc @ R.T


def superpose(P: np.ndarray, Q: np.ndarray, cutoff: float = 4.0, rounds: int = 10):
    """Iterative superposition: fit, drop pairs further apart than `cutoff` angstrom, refit.
    Returns per-pair distances after the final fit, the core mask, and the core RMSD."""
    core = np.ones(len(P), bool)
    for _ in range(rounds):
        R, t = kabsch(P[core], Q[core])
        d = np.linalg.norm(P @ R.T + t - Q, axis=1)
        new = d < cutoff
        if new.sum() < 3 or (new == core).all():
            break
        core = new
    rmsd = float(np.sqrt(np.mean(d[core] ** 2)))
    return d, core, rmsd


def rmsd_all(P: np.ndarray, Q: np.ndarray) -> float:
    R, t = kabsch(P, Q)
    return float(np.sqrt(np.mean(np.sum((P @ R.T + t - Q) ** 2, axis=1))))
