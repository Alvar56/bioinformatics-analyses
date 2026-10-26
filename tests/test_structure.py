import numpy as np

from common.structure import kabsch, match_by_sequence, rmsd_all, superpose


def random_rotation(rng):
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q * np.sign(np.linalg.det(q))


def test_kabsch_recovers_a_rigid_motion():
    rng = np.random.default_rng(0)
    P = rng.normal(size=(50, 3)) * 10
    R = random_rotation(rng)
    Q = P @ R.T + np.array([3.0, -2.0, 7.5])
    R_fit, t_fit = kabsch(P, Q)
    assert np.allclose(R_fit, R, atol=1e-8)
    assert rmsd_all(P, Q) < 1e-8


def test_superpose_rejects_a_moved_loop():
    rng = np.random.default_rng(1)
    P = rng.normal(size=(100, 3)) * 10
    Q = P.copy()
    Q[90:] += 15.0                       # a loop that moved
    d, core, rmsd = superpose(P, Q)
    assert core[:90].all() and not core[90:].any()
    assert rmsd < 1e-6


def test_alignment_identity():
    pairs, ident = match_by_sequence("MKTAYIAKQR", "MKTAYLAKQR")
    assert len(pairs) == 10
    assert ident == 90.0
