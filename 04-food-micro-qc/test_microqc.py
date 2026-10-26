import math

import pytest

from microqc import Plate, colony_count, mpn, mpn_interval, round_iso

TUBES = (3, 3, 3)
AMOUNTS = (0.1, 0.01, 0.001)


def test_iso7218_worked_example():
    # ISO 7218 worked example: 168 and 215 colonies at 10^-2, 14 and 25 at 10^-3
    plates = [Plate(1e-2, 168), Plate(1e-2, 215), Plate(1e-3, 14), Plate(1e-3, 25)]
    n = colony_count(plates)
    assert n == pytest.approx(422 / 0.022)
    assert round_iso(n) == "1.9 × 10^4"


def test_plates_outside_limits_are_ignored():
    plates = [Plate(1e-1, 1200), Plate(1e-2, 150), Plate(1e-3, 6)]
    assert colony_count(plates) == pytest.approx(150 / 1e-2)


def test_no_countable_plate_raises():
    with pytest.raises(ValueError):
        colony_count([Plate(1e-1, 4), Plate(1e-2, 0)])


@pytest.mark.parametrize("positives, expected", [
    # FDA Bacteriological Analytical Manual, Appendix 2, Table 1 (3 tubes x 0.1/0.01/0.001 g)
    ((1, 0, 0), 3.6),
    ((2, 0, 0), 9.2),
    ((3, 0, 0), 23),
    ((3, 1, 0), 43),
    ((3, 2, 0), 93),
    ((3, 3, 0), 240),
    ((3, 3, 1), 460),
    ((3, 3, 2), 1100),
])
def test_mpn_matches_bam_table(positives, expected):
    assert mpn(positives, TUBES, AMOUNTS) == pytest.approx(expected, rel=0.01)


def test_mpn_edge_cases():
    assert mpn((0, 0, 0), TUBES, AMOUNTS) == 0
    assert math.isinf(mpn((3, 3, 3), TUBES, AMOUNTS))


def test_interval_brackets_estimate():
    lo, hi = mpn_interval((3, 2, 0), TUBES, AMOUNTS)
    m = mpn((3, 2, 0), TUBES, AMOUNTS)
    assert lo < m < hi
