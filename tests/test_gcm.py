# /mpcfit/tests/test_gcm.py

"""
# --------------------------------------------------------------
# Oct 2018
# Payne
#
# Test the functions in /mpcfit/mpcfit/gcm.py 
# - GCM class checks for deviation from GCM for a tracket
# - "GCM" == Great Circle Motion
#
# Straight port from work by Sonia Keys
# - Some of Sonia's notes are in /mpcfit/mpcfit/gcm.md
# - At present, this is a striaght C-port, so NOT pythonic
#   --- Could use a lot more numpy & get a lot of speed-up
#
# --------------------------------------------------------------
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import math
import pytest

# Importing of local modules/packages required for this test
# --------------------------------------------------------------
import mpcfit.phys_const as PHYS

# Import the specific package/module/function we are testing
# --------------------------------------------------------------
from mpcfit.gcm import gcm


# Basic test(s) of the GCM module
# ----------------------------------------------------------

def near(a, b, eps):
    return abs(a - b) <= max(eps, eps * max(abs(a), abs(b)))

def cartNear(a, b, eps):
    return near(a.x, b.x, eps) and near(a.y, b.y, eps) and near(a.z, b.z, eps)

def sphrNear(a, b, eps):
    return near(a.ra, b.ra, eps) and near(a.dec, b.dec, eps)


def test_cartCross():
    a = gcm.Cart(1, 0, 0)
    b = gcm.Cart(0, 1, 0)
    got = gcm.cartCross(a, b)
    want = gcm.Cart(0, 0, 1)
    eps = 1e-14
    assert cartNear(got, want, eps)


def test_cartDot():
    a = gcm.Cart(1, 2, 3)
    b = gcm.Cart(7, 2, 0)
    assert gcm.cartDot(a, b) == 11


def test_cartSquare():
    a = gcm.Cart(1, -2, 3)
    assert gcm.cartSquare(a) == 14


def test_cartRot():
    a = gcm.Cart(0, 1, 0)
    d30 = 30 * math.pi / 180
    s = math.sin(d30)
    c = math.cos(d30)
    m = [  # rotate around x axis
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]]
    got = gcm.cartRot(m, [a])[0]
    want = gcm.Cart(0, math.sqrt(3) / 2, .5)
    eps = 1e-14
    assert cartNear(got, want, eps)


def test_cartToSphr():
    s = [gcm.Sphr(0, 0),
        gcm.Sphr(math.pi / 2, 0),
        gcm.Sphr(0, math.pi / 2),
        gcm.Sphr(math.pi / 4, 0),
        gcm.Sphr(math.pi / 6, -math.pi / 4)]
    c = [gcm.Cart(1, 0, 0),
        gcm.Cart(0, 1, 0),
        gcm.Cart(0, 0, 1),
        gcm.Cart(math.sqrt(2) / 2, math.sqrt(2) / 2, 0),
        gcm.Cart(math.sqrt(6) / 4, math.sqrt(2) / 4, -math.sqrt(2) / 2)]
    eps = 1e-14
    for got, want in zip(gcm.sphrToCart(s), c):
        assert cartNear(got, want, eps)


def test_sphrToCart():
    s = [gcm.Sphr(0, 0),
        gcm.Sphr(math.pi / 2, 0),
        gcm.Sphr(0, math.pi / 2),
        gcm.Sphr(math.pi / 4, 0),
        gcm.Sphr(math.pi / 6, -math.pi / 4)]
    c = [gcm.Cart(1, 0, 0),
        gcm.Cart(0, 1, 0),
        gcm.Cart(0, 0, 1),
        gcm.Cart(math.sqrt(2) / 2, math.sqrt(2) / 2, 0),
        gcm.Cart(math.sqrt(6) / 4, math.sqrt(2) / 4, -math.sqrt(2) / 2)]
    eps = 1e-14
    for got, want in zip(gcm.cartToSphr(c), s):
        assert sphrNear(got, want, eps)


def test_transpose3():
    a = ((1, 2, 3), (4, 5, 6), (7, 8, 9))
    got = gcm.transpose3(a)
    want = ((1, 4, 7), (2, 5, 8), (3, 6, 9))
    assert got == want


def dmsToRad(neg, d, m, s):
    r = ((d * 60 + m) / (180. * 60) + s / (180. * 3600)) * math.pi
    if neg == '-':
        return -r
    return r


def test_two():
    # two position test data
    mjd = [56123, 56123.01]
    pos = [
        gcm.Sphr(0, dmsToRad(' ', 89, 59, 40)),
        gcm.Sphr(0, dmsToRad(' ', 90, 0, 0))]

    # gcm fit
    g = gcm.GCM(mjd, pos)

    # test the method that returns just residuals
    assert g.res() == [[0.0, 0.0], [0.0, 0.0]]

    # test the method that returns just rms
    assert g.rms() == 0

    # test the method that returns just both
    assert g.rmsRes() == (0, [[0.0, 0.0], [0.0, 0.0]])


def test_three():
    '''Three observations, the middle one an arc second off in both RA
     and Dec.'''
    mjd = [56123, 56123.01, 56123.02]
    pos = [
        gcm.Sphr(dmsToRad(' ', 0, 0, -15), 0),
        gcm.Sphr(dmsToRad(' ', 0, 0, 1), dmsToRad(' ', 0, 0, 1)),
        gcm.Sphr(dmsToRad(' ', 0, 0, 15), 0)]
    eps = 1e-6

    # gcm fit
    g = gcm.GCM(mjd, pos)

    # test the method that returns just residuals
    res = g.res()
    assert near(res[0][0], -1. / 3, eps)
    assert near(res[0][1], -1. / 3, eps)
    assert near(res[1][0], 2. / 3, eps)
    assert near(res[1][1], 2. / 3, eps)
    assert near(res[2][0], -1. / 3, eps)
    assert near(res[2][1], -1. / 3, eps)

    # test the method that returns just rms
    assert near(g.rms(), 2 / 3, eps)

    # test the method that returns just both
    rms, res = g.rmsRes()
    assert near(res[0][0], -1. / 3, eps)
    assert near(res[0][1], -1. / 3, eps)
    assert near(res[1][0], 2. / 3, eps)
    assert near(res[1][1], 2. / 3, eps)
    assert near(res[2][0], -1. / 3, eps)
    assert near(res[2][1], -1. / 3, eps)
    assert near(rms, 2 / 3, eps)


def test_sphrStr():
    assert str(gcm.Sphr(
        dmsToRad(' ', 1, 20, 0),
        dmsToRad('-', 14, 17, 8.5))) == '( 1 20 0.0, -14 17 8.5)'


def test_exceptions():
    mjd = [56123, 56123.01, 56123.02]
    pos = [
        gcm.Sphr(dmsToRad(' ', 0, 0, -15), 0),
        gcm.Sphr(dmsToRad(' ', 0, 0, 1), dmsToRad(' ', 0, 0, 1)),
        gcm.Sphr(dmsToRad(' ', 0, 0, 15), 0)]

    with pytest.raises(ValueError) as e:
        gcm.GCM(mjd[:-1], pos)
    assert str(e.value) == 'date and pos different lengths'

    with pytest.raises(ValueError) as e:
        gcm.GCM(mjd[:1], pos[:1])
    assert str(e.value) == 'at least two positions needed'

    with pytest.raises(ValueError) as e:
        gcm.GCM([mjd[0], mjd[0]], pos[:2])
    assert str(e.value) == 'positive elapsed time needed'

    with pytest.raises(ValueError) as e:
        gcm.GCM(mjd[:2], [pos[0], pos[0]])
    assert str(e.value) == 'motion across sky needed'
