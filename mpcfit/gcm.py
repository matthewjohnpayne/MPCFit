# /mpcfit/mpcfit/gcm.py
"""
# --------------------------------------------------------------
# Oct 2018
# Payne
#
# Fit great circle motion to a list of sky positions.
# The fit is not only to a great circle ...
# ... but also to linear motion along the great circle.
#
# Straight port from work by Sonia Keys
# - Some of Sonia's notes are in /mpcfit/mpcfit/gcm.md
# - At present, this is a striaght C-port, so NOT pythonic
#   --- Could use a lot more numpy & get a lot of speed-up
# ...
#
# --------------------------------------------------------------
"""


# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import collections
import math

# Import neighboring packages
# --------------------------------------------------------------
#import mpcutilities.phys_const as PHYS




# GCM functions & classes
# --------------------------------------------------------------

math.tau = math.pi * 2.0

def R2D(r):
    '''Radians to degrees.'''
    return r * 180 / math.pi


Cart = collections.namedtuple('Cart', ['x', 'y', 'z'])
'''3D X-Y-Z coordinates.'''

Sphr = collections.namedtuple('Sphr', ['ra', 'dec'])
'''RA and dec in radians.'''


def dms(r):
    '''Utility formatter for Sphr.__str__.'''
    neg = ' '
    if r < 0:
        neg = '-'
        r = -r
    s = R2D(r) * 3600
    m = int(s) // 60
    s -= m * 60
    d = m // 60
    m -= d * 60
    return '{}{} {} {}'.format(neg, d, m, s)

Sphr.__str__ = lambda self: '({}, {})'.format(dms(self.ra), dms(self.dec))
'''Sexagesimal stringer useful for debugging.'''




def cartToSphr(cl):
    '''Convert 3D cartesian coordinates to spherical (RA, Dec) in radians.

    Parameters
    ----------
    cl : iterable of objects with x, y, z attributes

    Returns
    -------
    list of Sphr
    '''
    return [Sphr(math.fmod(math.atan2(c.y, c.x) + math.tau, math.tau),
        math.asin(c.z)) for c in cl]


def sphrToCart(sl):
    '''Convert spherical coordinates to 3D cartesian unit vectors.

    Parameters
    ----------
    sl : iterable of objects with ra, dec attributes

    Returns
    -------
    list of Cart

    '''
    return [Cart(t * math.cos(s.ra), t * math.sin(s.ra), math.sin(s.dec))
        for s, t in zip(sl, [math.cos(s.dec) for s in sl])]


def cartCross(a, b):
    '''3D vector cross product, a x b.

    Parameters
    ----------
    a, b : objects with x, y, z attributes
    '''
    return Cart(*np.cross(a,b).tolist())   # MJP: Forcing into a Cart for now ...
    """
    return Cart(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x)
    """

def cartDot(a, b):
    '''3D vector dot product.
        
        Parameters
        ----------
        a, b : objects with x, y, z attributes
        '''
    return np.inner(a,b)
    #return a.x * b.x + a.y * b.y + a.z * b.z

def cartSquare(a):
    '''a dot a

    Parameters
    ----------
    a : object with x, y, z attributes
    '''
    return cartDot(a, a)




def cartRot(rm, aa):
    '''Rotate cartesian coordinates.

    It does matrix multiplication, specialized to broadcast multiplication of
    a rotation matrix to a list of carts.'''
    return [Cart(
        rm[0][0] * a.x + rm[0][1] * a.y + rm[0][2] * a.z,
        rm[1][0] * a.x + rm[1][1] * a.y + rm[1][2] * a.z,
        rm[2][0] * a.x + rm[2][1] * a.y + rm[2][2] * a.z) for a in aa]


def transpose3(a):
    '''Transpose of a 3x3 matrix.

    Parameters
    ----------
    a : 3x3 indexable

    Returns
    -------
    3x3 tuple of tuples
    '''
    return (
        (a[0][0], a[1][0], a[2][0]),
        (a[0][1], a[1][1], a[2][1]),
        (a[0][2], a[1][2], a[2][2]))


eye3 = ((1., 0., 0.), (0., 1., 0.), (0., 0., 1.))
'''3x3 identity matrix as tuple of tuples of floats'''


class GCM:
    def __init__(self, date, pos):
        '''Fit positions and times to linear motion along a great circle.

        Parameters
        ----------
        date : iterable of floats
            Iterable parallel to pos, representing times at positions.
            Times must be on some linear time scale, for example JD or MJD.
        pos : iterable of objects with 'ra' and 'dec' attributes.
            RA and dec values must be in radians.
        '''
        if len(date) != len(pos):
            raise ValueError('date and pos different lengths')
        if len(date) < 2:
            raise ValueError('at least two positions needed')
        if date[0] >= date[-1]:
            raise ValueError('positive elapsed time needed')
        if pos[0] == pos[-1]:
            raise ValueError('motion across sky needed')

        # convert obs to cartesian
        c = sphrToCart(pos)

        # vector normal to motion
        norm = cartCross(c[0], c[-1])
        nmag2 = cartSquare(norm)
        nmag = math.sqrt(nmag2)

        # rotation angle is angle from norm to +z
        xy = math.hypot(norm.x, norm.y)

        try:
            tana = xy / norm.z
        except ZeroDivisionError:
            tana = math.inf
        if not (abs(tana) > .0003):
            # if norm is close to the pole already, don't mess with rotation.
            self.rm = eye3
            rs = list(pos)  # copy
        else:
            sina = xy / nmag
            cosa = norm.z / nmag
            lmix = norm.y / xy
            lmiy = -norm.x / xy

            '''build rm, rotation matrix that will rotate the coordinate system
            to the obs, so that a cylindrical projection will have negligible
            distortion.'''

            sinagx = sina * lmix
            sinagy = sina * lmiy
            onemcosa = 1 - cosa
            onemcosagx = onemcosa * lmix
            onemcosagxgy = onemcosagx * lmiy
            self.rm = (
                (cosa + onemcosagx * lmix, onemcosagxgy, sinagy),
                (onemcosagxgy, cosa + onemcosa * lmiy * lmiy, -sinagx),
                (-sinagy, sinagx, cosa))

            # rotate all of cart
            rotated = cartRot(self.rm, c)

            # transpose rotation array so it will derotate, after
            # least-squares fit.
            transpose3(self.rm)

            # convert back to spherical coordinates for adjustment.
            # this does a cylindrical projection.
            rs = cartToSphr(rotated)

        # normalize ra to near 0 to avoid wraparound problems
        self.ra0 = ra0 = rs[0].ra
        for i, s in enumerate(rs):
            rs[i] = Sphr(
                math.fmod(s.ra + 3 * math.pi - ra0, 2 * math.pi) - math.pi,
                rs[i].dec)
        self.rs = rs

        # normalize time to near 0 to maintain precision
        self.t0 = t0 = date[0]
        nt = [t - t0 for t in date]
        self.nt = nt

        if len(nt) == 2:
            # least squares fit not needed with just two points.
            self.r0 = 0.
            self.rr = rs[1].ra / nt[1]
            self.d0 = 0.
            self.dr = 0.
        else:
            # here's the least squares stuff
            sumt = 0.
            sumra = 0.
            sumdec = 0.
            sumt2 = 0.
            sumtra = 0.
            sumtdec = 0.
            for i in range(len(nt)):
                sumt += nt[i]
                sumra += rs[i].ra
                sumdec += rs[i].dec
                sumt2 += nt[i] * nt[i]
                sumtra += nt[i] * rs[i].ra
                sumtdec += nt[i] * rs[i].dec
            invd = 1 / (len(nt) * sumt2 - sumt * sumt)
            self.r0 = invd * (sumra * sumt2 - sumtra * sumt)
            self.rr = invd * (len(nt) * sumtra - sumra * sumt)
            self.d0 = invd * (sumdec * sumt2 - sumtdec * sumt)
            self.dr = invd * (len(nt) * sumtdec - sumdec * sumt)

    def rms(self):
        '''RMS of residuals, in arc seconds.'''
        rms, _ = self.rmsRes()
        return rms

    def rmsRes(self):
        '''RMS and residuals, in arc seconds.'''
        res = self.res()
        return math.sqrt(sum(r[0] * r[0] + r[1] * r[1] for r in res) /
            len(res)), res

    def res(self):
        '''Residuals in arc seconds.'''

        rsc = [Sphr(self.r0 + self.rr * t, self.d0 + self.dr * t)
            for t in self.nt]

        ''' residuals are computed on rotated-and-derotated observed to
        minimize any systematic errors from rotating and derotating the
        computed.'''

        # fix up ra on both observed and computed
        rso = [Sphr(o.ra + self.ra0, o.dec) for o in self.rs]
        rsc = [Sphr(c.ra + self.ra0, c.dec) for c in rsc]

        # rotate both back up to original place in the sky
        rco = sphrToCart(rso)
        rcc = sphrToCart(rsc)
        co = cartRot(self.rm, rco)
        cc = cartRot(self.rm, rcc)
        so = cartToSphr(co)
        sc = cartToSphr(cc)

        return [[  # residuals in arc seconds
            R2D(o.dec - c.dec) * 3600,
            R2D(o.ra - c.ra) * 3600 * math.cos(c.dec)]
            for o, c in zip(so, sc)]


class GCMNEW:
    def __init__(self, date, RA, DEC, pos):
        '''Fit positions and times to linear motion along a great circle.
            
        Parameters
        ----------
        date : iterable of floats
            Iterable parallel to RA & DEC, representing times at positions.
            Times must be on some linear time scale, for example JD or MJD.
        RA : iterable of objects
            Iterable parallel to date & DEC, representing right-ascention
            RA values must be in radians.
        DEC: iterable of objects
            Iterable parallel to date & RA, representing declination
        '''
        if len(date) != len(RA):
            raise ValueError('date and RA different lengths')
        if len(date) != len(DEC):
            raise ValueError('date and DEC different lengths')
        if len(date) < 2:
            raise ValueError('at least two positions needed')
        if date[0] >= date[-1]:
            raise ValueError('positive elapsed time needed')
        if RA[0] == RA[-1] and DEC[0] == DEC[-1] :
            raise ValueError('motion across sky needed')
        
        # convert RA,DEC observations to cartesian unit vectors
        self.cartUV = np.reshape(
                                 np.array( [np.cos(DEC) * np.cos(RA), np.cos(DEC) * np.sin(RA), np.sin(DEC)]) ,
                                 (-1,3)
                                 )

        # vector normal to motion
        # *** MJP : Note that at this point Sonia uses first and last points in the tracklet ( [0] & [-1] )
        # *** MJP : May well be OK to do this, but I am slightly suspicious: why not fit generally ?
        c = sphrToCart(pos)
        print("pos=", pos)
        print("c=", c)

        # rotation angle is angle from norm to +z

        # if norm is close to the pole already, don't mess with rotation.

        # Build "rm", the rotation matrix that will rotate the coordinate system
        # to the obs, so that a cylindrical projection will have negligible
        # distortion.'''

        # rotate all of cart

        # transpose rotation array so it will derotate,
        # after least-squares fit.

        # convert back to spherical coordinates for adjustment.
        # this does a cylindrical projection.

        # normalize ra to near 0 to avoid wraparound problems

        # normalize time to near 0 to maintain precision


        # Calculate least squares

        # Least squares fit not needed with just two points.


