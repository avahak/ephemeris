# Miscellaneous helper functions

import contextlib
import numpy as np
import json
from jplephem.spk import SPK

DEG = np.pi / 180.0
ARCSEC = np.pi / 180.0 / 3600

def simplest_float_in(x: float, a: float, b: float):
    """
    Finds number x0 in (a,b) with simplest decimal expansion, i.e. least number 
    of significant digits. If there are multiple, picks the one closest to x.
    BUG Unexpected behavior in situations like (x,a,b)=(1,0.5,10000), needs adjustment.
    """
    if (x == 0.0) or (a <= 0.0 and b >= 0.0):
        return 0.0
    # If too much precision is needed, just return x
    if (b-a) / abs(x) < 1.0e-10:
        return x
    if x < 0.0:
        return -simplest_float_in(-x, -b, -a)
    # Now: 0 < a < x < b

    mag = np.floor(np.log10(b - a))
    scale = np.power(10, mag)
    scale10 = 10 * scale
    # Now: scale <= b-a < scale10
    sig = 12
    x1 = round_sig(np.floor(x / scale10) * scale10, sig)
    x2 = round_sig(np.ceil(x / scale10) * scale10, sig)
    x3 = round_sig(np.floor(x / scale) * scale, sig)
    x4 = round_sig(np.ceil(x / scale) * scale, sig)
    if x1 >= a and x1 <= b:
        return x1
    if x2 >= a and x2 <= b:
        return x2
    if x3 >= a and x3 <= b and x4 >= a and x4 <= b:
        return x3 if abs(x-x3) <= abs(x-x4) else x4     # x0 is closer of x3, x4 to x
    if x3 >= a and x3 <= b:
        return x3
    return x4

def cartesian_from_spherical(r, theta, phi):
    """
    Returns Cartesian coordinates given spherical coordinates.
    """
    q = r*np.cos(theta)
    return np.array([q*np.cos(phi), q*np.sin(phi), r*np.sin(theta)])

def spherical_from_cartesian(x, y, z):
    """
    Returns spherical coordinates given Cartesian coordinates.
    """
    r = np.sqrt(x*x + y*y + z*z)
    theta = np.arcsin(z/r)
    phi = np.arctan2(y, x)
    return np.array([r, theta, phi])

def round_sig(x: float, sig_figs: int) -> float:
    """
    Rounds number to given number of significant figures
    """
    if x == 0 or not np.isfinite(x):
        return x
    magnitude = np.floor(np.log10(np.abs(x)))
    scale = np.power(10, sig_figs - magnitude - 1)
    return np.round(x * scale) / scale

def convert_integer_floats(obj):
    """
    Converts floats that are also integers to integer type
    """
    if isinstance(obj, float) and obj.is_integer():
        return int(obj)
    return obj

def load_json(file_name: str):
    with open(file_name, 'r') as f:
        data = json.loads(f.read())
    return data

def stats(v: np.ndarray):
    if v.size == 0:
        raise ValueError('Input array must not be empty')
    return {
        'n': v.size,
        'mean': v.mean(),
        'stdDev': v.std(),
        'min': v.min(),
        'max': v.max()
    }

def dms_string(x: float, sig_figs: int=2) -> str:
    # NOTE significant figures is not really well defined here.
    if x < 0.0:
        return f'-{dms_string(-x, sig_figs)}'
    total_degrees = x * 180.0 / np.pi
    degrees = int(total_degrees)
    remainder = (total_degrees - degrees) * 60 
    arcminutes = int(remainder)
    arcseconds = (remainder - arcminutes) * 60
    if degrees != 0:
        deg_str = f'{degrees}\u00B0'
    else:
        deg_str = ''
    if arcminutes != 0 or degrees != 0:
        arcmin_str = f"{arcminutes}'"
    else:
        arcmin_str = ''
    arcsec_str = f'{round_sig(arcseconds, sig_figs)}\"'
    return ''.join([deg_str, arcmin_str, arcsec_str])
    # return ' '.join(filter(None, [deg_str, arcmin_str, arcsec_str]))

def jpl_segment_getter(kernel):
    """
    Returns function that chooses correct segment in the bsp file.
    This function is only needed for de441.bsp because jplephem segment lookup mechanism
    "kernel[3,399]" does not support multiple segments for one pair, which happens in de441.
    Otherwise you can just use "kernel[center,target]" to get the segment.
    """
    segment_index = {}
    for segment in kernel.segments:
        segment_index.setdefault((segment.center, segment.target), []).append(segment)

    def get_segment(center, target, jd):
        for segment in segment_index[center, target]:
            if segment.start_jd <= jd and jd <= segment.end_jd:
                return segment
        return None
    return get_segment

def jpl_get_pos_vel(segment_getter, route, t):
    """ 
    Here segment_getter can be jplephem SPK object or function obtained from jpl_segment_getter.
    Units: km, julian day.
    """
    jd = 2451545.0 + t*36525.0
    pos = np.zeros(3)
    vel = np.zeros(3)
    for k in range(len(route)-1):
        tr = (route[k], route[k+1])
        segment = segment_getter[min(tr), max(tr)] if type(segment_getter).__name__ == 'SPK' \
                else segment_getter(min(tr), max(tr), jd)
        rpos, rvel = segment.compute_and_differentiate(jd)
        pos += np.sign(tr[1]-tr[0]) * rpos
        vel += np.sign(tr[1]-tr[0]) * rvel
    return pos, vel

@contextlib.contextmanager
def jplephem_pos_vel(jpl_ephemeris_file_path):
    """
    Convenience context manager that on enter provides a function 
    for computing position and velocity.
    """
    with SPK.open(jpl_ephemeris_file_path) as kernel:
        segment_getter = jpl_segment_getter(kernel)
        fn = lambda route, t: jpl_get_pos_vel(segment_getter, route, t)
        yield fn

def rotation_matrix(k, theta):
    # Always CHECK the sign with rotations! Source of sign differences:
    # https://en.wikipedia.org/wiki/Active_and_passive_transformation
    c, s = np.cos(theta), np.sin(theta)
    rot = np.zeros((3, 3))
    k1, k2 = (k+1)%3, (k+2)%3
    rot[k,k] = 1
    rot[k1,k1] = c
    rot[k1,k2] = -s
    rot[k2,k1] = s
    rot[k2,k2] = c
    return rot

def precession(t):
    # from IAUCircular179.pdf, p.44
    ze = (2.650545+t*(2306.083227+t*(0.2988499+t*(0.01801828 \
        +t*(-0.000005971+t*(-0.0000003173))))))*ARCSEC
    z = (-2.650545+t*(2306.077181+t*(1.0927348+t*(0.01826837 \
        +t*(-0.000028596+t*(-0.0000002904))))))*ARCSEC
    th = (t*(2004.191903+t*(-0.4294934+t*(-0.04182264 \
        +t*(-0.000007089+t*(-0.0000001274))))))*ARCSEC
    return ze, z, th

def precession_matrix(t):
    # Precessions matrix only
    ze, z, th = precession(t)
    return rotation_matrix(2, z) @ rotation_matrix(1, -th) @ rotation_matrix(2, ze)

def nutation(t):
    # Source: Fundamental Astronomy, p39
    d = t*36525
    ep = (23.439291111+t*(-0.013004167+t*(-0.000000164+t*(0.000000504))))*DEG
    a1 = (125.0-0.05295*d)*DEG
    a2 = (200.9+1.97129*d)*DEG
    d_ep = (0.0026*np.cos(a1)+0.0002*np.cos(a2))*DEG
    d_psi = (-0.0048*np.sin(a1)-0.0004*np.sin(a2))*DEG
    return ep, d_ep, d_psi

def nutation_matrix(t):
    # Nutation matrix only
    ep, d_ep, d_psi = nutation(t)
    return rotation_matrix(0, ep+d_ep) @ rotation_matrix(2, d_psi) @ rotation_matrix(0, -ep)

def nutation_precession_matrix(t):
    # Precession and nutation combined
    return nutation_matrix(t) @ precession_matrix(t)

if __name__ == '__main__':
    pass