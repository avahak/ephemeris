# Miscellaneous helper functions

import contextlib
import numpy as np
import json
import time
import functools
from jplephem.spk import SPK

DEG = np.pi / 180.0
ARCSEC = np.pi / 180.0 / 3600

# def _json_structure(data, indent, max_iterable_entries=4):
#     # Helper function for json_structure_string
#     lines = []
#     pad = ' ' * 4
#     prefix = pad * indent

#     if isinstance(data, dict):
#         lines.append(f'{prefix}- dict of length {len(data)}')
#         for count, (key, value) in enumerate(data.items()):
#             if (count >= max_iterable_entries-1) and (count < len(data)-1):
#                 if count == max_iterable_entries-1:
#                     lines.append(f'{prefix}- ...(items {max_iterable_entries-1}-{len(data)-2})')
#                 continue
#             sublines = _json_structure(value, indent + 1, max_iterable_entries)
#             if isinstance(sublines, str):
#                 lines.append(f'{prefix}- {key} {sublines}')
#             else:
#                 lines.append(f'{prefix}- {key}')
#                 lines.extend(sublines)

#     elif isinstance(data, list) and data:
#         if not data:
#             return 'list (empty)'
        
#         is_homogenous = True
#         first = json_structure_string(data[0])
#         for entry in data[1:]:
#             mask = json_structure_string(entry)
#             if mask != first:
#                 is_homogenous = False
#                 break
#         if is_homogenous:
#             sublines = _json_structure(data[0], indent + 1, max_iterable_entries)
#             s = f'{prefix}- list (homogenous of length {len(data)}) of form'
#             if isinstance(sublines, str):
#                 lines.append(f'{s} {sublines}')
#             else:
#                 lines.append(f'{s}')
#                 lines.extend(sublines)
#         else:
#             lines.append(f'{prefix}- list (length {len(data)})')
#             for index, entry in enumerate(data):
#                 if (index >= max_iterable_entries-1) and (index < len(data)-1):
#                     if index == max_iterable_entries-1:
#                         lines.append(f'{prefix}- ...(items {max_iterable_entries-1}-{len(data)-2})')
#                     continue
#                 s = f'{prefix}- item {index}'
#                 sublines = _json_structure(entry, indent + 1, max_iterable_entries)
#                 if isinstance(sublines, str):
#                     lines.append(f'{s} {sublines}')
#                 else:
#                     lines.append(f'{s}')
#                     lines.extend(sublines)

#     elif isinstance(data, str):
#         return 'str'
#     elif isinstance(data, int):
#         return 'int'
#     elif isinstance(data, float):
#         return 'float'
#     elif isinstance(data, bool):
#         return 'bool'
#     elif data is None:
#         return 'null'
#     else:
#         return 'unknown'
#     return lines

# def json_structure_string(data):
#     """
#     Writes out the general structure of a (json) object. 
#     Scuffy, for entertainment purposes only!
#     """
#     structure = _json_structure(data, 0, 4)
#     if isinstance(structure, list):
#         return '\n'.join(structure)
#     return str(structure)

def format_time(t: float) -> str:
    if t >= 3600.0:
        hours = int(t / 3600.0)
        minutes = int((t % 3600.0) / 60.0)
        return f'{hours}h {minutes}min'
    elif t >= 60.0:
        minutes = int(t / 60.0)
        seconds = int(t % 60.0)
        return f'{minutes}min {seconds}s'
    return f'{t:.2g}s'

def time_it(f):
    """
    Writes elapsed time for function execution.
    """
    @functools.wraps(f)
    def wrapper(*pos_args, **keyw_args):
        time0 = time.perf_counter()
        return_value = f(*pos_args, **keyw_args)
        print(f'Function {f.__name__} took {format_time(time.perf_counter()-time0)}.')
        return return_value
    return wrapper

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
        'std': v.std(),
        'min': v.min(),
        'max': v.max()
    }

def dms_string(x: float, arcsec_digits=0) -> str:
    if x < 0.0:
        return f'-{dms_string(-x, arcsec_digits)}'
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
    arcsec_str = f'{arcseconds:.{arcsec_digits}f}\"'
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