"""
Truncates ELP/MPP02 series by rounding coefficients and dropping small terms 
to take less space in json file.
"""

import numpy as np
import json

import tools.misc as misc
from mpp02_ephemeris import MPP02Ephemeris
from tools.round_compact import FormattedFloat, FormattedFloatEncoder, round_compact
from testing import compare_pos_vel_functions

# Represents the assumption for maximum |t| when simplifying the coefficients.
T_MAX = 10.0
# Key refers to file, value is maximum allowed error per character
MAX_ERROR_PER_CHAR = {
    'small': 5.0e-8,
    'medium': 2.0e-9,
    'large': 4.0e-11
}

# The series gives lon (coord 0) and lat (coord 1) in arcsec, and r (coord 2) in km.
# We normalize by scaling to unit sphere. Since inclination of the Moon is small,
# latitude stays close to 0. Therefore, after normalization the coordinate
# basis vectors of (lon,lat,r) are nearly orthonormal. This means that the relative 
# errors in Cartesian coordinates can be approximated by the errors in 
# these normalized (lon,lat,r) coordinates.
COORD_NORMALIZATION = [misc.ARCSEC, misc.ARCSEC, 1.0 / 385000.0]

INPUT_RAW_JSON_PATH = R'./json/mpp02_llr_raw.json'
OUTPUT_JSON_PATH = lambda size: f'./json/mpp02_llr_truncated_{size}.json'
JPL_DE_EPHEMERIS_PATH = R'd:/resources/astro/de/de441.bsp'

def simplify(coeffs: list[float], alpha: int, coord: int, size: str):
    """
    Truncate digits in coefficients where the estimated error per dropped digit 
    is below max_error_per_char. After truncation, drop the whole term if the 
    total error per character is below max_error_per_char.
    """
    max_error_per_char = MAX_ERROR_PER_CHAR[size]
    ntma = COORD_NORMALIZATION[coord] * np.power(T_MAX, alpha)
    a = coeffs[0]   # amplitude

    # First truncate coefficients
    error_multipliers = [ntma, *[abs(a)*ntma*np.power(T_MAX, k) for k in range(5)]]
    coeffs_truncated = []
    for k, x in enumerate(coeffs):
        x_current = FormattedFloat(str(x))
        for sig_figs in range(17, -1, -1):
            x_proposed = round_compact(x, sig_figs)
            chars_saved = len(x_current) - len(x_proposed)
            error_added = error_multipliers[k] * (abs(x-float(x_proposed)) - abs(x-float(x_current)))
            if (chars_saved > 0) and (error_added / chars_saved < max_error_per_char):
                x_current = x_proposed
        coeffs_truncated.append(x_current)

    # Then consider dropping the term
    drop_chars_saved = 6 + sum([len(x) for x in coeffs_truncated])
    drop_error_added = ntma*abs(a)
    # NOTE Here we could add a fudge factor to the right hand side, like 0.5*max_error_per_char,
    # to try to balance the two different kinds of truncation.
    if drop_error_added / drop_chars_saved < max_error_per_char:
        # The term is dropped
        return [FormattedFloat.ZERO] * 6
    
    return coeffs_truncated

def count_terms(obj, report_all=False):
    counts = {}
    for group in obj['groups']:
        key = f'({group['coord']},{group['alpha']})'
        counts[key] = len(group['coeffs']) // 6
    total_count = sum(counts.values())
    return counts if report_all else total_count

def truncate_series(obj_raw, size): 
    groups_truncated = []
    for group in obj_raw['groups']:
        coord = group['coord']
        alpha = group['alpha']
        coeffs = np.array(group['coeffs']).reshape(-1, 6)
        coeffs_truncated = []
        for c in coeffs:
            c_new = simplify(c, alpha, coord, size)
            if float(c_new[0]) == 0.0:
                continue
            coeffs_truncated.extend(c_new)
        if coeffs_truncated:
            groups_truncated.append({ 'coord': coord, 'alpha': alpha, 'coeffs': coeffs_truncated })
    comment = obj_raw['_comment'] + f', truncated with {T_MAX=}, {MAX_ERROR_PER_CHAR[size]=}'
    obj = { 
        '_comment': comment, 
        'W': obj_raw['W'], 
        'PC': obj_raw['PC'], 
        'QC': obj_raw['QC'], 
        'groups': groups_truncated 
    }
    print(f'Truncated from {count_terms(obj_raw)} ({count_terms(obj_raw, True)})\n to {count_terms(obj)} ({count_terms(obj, True)})')
    return obj

def write_data(obj, size):
    # Write the truncated json file
    file_path = OUTPUT_JSON_PATH(size)
    with open(file_path, 'w') as f:
        json_str = json.dumps(obj, cls=FormattedFloatEncoder, indent=None, separators=(',', ':'))
        json_str = FormattedFloatEncoder.apply_formatting(json_str)
        f.write(json_str)
    print(f'Wrote file: {file_path} ({round(np.ceil(len(json_str)/1024))} kb).')
    return json.loads(json_str)

def show_truncation_error(obj_raw, obj_truncated): 
    # Show errors in truncation
    test_num = 2000
    with misc.jplephem_pos_vel(JPL_DE_EPHEMERIS_PATH) as jpl_pos_vel:
        mpp02_raw = MPP02Ephemeris(obj_raw)
        mpp02_truncated = MPP02Ephemeris(obj_truncated)
        jpl_moon_pos_vel = lambda t: jpl_pos_vel([399, 3, 301], t)
        print('-'*20, 'ERRORS AGAINST RAW SERIES', '-'*20)
        compare_pos_vel_functions(mpp02_truncated.get_pos_vel, mpp02_raw.get_pos_vel, test_num)
        print('-'*20, 'ERRORS AGAINST JPL DE 441', '-'*20)
        compare_pos_vel_functions(mpp02_truncated.get_pos_vel, jpl_moon_pos_vel, test_num)

def truncate_and_write_file(size):
    print(f'{size = }')
    obj_raw = misc.load_json(INPUT_RAW_JSON_PATH)
    obj_truncated = truncate_series(obj_raw, size)
    obj_truncated = write_data(obj_truncated, size)
    
    show_truncation_error(obj_raw, obj_truncated)

if __name__ == '__main__':
    truncate_and_write_file('small')
    truncate_and_write_file('medium')
    truncate_and_write_file('large')