"""
Truncates ELP/MPP02 series by dropping small terms and rounds the coefficients
to take less space in json file.
"""

import numpy as np
import json
import re

import tools
from mpp02_ephemeris import MPP02Ephemeris
from testing import compare_pos_vel_functions

T_MAX = 30.0        # Assumption for max|t| to use in truncation (default 30.0)
THRESHOLD_EXP = 7
THRESHOLD = np.power(10.0, -THRESHOLD_EXP)  # maximum threshold for each modification
# (T_MAX,THRESHOLD_EXP)=(30,7): ~1000 terms, good accuracy (often error<1")

def simplify(c: list[float], alpha: int, limit: float):
    """
    Question to ask here is: how much can we simplify a component (like c[3]) while keeping
    the resulting change in the term below given limit. Each component is changed this way
    and the simplified coefficients are returned.

    Assumption for t is |t| <= T_MAX.
    Term: c0 * t^alpha * sin(c1 + c2*t + c3*t^2 + c4*t^3 + c5*t^4).
    Maximum leeway to change the components:
    c0: limit * T_MAX^(-alpha)
    c1: limit * T_MAX^(-alpha) / |c0|
    c2: limit * T_MAX^(-1-alpha) / |c0|
    ...
    c5: limit * T_MAX^(-4-alpha) / |c0|
    """
    if c[0] == 0.0:
        return [0, 0, 0, 0, 0, 0]
    leeway0 = limit * np.power(T_MAX, -alpha)
    c_new = [tools.simplest_float_in(c[0], c[0]-leeway0, c[0]+leeway0)]
    for k in range(5):
        leeway = leeway0 * np.power(T_MAX, -k) / abs(c[0])
        c_new.append(tools.simplest_float_in(c[k+1], c[k+1]-leeway, c[k+1]+leeway))
    return [tools.convert_integer_floats(c) for c in c_new]

def count_coeffs(obj, report_all=False):
    counts = {}
    for group in obj['groups']:
        key = f'({group['coord']},{group['alpha']})'
        counts[key] = len(group['coeffs']) // 6
    total_count = sum(counts.values())
    return counts if report_all else total_count

def truncate_series(obj_raw): 
    groups_truncated = []
    for group in obj_raw['groups']:
        coord = group['coord']
        alpha = group['alpha']
        coeffs = np.array(group['coeffs']).reshape(-1, 6)
        coeffs_truncated = []
        for c in coeffs:
            # For coords 0,1 (lon,lat) the unit is arcsec and for 2 (r) it is km.
            # We normalize here by reducing to unit sphere.
            limit = THRESHOLD/tools.ARCSEC if coord != 2 else THRESHOLD*384399.0
            c_new = simplify(c, alpha, limit)
            if c_new[0] == 0:
                continue
            coeffs_truncated.extend(c_new)
        if coeffs_truncated:
            groups_truncated.append({ 'coord': coord, 'alpha': alpha, 'coeffs': coeffs_truncated })
    comment = obj_raw['_comment'] + f' :: Truncated ({T_MAX=}, {THRESHOLD=})'
    obj = { '_comment': comment, 'W': obj_raw['W'], 'groups': groups_truncated }
    print(f'Truncated from {count_coeffs(obj_raw)} ({count_coeffs(obj_raw, True)})\n to {count_coeffs(obj)} ({count_coeffs(obj, True)})')
    return obj

def write_data(obj, file_path: str):
    # Write the truncated json file
    with open(file_path, 'w') as f:
        json_string = json.dumps(obj, indent=None, separators=(',', ':'))
        json_string = re.sub(r'(\d+\.?\d*)e-0(\d+)', r'\1e-\2', json_string)
        f.write(json_string)
    print(f'Wrote file: {file_path} ({round(np.ceil(len(json_string)/1024))} kb).')

def truncate_and_show_error():
    test_num = 200
    with tools.jplephem_pos_vel(R'd:/resources/astro/de/de441.bsp') as jpl_pos_vel:
        json_raw_path = './json/mpp02_llr_raw.json'
        obj_raw = tools.load_json(json_raw_path)
        obj_truncated = truncate_series(obj_raw)
        output_json_path = f'./json/mpp02_llr_truncated_{THRESHOLD_EXP}.json'
        write_data(obj_truncated, output_json_path)

        # Show errors in truncation
        mpp02_raw = MPP02Ephemeris(obj_raw)
        mpp02_truncated = MPP02Ephemeris(obj_truncated)
        jpl_moon_pos_vel = lambda t: jpl_pos_vel([399, 3, 301], t)
        print('-'*20, 'ERRORS AGAINST RAW SERIES', '-'*20)
        compare_pos_vel_functions(mpp02_truncated.get_pos_vel, mpp02_raw.get_pos_vel, test_num)
        print('-'*20, 'ERRORS AGAINST JPL DE 441', '-'*20)
        compare_pos_vel_functions(mpp02_truncated.get_pos_vel, jpl_moon_pos_vel, test_num)

if __name__ == '__main__':
    truncate_and_show_error()
