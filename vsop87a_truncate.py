"""
Discards and rounds the coefficients in the raw VSOP87 json file based on a threshold,
and stores the result in a "truncated" json file. 

NOTE Precision of position for the Earth and the Moon based on difference 
of EARTH-MOON and EARTH is low compared to alternative methods of computing them, 
so we ignore EARTH and do not write it in the truncated json file.
"""
import numpy as np
import json
import re

import tools
from testing import JPL_DE_ROUTES, compare_pos_vel_functions, test_planet_ephemeris_against_jpl_de
from vsop87a_ephemeris import VSOP87Ephemeris

MEAN_DISTANCE_FROM_SUN = {
    'MERCURY': 0.39, 
    'VENUS': 0.72, 
    'MARS': 1.52, 
    'JUPITER': 5.2, 
    'SATURN': 9.54, 
    'URANUS': 19.2, 
    'NEPTUNE': 30.06, 
    'EARTH-MOON': 1.0,
}

# Choosing min(mean_dist, 1-mean_dist) means focusing 
# on relative error w.r.t. both the Sun and the Earth.
THRESHOLD_WEIGHT = { key: min(value, abs(1.0-value) if key != 'EARTH-MOON' else 1.0) for key, value in MEAN_DISTANCE_FROM_SUN.items()}
THRESHOLD_EXP = 7
# THRESHOLD_BASE represents the general amount of wiggleroom to round or truncate the coeffs.
# Here 1.0e-7 is a good starting point resulting in relatively small error and file size
THRESHOLD_BASE = np.power(10.0, -THRESHOLD_EXP)
# T_MAX represents the assumption for maximum |t| when simplifying the coefficients.
T_MAX = 5.0

INPUT_RAW_JSON_PATH = R'./json/vsop87a_raw.json'
OUTPUT_JSON_PATH = f'./json/vsop87a_truncated_{THRESHOLD_EXP}.json'
JPL_DE_EPHEMERIS_PATH = R'd:/resources/astro/de/de441.bsp'

def simplify(a: float, b: float, c: float, limit: float, alpha: int):
    """
    Theory here: Assume |t|<=T_MAX. The terms in the series have form 
    t^alpha * a * cos(b + c*t),
    where a, b, c can be simplified. How much can each of the terms a, b, c be modified
    (replacing a with a' with |a-a'|<delta_a etc.) while affecting the term less than the limit? 
    Over all possible |t|<=T_MAX the equations for upper bounds for modifications delta are 
        delta_a <= limit / T_MAX^alpha
        delta_b * |a| <= limit / T_MAX^alpha
        delta_c * |a| <= limit / T_MAX^alpha
    So we can take delta_a=limit/T_MAX^alpha, delta_b=delta_c=THRESHOLD/|a|/T_MAX^alpha.
    This function returns the simplest a',b',c' found in the allowed intervals.
    """
    if a == 0.0:
        return 0, 0, 0
    delta_a = limit / np.power(T_MAX, alpha)
    delta_bc = delta_a / abs(a)
    ap = tools.convert_integer_floats(tools.simplest_float_in(a, a-delta_a, a+delta_a))
    bp = tools.convert_integer_floats(tools.simplest_float_in(b, b-delta_bc, b+delta_bc))
    cp = tools.convert_integer_floats(tools.simplest_float_in(c, c-delta_bc, c+delta_bc))
    return ap, bp, cp

def truncate_series(obj_raw): 
    # Considers each term in the series and drops it or simplifies the coefficients
    # based on threshold (THRESHOLD_WEIGHT[body_name] * THRESHOLD_BASE).
    bodies_truncated = {}
    term_count = {}
    for body_name, groups in obj_raw['bodies'].items():
        if body_name == 'EARTH':     # Skip EARTH for being too inaccurate
            continue
        body_limit = THRESHOLD_WEIGHT[body_name] * THRESHOLD_BASE

        for group in groups:
            alpha = group['alpha']
            coord = group['coord']
            coeffs = np.array(group['coeffs']).reshape(-1, 3)
            coeffs_truncated = []
            for a, b, c in coeffs:
                ap, bp, cp = simplify(a, b, c, body_limit, alpha)
                if ap == 0:
                    continue

                coeffs_truncated.extend([ap, bp, cp])
                term_count[body_name] = term_count.get(body_name, 0) + 1

            if coeffs_truncated:
                groups_truncated = bodies_truncated.setdefault(body_name, [])
                groups_truncated.append({ 'coord': coord, 'alpha': alpha, 'coeffs': coeffs_truncated })

    print(f'total # of terms = {sum(term_count.values())},\n{term_count = }')

    return { 
        '_comment': f'VSOP87A truncated series with {T_MAX=}, {THRESHOLD_EXP=}',
        # 'lambda': obj_raw['lambda'], 
        'matrix': obj_raw['matrix'], 
        'bodies': bodies_truncated
    }

def write_truncated_json(obj_truncated):  
    # Writes the truncated json file
    with open(OUTPUT_JSON_PATH, 'w') as f:
        json_string = json.dumps(obj_truncated, indent=None, separators=(',', ':'))
        # Remove leading zeros in negative exponents of scientific notation floats, HACKY!
        json_string = re.sub(r'(\d+\.?\d*)e-0(\d+)', r'\1e-\2', json_string)
        f.write(json_string)
    print(f'Wrote file: {OUTPUT_JSON_PATH} ({round(np.ceil(len(json_string)/1024))} kb).')

def show_truncation_error(obj_truncated, obj_raw): 
    # Show error between truncated series and the raw series, and error between
    # the truncated series and JPL DE.
    test_num = 500
    with tools.jplephem_pos_vel(JPL_DE_EPHEMERIS_PATH) as jpl_pos_vel:
        vsop87_raw = VSOP87Ephemeris(obj_raw)
        vsop87_truncated = VSOP87Ephemeris(obj_truncated)
        for body_name in JPL_DE_ROUTES:
            if body_name not in obj_raw['bodies'].keys():
                continue
            pos_vel_raw = lambda t: vsop87_raw.get_pos_vel(body_name, t)
            pos_vel_truncated = lambda t: vsop87_truncated.get_pos_vel(body_name, t)
            print('-'*20, f'{body_name}: truncation errors against raw series', '-'*20)
            compare_pos_vel_functions(pos_vel_truncated, pos_vel_raw, test_num)

        print('-'*20, 'Errors against JPL DE441', '-'*20)
        test_planet_ephemeris_against_jpl_de(vsop87_truncated.get_pos_vel, jpl_pos_vel, test_num)   

def truncate_and_test():
    print('-'*20, 'Truncating series, writing truncated json', '-'*20)
    obj_raw = tools.load_json(INPUT_RAW_JSON_PATH)
    obj_truncated = truncate_series(obj_raw)
    write_truncated_json(obj_truncated)

    show_truncation_error(obj_truncated, obj_raw)


if __name__ == '__main__':
    truncate_and_test()