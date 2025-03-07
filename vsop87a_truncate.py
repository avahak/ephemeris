"""
Rounds and discards coefficients in the raw VSOP87 json file and stores 
the result in another json file. 

NOTE Precision of position for the Earth and the Moon based on difference 
of EARTH-MOON and EARTH is low compared to alternative methods of computing them 
(like ELP/MPP02), so we ignore EARTH and do not write it in the truncated json file.
EARTH-MOON is still written.
"""
from dataclasses import dataclass
import numpy as np
import json

import tools.misc as misc
from testing import JPL_DE_ROUTES, compare_pos_vel_functions, test_planet_ephemeris_against_jpl_de
from tools.round_compact import FormattedFloat, FormattedFloatEncoder, round_compact
from vsop87a_ephemeris import VSOP87AEphemeris

@dataclass
class ErrorConfig:
    """
    t_max represents the assumption for maximum |t| when simplifying the coefficients.
    max_error_per_char is maximum allowed error per truncated character.
    """
    t_max: float
    max_error_per_char: float

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
# t_max and max_error_per_char for each planet
ERROR_CONFIGS = {
    'MERCURY': ErrorConfig(50.0, 0.5),
    'VENUS': ErrorConfig(50.0, 0.2),
    'EARTH-MOON': ErrorConfig(50.0, 0.25),
    'MARS': ErrorConfig(30.0, 0.75),
    'JUPITER': ErrorConfig(30.0, 1.0),
    'SATURN': ErrorConfig(20.0, 2.0),
    'URANUS': ErrorConfig(30.0, 2.5),
    'NEPTUNE': ErrorConfig(50.0, 1.0),
}
# Key refers to file, value is multiplier for each max_error_per_char
MAX_ERROR_MULTIPLIER = {
    'small': 5.0e-7,
    'medium': 1.0e-8,
    'large': 3.0e-10
}

INPUT_RAW_JSON_PATH = R'./json/vsop87a_raw.json'
OUTPUT_JSON_PATH = lambda size: f'./json/vsop87a_truncated_{size}.json'
JPL_DE_EPHEMERIS_PATH = R'd:/resources/astro/de/de441.bsp'

def simplify(a: float, b: float, c: float, alpha: int, body_name: str, size: str):
    """
    Truncate digits in coefficients where the estimated error per dropped digit 
    is below max_error_per_char. After truncation, drop the whole term if the 
    total error per character is below max_error_per_char.
    """
    t_max = ERROR_CONFIGS[body_name].t_max
    max_error_per_char = MAX_ERROR_MULTIPLIER[size] * ERROR_CONFIGS[body_name].max_error_per_char * MEAN_DISTANCE_FROM_SUN[body_name]
    tma = np.power(t_max, alpha)

    # First truncate coefficients
    coeffs = [a, b, c]
    error_multipliers = [tma, abs(a)*tma, abs(a)*tma*t_max]
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
    drop_chars_saved = 3 + sum([len(x) for x in coeffs_truncated])  # 3 commas
    drop_error_added = tma*abs(a)
    # NOTE Here we could add a fudge factor to the right hand side, like 0.5*max_error_per_char,
    # to try to balance the two different kinds of truncation.
    if drop_error_added / drop_chars_saved < max_error_per_char:
        # The term is dropped
        return [FormattedFloat.ZERO] * 3
    
    return coeffs_truncated
   

def truncate_series(obj_raw, size): 
    # Considers each term in the series and drops it or simplifies the coefficients
    bodies_truncated = {}
    term_count = {}
    for body_name, groups in obj_raw['bodies'].items():
        if body_name == 'EARTH':     # Skip EARTH for being too inaccurate
            continue

        for group in groups:
            alpha = group['alpha']
            coord = group['coord']
            coeffs = np.array(group['coeffs']).reshape(-1, 3)
            coeffs_truncated = []
            for a, b, c in coeffs:
                ap, bp, cp = simplify(a, b, c, alpha, body_name, size)
                if float(ap) == 0.0:
                    continue

                coeffs_truncated.extend([ap, bp, cp])
                term_count[body_name] = term_count.get(body_name, 0) + 1

            if coeffs_truncated:
                groups_truncated = bodies_truncated.setdefault(body_name, [])
                groups_truncated.append({ 'coord': coord, 'alpha': alpha, 'coeffs': coeffs_truncated })

    print(f'total # of terms = {sum(term_count.values())},\n{term_count = }')

    return { 
        '_comment': f'VSOP87A series truncated, {MAX_ERROR_MULTIPLIER[size]=}',
        # 'lambda': obj_raw['lambda'], 
        'matrix': obj_raw['matrix'], 
        'bodies': bodies_truncated
    }

def write_truncated_json(obj_truncated, size):  
    # Writes the truncated json file
    file_path = OUTPUT_JSON_PATH(size)
    with open(file_path, 'w') as f:
        json_str = json.dumps(obj_truncated, cls=FormattedFloatEncoder, indent=None, separators=(',', ':'))
        json_str = FormattedFloatEncoder.apply_formatting(json_str)
        f.write(json_str)
    print(f'Wrote file: {file_path} ({round(np.ceil(len(json_str)/1024))} kb).')
    return json.loads(json_str)

def show_truncation_error(obj_truncated, obj_raw): 
    # Show error between truncated series and the raw series, and error between
    # the truncated series and JPL DE.
    test_num = 500
    with misc.jplephem_pos_vel(JPL_DE_EPHEMERIS_PATH) as jpl_pos_vel:
        vsop87_raw = VSOP87AEphemeris(obj_raw)
        vsop87_truncated = VSOP87AEphemeris(obj_truncated)
        for body_name in JPL_DE_ROUTES:
            if body_name not in obj_raw['bodies'].keys():
                continue
            pos_vel_raw = lambda t: vsop87_raw.get_pos_vel(body_name, t)
            pos_vel_truncated = lambda t: vsop87_truncated.get_pos_vel(body_name, t)
            print('-'*20, f'{body_name}: truncation errors against raw series', '-'*20)
            compare_pos_vel_functions(pos_vel_truncated, pos_vel_raw, test_num)

        print('-'*20, 'Errors against JPL DE441', '-'*20)
        test_planet_ephemeris_against_jpl_de(vsop87_truncated.get_pos_vel, jpl_pos_vel, test_num)   

def truncate_and_write_file(size):
    print(f'{size = }')
    print('-'*20, 'Truncating series, writing truncated json', '-'*20)
    obj_raw = misc.load_json(INPUT_RAW_JSON_PATH)
    obj_truncated = truncate_series(obj_raw, size)
    obj_truncated = write_truncated_json(obj_truncated, size)

    show_truncation_error(obj_truncated, obj_raw)

if __name__ == '__main__':
    truncate_and_write_file('small')
    truncate_and_write_file('medium')
    truncate_and_write_file('large')