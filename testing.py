"""
Checks VSOP87Ephemeris and MPP02Ephemeris implementations against
reference output from vsop87.chk and ELP/MPP02 Fortran code. Also computes
errors against JPL DE ephemeris.
"""
import numpy as np

import tools
from vsop87a_ephemeris import VSOP87Ephemeris
from mpp02_ephemeris import MPP02Ephemeris

JPL_DE_EPHEMERIS_PATH = R'd:/resources/astro/de/de441.bsp'

# ecliptic to equatorial J2000.0
ROT_ECL_EQU = tools.rotation_matrix(0, 84381.448/3600*np.pi/180.0)

JPL_DE_ROUTES = {
    'MERCURY': [10, 0, 1, 199],
    'VENUS': [10, 0, 2, 299],
    'EARTH-MOON': [10, 0, 3],
    'EMB-EARTH': [3, 399],
    'EMB-MOON': [3, 301],
    'MARS': [10, 0, 4],
    'JUPITER': [10, 0, 5],
    'SATURN': [10, 0, 6],
    'URANUS': [10, 0, 7],
    'NEPTUNE': [10, 0, 8],
    'PLUTO': [10, 0, 9],
}

def generate_test_times(num):
    """
    Returns time intervals (J2000.0-r, J2000.0+r) and a list of num evenly spaced times 
    (Julian centuries since J2000.0) for each.
    """
    times = {}
    interval_radii = [15000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 30]
    for r in interval_radii:
        times_r = []
        for k in range(num):
            s = k / (num-1)
            t = (-r + 2.0*s*r) / 100.0
            times_r.append(t)
        times[f'(-{r},{r})'] = times_r
    return times

def run_vsop87_checks(vsop87: VSOP87Ephemeris):
    """
    Runs checks from vsop87.chk for VSOP87A.
    """
    tests = tools.load_json(R'./json/vsop87a_tests.json')
    errors_pos = []
    errors_vel = []
    for test in tests:
        body = test['body']
        jd = test['jd']
        t = (jd - 2451545.0) / 36525.0
        correct_pos = np.array(test['p'])
        correct_vel = np.array(test['v'])
        pos, vel = vsop87.get_pos_vel(body, t)
        # The matrix is not applied in the test data so unapply it. 
        # The length unit in the test data is AU so convert from km to AU.
        pos = np.linalg.inv(vsop87.matrix) @ pos / VSOP87Ephemeris.AU
        vel = np.linalg.inv(vsop87.matrix) @ vel / VSOP87Ephemeris.AU
        error_pos = np.linalg.norm(pos-correct_pos) / np.linalg.norm(correct_pos)
        error_vel = np.linalg.norm(vel-correct_vel) / np.linalg.norm(correct_vel)
        errors_pos.append(error_pos)
        errors_vel.append(error_vel)
    print('Code port error in pos:', tools.stats(np.array(errors_pos)))
    print('Code port error in vel:', tools.stats(np.array(errors_vel)))

def run_mpp02_tests(mpp02_llr, mpp02_405):
    """
    The output of MPP02Ephemeris should match the sample output of the original Fortran code
    when using nontruncated series. This is tested here.
    """
    FORTRAN_OUTPUT = [
        {'mode': 'LLR', 'jd': 2444239.5, 'p': [43890.282400519, 381188.727452277, -31633.381652398], 'v': [-87516.197484182, 13707.664435129, 2754.221242424]}, 
        {'mode': 'LLR', 'jd': 2446239.5, 'p': [-313664.596449897, 212007.266738547, 33744.751203895], 'v': [-47315.912806896, -75710.875009121, -1475.628688714]}, 
        {'mode': 'LLR', 'jd': 2448239.5, 'p': [-273220.060671398, -296859.768222889, -34604.356996204], 'v': [60542.327590279, -58162.316742506, 2270.886905259]}, 
        {'mode': 'LLR', 'jd': 2450239.5, 'p': [171613.142799329, -318097.337502489, 31293.548240386], 'v': [83266.779903589, 42585.830284948, -1695.826110683]}, 
        {'mode': 'LLR', 'jd': 2452239.5, 'p': [396530.006351246, 47487.922488616, -36085.309034347], 'v': [-12664.286936043, 83512.757191241, 1507.367556647]}, 
        {'mode': '405', 'jd': 2500000.5, 'p': [274034.591031876, 252067.536890015, -18998.755190615], 'v': [-62463.613377925, 65693.963923377, 6595.328896382]}, 
        {'mode': '405', 'jd': 2300000.5, 'p': [353104.31359001, -195254.11808254, 34943.545919834], 'v': [39543.136778978, 74373.180694451, -700.653514055]}, 
        {'mode': '405', 'jd': 2100000.5, 'p': [-19851.276740747, -385646.177171305, -27597.661337993], 'v': [87539.407435142, -7599.684836443, -4960.443599351]}, 
        {'mode': '405', 'jd': 1900000.5, 'p': [-370342.792543067, -37574.255320663, -4527.918405246], 'v': [12255.287457433, -89710.975082649, 7649.44284506]}, 
        {'mode': '405', 'jd': 1700000.5, 'p': [-164673.047197812, 367791.713293163, 31603.98027061], 'v': [-75884.688148488, -35802.265583067, -4239.598946665]}
    ]

    errors_p = []
    errors_v = []
    for test in FORTRAN_OUTPUT:
        mpp02 = mpp02_llr if test['mode'] == 'LLR' else mpp02_405
        t = (test['jd'] - 2451545.0) / 36525.0
        p0, v0 = test['p'], test['v']
        p, v = mpp02.get_pos_vel(t)
        p = ROT_ECL_EQU.T @ p   # back to ecliptic
        v = ROT_ECL_EQU.T @ v
        err_p = np.linalg.norm(p-p0) / np.linalg.norm(p0)
        err_v = np.linalg.norm(v-v0) / np.linalg.norm(v0)
        errors_p.append(err_p)
        errors_v.append(err_v)
    print('Code port error in pos:', tools.stats(np.array(errors_p)))
    print('Code port error in vel:', tools.stats(np.array(errors_v)))

def test_planet_ephemeris_against_jpl_de(planet_ephemeris, jpl_pos_vel, num):
    """
    Compares output of the given ephemeris function `planet_ephemeris` to the JPL DE ephemeris 
    function given by `jpl_pos_vel`. Each planet is tested with times evenly sampled from
    multiple time intervals.
    """
    test_times = generate_test_times(num)
    print(f'{'BODY':>10}{'INTERVAL':>16}{'MEAN_POS_ERR':>15}{'STDDEV_POS_ERR':>15}{'MAX_POS_ERR':>15}{'MAX_VEL_ERR':>15}')
    for body_name in JPL_DE_ROUTES.keys():
        if body_name in ['PLUTO', 'EMB-MOON', 'EMB-EARTH']:
            continue
        for interval in test_times.keys():
            errors_p = []
            errors_v = []
            for t in test_times[interval]:
                p0_earth, _ = jpl_pos_vel([0, 3, 399], t)
                p0, v0 = jpl_pos_vel(JPL_DE_ROUTES[body_name], t)
                p, v = planet_ephemeris(body_name, t)
                error_p_sun = np.linalg.norm(p-p0) / np.linalg.norm(p0)
                error_p_earth = np.linalg.norm(p-p0) / np.linalg.norm(p0-p0_earth)
                error_v = np.linalg.norm(v-v0) / np.linalg.norm(v0)
                errors_p.extend([error_p_sun, error_p_earth] if body_name != 'EARTH-MOON' else [error_p_sun])
                errors_v.append(error_v)
            err_0 = tools.dms_string(tools.stats(np.array(errors_p))['mean'], 2)
            err_stddev = tools.dms_string(tools.stats(np.array(errors_p))['stdDev'], 2)
            err_1 = tools.dms_string(tools.stats(np.array(errors_p))['max'], 2)
            err_2 = tools.stats(np.array(errors_v))['max']
            print(f'{body_name:>10}{interval:>16}{str(err_0):>15}{str(err_stddev):>15}{str(err_1):>15}{f'{err_2:.0e}':>15}')

def compare_pos_vel_functions(pos_vel, pos_vel_ref, num, print_header=True):
    """
    Compares output of two functions that produce pos, vel vectors as a function of time.
    """
    test_times = generate_test_times(num)
    if print_header:
        print(f'{'INTERVAL':>16}{'MEAN_POS_ERR':>15}{'STDDEV_POS_ERR':>15}{'MAX_POS_ERR':>15}{'MAX_VEL_ERR':>15}')
    for interval in test_times.keys():
        errors_p = []
        errors_v = []
        for t in test_times[interval]:
            p, v = pos_vel(t)
            p0, v0 = pos_vel_ref(t)
            err_p = np.linalg.norm(p-p0) / np.linalg.norm(p0)
            err_v = np.linalg.norm(v-v0) / np.linalg.norm(v0)
            errors_p.append(err_p)
            errors_v.append(err_v)
        err_0 = tools.dms_string(tools.stats(np.array(errors_p))['mean'], 2)
        err_stddev = tools.dms_string(tools.stats(np.array(errors_p))['stdDev'], 2)
        err_1 = tools.dms_string(tools.stats(np.array(errors_p))['max'], 2)
        err_2 = tools.stats(np.array(errors_v))['max']
        print(f'{interval:>16}{str(err_0):>15}{str(err_stddev):>15}{str(err_1):>15}{f'{err_2:.0e}':>15}')

def run_tests(test_num):
    vsop87 = VSOP87Ephemeris(R'./json/vsop87a_raw.json')
    mpp02_llr = MPP02Ephemeris(R'./json/mpp02_llr_raw.json')
    mpp02_405 = MPP02Ephemeris(R'./json/mpp02_405_raw.json')

    print('\n', '-'*20, 'VSOP87A: CODE PORT TESTS AGAINST vsop87.chk', '-'*20)
    run_vsop87_checks(vsop87)

    print('\n', '-'*20, 'ELP/MPP02: CODE PORT TESTS AGAINST FORTRAN OUTPUT', '-'*20)
    run_mpp02_tests(mpp02_llr, mpp02_405)

    with tools.jplephem_pos_vel(JPL_DE_EPHEMERIS_PATH) as jpl_pos_vel:
        print('\n', '-'*20, 'RAW VSOP87A ERROR AGAINST DE441', '-'*20)
        test_planet_ephemeris_against_jpl_de(vsop87.get_pos_vel, jpl_pos_vel, test_num)

        print('\n', '-'*20, 'RAW ELP/MPP02(LLR) ERROR AGAINST DE441', '-'*20)
        jpl_moon_pos_vel = lambda t: jpl_pos_vel([399, 3, 301], t)
        compare_pos_vel_functions(mpp02_llr.get_pos_vel, jpl_moon_pos_vel, test_num)

if __name__ == '__main__':
    run_tests(100)