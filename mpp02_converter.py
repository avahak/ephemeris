"""
Creates mpp02_405_raw.json and mpp02_llr_raw.json files that contain nontruncated
coefficients for ELP/MPP02 series (fitted to DE405 or LLR data).
This is mostly just a direct code port of the Fortran code for initialization and 
file reading from ELPMPP02.for. Fortran code has comments that are not copied here.

ELP/MMP02 article: Lunar solution ELP, version ELP/MPP02, Jean Chapront and Gerard Francou. 
Fortran code and original series data files from: 
ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/
"""
import numpy as np
import json

from tools.fixed_length_reader import FixedLengthReader

RAD = 648000.0 / np.pi
DEG = np.pi / 180.0
SERIES_DIRECTORY = R'd:/resources/astro/elp_mpp02'   # directory for ELP_MAIN.S1, etc.
OUTPUT_DIRECTORY = R'./json'

def dms(ideg, imin, sec):
    return (ideg + imin/60.0 + sec/3600.0) * DEG

def load_series_and_write_json(mode):
    # Code here is direct Fortran port with almost same variable names.
    if mode != '405':
        mode = 'LLR'

    max1 = 2645
    max2 = 33256
    nmpb = np.zeros((3, 3), dtype=int)
    cmpb = np.zeros(max1)
    fmpb = np.zeros((5, max1))
    nper = np.zeros((3, 4, 3), dtype=int)
    cper = np.zeros(max2)
    fper = np.zeros((5, max2))

    w = np.zeros((3, 5))
    eart = np.zeros(5)
    peri = np.zeros(5)
    zeta = np.zeros(5)
    dela = np.zeros((4, 5))
    p = np.zeros((8, 5))

    # NOTE Careful with code port from Fortran here! 
    # In Fortran arrays arrays are stored in column-major order.
    bp = np.array([
        [0.311079095, -0.103837907],
        [-0.4482398e-2, 0.668287e-3],
        [-0.1102485e-2, -0.1298072e-2],
        [0.1056062e-2, -0.178028e-3],
        [0.50928e-4, -0.37342e-4]
    ])

    Dprec = -0.29965
    am =  0.074801329
    alpha = 0.002571881
    dtasm = (2.0*alpha) / (3.0*am)
    xa = (2.0*alpha) / 3.0

    if mode == 'LLR':
        Dw1_0   = -0.10525
        Dw2_0   =  0.16826
        Dw3_0   = -0.10760
        Deart_0 = -0.04012
        Dperi   = -0.04854
        Dw1_1   = -0.32311
        Dgam    =  0.00069
        De      =  0.00005
        Deart_1 =  0.01442
        Dep     =  0.00226
        Dw2_1   =  0.08017
        Dw3_1   = -0.04317
        Dw1_2   = -0.03794
    else:
        Dw1_0   = -0.07008
        Dw2_0   =  0.20794
        Dw3_0   = -0.07215
        Deart_0 = -0.00033
        Dperi   = -0.00749
        Dw1_1   = -0.35106
        Dgam    =  0.00085
        De      = -0.00006
        Deart_1 =  0.00732
        Dep     =  0.00224
        Dw2_1   =  0.08017
        Dw3_1   = -0.04317
        Dw1_2   = -0.03743

    w[0, 0] = dms(218, 18, 59.95571 + Dw1_0)
    w[0, 1] = (1732559343.73604 + Dw1_1) / RAD
    w[0, 2] = (-6.8084 + Dw1_2) / RAD
    w[0, 3] = 0.66040e-2 / RAD
    w[0, 4] = -0.31690e-4 / RAD

    w[1, 0] = dms(83, 21, 11.67475 + Dw2_0)
    w[1, 1] = (14643420.3171 + Dw2_1) / RAD
    w[1, 2] = -38.2631 / RAD
    w[1, 3] = -0.45047e-1 / RAD
    w[1, 4] = 0.21301e-3 / RAD

    w[2, 0] = dms(125, 2, 40.39816 + Dw3_0)
    w[2, 1] = (-6967919.5383 + Dw3_1) / RAD
    w[2, 2] = 6.3590 / RAD
    w[2, 3] = 0.76250e-2 / RAD
    w[2, 4] = -0.35860e-4 / RAD

    eart[0] = dms(100, 27, 59.13885 + Deart_0)
    eart[1] = (129597742.29300 + Deart_1) / RAD
    eart[2] = -0.020200 / RAD
    eart[3] = 0.90000e-5 / RAD
    eart[4] = 0.15000e-6 / RAD

    peri[0] = dms(102, 56, 14.45766 + Dperi)
    peri[1] = 1161.24342 / RAD
    peri[2] = 0.529265 / RAD
    peri[3] = -0.11814e-3 / RAD
    peri[4] = 0.11379e-4 / RAD

    if mode == '405':
        w[0, 3] += -0.00018865 / RAD
        w[0, 4] += -0.00001024 / RAD
        w[1, 2] +=  0.00470602 / RAD
        w[1, 3] += -0.00025213 / RAD
        w[2, 2] += -0.00261070 / RAD
        w[2, 3] += -0.00010712 / RAD

    x2 = w[1, 1] / w[0, 1]
    x3 = w[2, 1] / w[0, 1]
    y2 = am*bp[0, 0] + xa*bp[4, 0]
    y3 = am*bp[0, 1] + xa*bp[4, 1]

    d21 = x2 - y2
    d22 = w[0, 1] * bp[1, 0]
    d23 = w[0, 1] * bp[2, 0]
    d24 = w[0, 1] * bp[3, 0]
    d25 = y2 / am

    d31 = x3 - y3
    d32 = w[0, 1] * bp[1, 1]
    d33 = w[0, 1] * bp[2, 1]
    d34 = w[0, 1] * bp[3, 1]
    d35 = y3 / am

    Cw2_1 = d21*Dw1_1 + d25*Deart_1 + d22*Dgam + d23*De + d24*Dep
    Cw3_1 = d31*Dw1_1 + d35*Deart_1 + d32*Dgam + d33*De + d34*Dep

    w[1, 1] += Cw2_1 / RAD
    w[2, 1] += Cw3_1 / RAD
        
    for i in range(5):
        dela[0, i] = w[0, i] - eart[i]
        dela[1, i] = w[0, i] - w[2, i]
        dela[2, i] = w[0, i] - w[1, i]
        dela[3, i] = eart[i] - peri[i]
    dela[0, 0] += np.pi

    p[0, 0] = dms(252, 15, 3.216919)
    p[1, 0] = dms(181, 58, 44.758419)
    p[2, 0] = dms(100, 27, 59.138850)
    p[3, 0] = dms(355, 26, 3.642778)
    p[4, 0] = dms(34, 21, 5.379392)
    p[5, 0] = dms(50, 4, 38.902495)
    p[6, 0] = dms(314, 3, 4.354234)
    p[7, 0] = dms(304, 20, 56.808371)

    p[0, 1] = 538101628.66888 / RAD
    p[1, 1] = 210664136.45777 / RAD
    p[2, 1] = 129597742.29300 / RAD
    p[3, 1] = 68905077.65936 / RAD
    p[4, 1] = 10925660.57335 / RAD
    p[5, 1] = 4399609.33632 / RAD
    p[6, 1] = 1542482.57845 / RAD
    p[7, 1] = 786547.89700 / RAD

    zeta[0] = w[0, 0]
    zeta[1] = w[0, 1] + (5029.0966 + Dprec) / RAD
    zeta[2] = w[0, 2]
    zeta[3] = w[0, 3]
    zeta[4] = w[0, 4]

    delnu = (0.55604 + Dw1_1) / RAD / w[0, 1]
    dele = (0.01789 + De) / RAD
    delg = (-0.08066 + Dgam) / RAD
    delnp = (-0.06424 + Deart_1) / RAD / w[0, 1]
    delep = (-0.12879 + Dep) / RAD

    # File reading

    reader_main_header = FixedLengthReader('25,i10')
    reader_main_data = FixedLengthReader('4i3,x2,f13,5f12')
    reader_perturbations_header = FixedLengthReader('25,2i10')
    reader_perturbations_data = FixedLengthReader('i5,2f20,16i3')
    ilu = np.zeros(4, dtype=int)
    b = np.zeros(5)
    ifi = np.zeros(16, dtype=int)

    ir = 0
    for iv in range(3):
        file_path = f'{SERIES_DIRECTORY}/ELP_MAIN.S{iv+1}'
        with open(file_path, 'r') as f:
            line_data = reader_main_header.read(f.readline())
            nmpb[iv, 0] = line_data[0]
            nmpb[iv, 1] = ir
            nmpb[iv, 2] = nmpb[iv, 0] + nmpb[iv, 1]
            for _ in range(nmpb[iv, 0]):
                line_data = reader_main_data.read(f.readline())
                ilu[0], ilu[1], ilu[2], ilu[3], a, *b = line_data
                tgv = b[0] + dtasm*b[4]
                if iv == 2:
                    a -= 2.0*a*delnu/3.0
                cmpb[ir] = a + tgv*(delnp - am*delnu) + b[1]*delg + b[2]*dele + b[3]*delep
                for k in range(5):
                    fmpb[k, ir] = 0.0
                    for i in range(4):
                        fmpb[k, ir] += ilu[i] * dela[i, k]
                if iv == 2: 
                    fmpb[0, ir] += np.pi/2.0
                ir += 1

    ir = 0
    for iv in range(3):
        file_path = f'{SERIES_DIRECTORY}/ELP_PERT.S{iv+1}'
        with open(file_path, 'r') as f:
            for it in range(4):
                line_data = reader_perturbations_header.read(f.readline())
                nper[iv, it, 0], _it = line_data
                nper[iv, it, 1] = ir
                nper[iv, it, 2] = nper[iv, it, 0] + nper[iv, it, 1]
                assert _it == it
                for _ in range(nper[iv, it, 0]):
                    line_data = reader_perturbations_data.read(f.readline())
                    _, s, c, *ifi = line_data
                    cper[ir] = np.sqrt(c*c + s*s)
                    pha = np.arctan2(c, s)
                    if pha < 0.0: 
                        pha += 2.0*np.pi
                    for k in range(5):
                        fper[k, ir] = 0.0
                        if k == 0: 
                            fper[k, ir] = pha
                        for i in range(4):
                            fper[k, ir] += ifi[i] * dela[i, k]
                        for i in range(4, 12):
                            fper[k, ir] += ifi[i] * p[i-4, k]
                        fper[k, ir] += ifi[12] * zeta[k]
                    ir += 1

    write_json(mode, w[0], nmpb, cmpb, fmpb, nper, cper, fper)


def write_json(mode, w0, nmpb, cmpb, fmpb, nper, cper, fper):
    # Writes the series data to a json file.
    print('----- write_json -----')
    groups = []
    for coord in range(3):
        # Main Problem
        coeffs_main = [(cmpb[ir], *[fmpb[k, ir] for k in range(5)]) for ir in range(nmpb[coord, 1], nmpb[coord, 2])]
        for alpha in range(4):
            # Perturbations
            coeffs_pert = [(cper[ir], *[fper[k, ir] for k in range(5)]) for ir in range(nper[coord, alpha, 1], nper[coord, alpha, 2])]
            # Combine the terms from both sources
            coeffs = [*coeffs_main, *coeffs_pert] if alpha == 0 else coeffs_pert
            if len(coeffs) > 0:
                groups.append({ 'coord': coord, 'alpha': alpha, 'coeffs': np.array(coeffs).flatten().tolist() })

    obj = { '_comment': f'ELP/MPP02({mode})', 'W': w0.tolist(), 'groups': groups }
    file_path = f'{OUTPUT_DIRECTORY}/mpp02_{mode.lower()}_raw.json'
    with open(file_path, 'w') as f:
        json.dump(obj, f, indent=None, separators=(',', ':'))
    print(f'Wrote to file {file_path}.')


if __name__ == '__main__':
    load_series_and_write_json('405')
    load_series_and_write_json('LLR')