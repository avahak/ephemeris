"""
Reads coefficients in VSOP87A files and stores the series in one json file. 
Also reads the accompanying vsop87.chk file that lists expected (x,y,z), (x',y',z') 
output of the series and stores the output in another json file.
"""
import json

from fixed_length_reader import FixedLengthReader

INPUT_DIRECTORY = 'd:/resources/astro/vsop87'
TEST_OUTPUT_PATH = './json/vsop87a_tests.json'
OUTPUT_RAW_JSON_PATH = './json/vsop87a_raw.json'

VSOP87_NAME_INDEXING = { 
    'MERCURY': (1, 'VSOP87A.mer'), 
    'VENUS': (2, 'VSOP87A.ven'), 
    'EARTH': (3, 'VSOP87A.ear'),
    'MARS': (4, 'VSOP87A.mar'), 
    'JUPITER': (5, 'VSOP87A.jup'), 
    'SATURN': (6, 'VSOP87A.sat'), 
    'URANUS': (7, 'VSOP87A.ura'), 
    'NEPTUNE': (8, 'VSOP87A.nep'), 
    'EARTH-MOON': (9, 'VSOP87A.emb')
}

# List of mean longitudes.
# LAMBDA = [
#     [4.40260884240, 26087.9031415742], 
#     [3.17614669689, 10213.2855462110],
#     [1.75347045953, 6283.0758499914],
#     [6.20347611291, 3340.6124266998],
#     [0.59954649739, 529.6909650946],
#     [0.87401675650, 213.2990954380],
#     [5.48129387159, 74.7815985673],
#     [5.31188628676, 38.1330356378],
#     [5.19846674103, 77713.7714681205],
#     [1.62790523337, 84334.6615813083],
#     [2.35555589827, 83286.9142695536],
#     [3.81034454697, 83997.0911355954]
# ]

# Transformation from coordinate system used in VSOP87A to FK5 (~ICRF).
VSOP87_MATRIX = [
    [ 1.0,            0.00000044036,  -0.000000190919],
    [-0.000000479966, 0.917482137087, -0.397776982902],
    [ 0.0,            0.397776982902,  0.917482137087]
]

def load_chk_tests(): 
    # Loads expected series output values from vsop87.chk.
    with open(f'{INPUT_DIRECTORY}/vsop87.chk', 'r') as f:
        lines = f.readlines()

    tests = []
    header_reader = FixedLengthReader('s8,x2,s12,x2,f11')
    vector_reader = FixedLengthReader('3,f14,30,f14,57,f14')
    for k, line in enumerate(lines):
        if line[1:8] != 'VSOP87A':
            continue
        _, body_name, jd = header_reader.read(line)
        body_index = VSOP87_NAME_INDEXING[body_name][0]
        p = vector_reader.read(lines[k+1])
        v = vector_reader.read(lines[k+2])
        tests.append({ 'body': body_name, 'body_index': body_index, 'jd': jd, 'p': p, 'v': v })
    return tests

def load_raw_data(): 
    # Reads coefficients from the input file
    bodies = {}
    term_count = {}
    reader = FixedLengthReader('3,i1,i1,10,12i3,x1,2f18,2f14,f20')
    for body_name, (_, file_name) in VSOP87_NAME_INDEXING.items():
        with open(f'{INPUT_DIRECTORY}/{file_name}', 'r') as f:
            for line in f:
                try: 
                    if line[1:7] == 'VSOP87':   # header lines (no coefficients)
                        continue
                    coord, alpha, *a_coeffs, s, k, a, b, c = reader.read(line)
                    # We only store a, b, c as they are enough for pos, vel
                    entry = [a, b, c]
                    groups: dict = bodies.setdefault(body_name, {})
                    groups.setdefault((coord, alpha), []).extend(entry)
                    term_count[body_name] = term_count.get(body_name, 0) + 1
                except ValueError as e:
                    print(e)

    print(f'total # of terms = {sum(term_count.values())},\n{term_count = }')

    # Change the structure of groups by making a list instead of dict under each body.
    # Intention is to make the object more easily readable.
    bodies_restructured = {}
    for body_name, groups in bodies.items():
        groups_restructured = [{ 'coord': key[0], 'alpha': key[1], 'coeffs': value } for key, value in groups.items()]
        bodies_restructured[body_name] = groups_restructured

    return { 
        '_comment': 'VSOP87A raw data with coefficients a, b, c',
        # 'lambda': LAMBDA, 
        'matrix': VSOP87_MATRIX, 
        'bodies': bodies_restructured 
    }

def write_raw_json(obj_raw):
    # Writes the raw data into a json file
    with open(OUTPUT_RAW_JSON_PATH, 'w') as f:
        json_string = json.dumps(obj_raw, indent=None, separators=(',', ':'))
        f.write(json_string)
    print(f'Wrote raw json to {OUTPUT_RAW_JSON_PATH}.')

def write_chk_tests(tests):
    # Writes the tests into a json file
    with open(TEST_OUTPUT_PATH, 'w') as f:
        json.dump(tests, f, indent=None, separators=(',', ':'))
    print(f'Wrote tests json to {TEST_OUTPUT_PATH}.')

if __name__ == '__main__':
    print('-'*20, 'Loading files, writing raw json', '-'*20)
    obj_raw = load_raw_data()
    write_raw_json(obj_raw)

    print('-'*20, 'Loading .chk tests, writing them to json', '-'*20)
    tests = load_chk_tests()
    write_chk_tests(tests)