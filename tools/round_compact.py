"""
Rounding floating point numbers and writing them in a compact form.

NOTE Implementations are slow so use Python rounding f'{x:.3e}' if compactness is not important.
"""
import json
import re
from decimal import Decimal
import numpy as np

class FormattedFloat:
    """
    Stores a specific string representation `as_str` for a float. 
    The value as float is still `float(a_str)`.
    """
    ZERO: "FormattedFloat" = None    # Assigned to '0' after

    def __init__(self, as_str: str):
        if not isinstance(as_str, str):
            raise ValueError('Argument should be formatted string representing the float.')
        self.__as_str = as_str
        self.__as_float = float(as_str)

    @property
    def as_str(self) -> str:
        return self.__as_str

    @property
    def as_float(self) -> float:
        return self.__as_float

    def __repr__(self):
        return f'FormattedFloat({self.__as_str})'
    
    def __len__(self):
        return len(self.__as_str)
    
    def __float__(self):
        return self.__as_float
    
FormattedFloat.ZERO = FormattedFloat('0')

class FormattedFloatEncoder(json.JSONEncoder):
    """
    JSON Encoder for FormattedFloat. Use apply_formatting after json encoding
    to apply the formatting. 

    Example:
    ```
    data = [8.1e6, FormattedFloat('8.1e6')]
    json_str = json.dumps(data, cls=FormattedFloatEncoder)
    json_str = FormattedFloatEncoder.apply_formatting(json_str)
    print(f'{json_str = }')     # json_str = '[8100000.0, 8.1e6]'
    ```
    """
    def default(self, obj):
        if isinstance(obj, FormattedFloat):
            return str(obj)
        return super().default(obj)
    
    @staticmethod
    def apply_formatting(json_str: str) -> str:
        """
        Replaces all occurrences of "FormattedFloat(...)" in the JSON string
        with the number inside the parentheses.
        """
        pattern = r'"FormattedFloat\(([^)]+)\)"'
        return re.sub(pattern, r'\1', json_str)

def round_standard_form_components(x, sig_figs):
    """
    Returns the components `(sign, mantissa, exponent)` of x rounded to sig_figs significant
    figures. The approximation for x is 
    `sign * (first digit of mantissa).(rest of the digits of mantissa) * 10^exponent`
    - `sign` is -1 if x<0 and 1 if x>0
    - `mantissa` is integer starting with digit 1..9 and has sig_figs digits
    - `exponent` is the exponent of the rounded number 
    Rounding to nearest, ties to even.
    """
    if x == 0:
        return 0, 0, 0
    _, digits, exponent = Decimal(x).as_tuple()
    digits = list(digits)
    if len(digits) < sig_figs + 1:
        # pad with zeros if we don't start with enough digits
        pad = sig_figs + 1 - len(digits)
        digits.extend([0] * pad)
        exponent -= pad
    exponent = exponent + len(digits) - 1

    mantissa = 0
    for k in range(sig_figs):
        mantissa = 10*mantissa + digits[k]

    # rounding rules
    if digits[sig_figs] > 5:
        mantissa += 1
    if digits[sig_figs] == 5:
        is_halfway = True
        for digit in digits[sig_figs+1:]:
            if digit != 0:
                is_halfway = False
        if not is_halfway:
            mantissa += 1
        elif digits[sig_figs-1] % 2 == 1:
            # ties to even
            mantissa += 1

    mantissa = str(mantissa)
    if len(mantissa) > sig_figs:
        # we rounded up and ended up with integer part 10, so normalize back to [1,10)
        mantissa = mantissa[:sig_figs]
        exponent += 1

    assert len(mantissa) == sig_figs
    return int(np.sign(x)), int(mantissa), exponent

def _compose_number(sign, mantissa_digits, decimal_place, exponent):
    # Writes out the number from its components
    sign_part = '-' if sign == -1 else ''
    if len(mantissa_digits) > decimal_place:
        decimal_part = f'{mantissa_digits[:decimal_place]}.{mantissa_digits[decimal_place:]}'
    else:
        decimal_part = f'{mantissa_digits[:decimal_place]}'
    exponent_part = f'e{exponent}' if exponent != 0 else ''
    return f'{sign_part}{decimal_part}{exponent_part}'

def round_compact(x: float, sig_figs: int) -> FormattedFloat:
    """
    Rounds x to sig_figs significant figures. Returns the result as a FormattedFloat that has 
    the most compact form found while still indicating exactly sig_figs significant figures.
    Returns `FormattedFloat.ZERO` if sig_figs is 0.
    NOTE Rounding to sig_figs=17 does not change the number.
    """
    if sig_figs < 0:
        raise ValueError('Invalid number of significant figures.')
    if sig_figs == 0:
        return FormattedFloat.ZERO
    sig_figs = min(sig_figs, 17)
    if x == 0:
        s = '0' if sig_figs == 1 else '0.' + ('0' * (sig_figs-1))
        return FormattedFloat(s)
    
    sign, mantissa, exponent = round_standard_form_components(x, sig_figs)
    mantissa = str(mantissa)
    exponent_length = len(str(exponent))

    decimal_place = 1

    # First attempt at rounding: standard form
    x_rounded_1 = _compose_number(sign, mantissa, decimal_place, exponent)

    # Attempts at writing the number without the exponent
    if exponent > 0 and exponent < sig_figs:
        # remove exponent and just move the decimal point
        decimal_place += exponent
        exponent = 0
    if exponent >= sig_figs and exponent <= sig_figs+exponent_length+6:
        # remove exponent and pad with zeros from right 
        decimal_place += exponent
        mantissa += '0' * (exponent-sig_figs+1)
        exponent = 0
    if exponent >= -6-exponent_length and exponent < 0:
        # add leading zeroes to get rid of the exponent
        mantissa = '0'*(-exponent) + mantissa
        exponent = 0

    # Second attempt at rounding
    x_rounded_2 = _compose_number(sign, mantissa, decimal_place, exponent)

    # Pick the more compact form with preference to avoiding exponent part
    x_rounded = x_rounded_1 if len(x_rounded_1) < len(x_rounded_2) else x_rounded_2

    # Check that rounding is correct for both tries
    assert float(x_rounded_1) == float(x_rounded_2)
    assert float(f'{x:.{sig_figs-1}e}') == float(x_rounded)

    return FormattedFloat(x_rounded)

def round_compact_in_interval(x: float, a: float, b: float) -> FormattedFloat:
    """
    Rounds `x` to varying significant digits and picks the rounded value that 
    has shortest string representation while still in the interval `[a,b]`. 
    Returns the picked string representation.
    """
    if (a > x) or (x > b):
        raise ValueError('Invalid parameters: a<=x<=b is not satisfied.')
    # Now a<=x<=b
    best_found = FormattedFloat(str(x))
    for sig_figs in range(17, -1, -1):
        x_round = round_compact(x, sig_figs)
        if (float(x_round) >= a) and (float(x_round) <= b) and (len(x_round) < len(best_found)):
            best_found = x_round
    return best_found

def _random_number(): 
    # This is just for testing
    exponent = np.random.rand()*30.0 - 15.0
    num = 2.0*np.random.rand() - 1.0
    if np.random.rand() < 0.02:
        return 0.0
    x = num*np.power(10.0, exponent)
    if np.random.rand() < 0.5:
        return float(f'{x:.{int(20*np.random.rand())}e}')
    return x

def _test_round_compact():
    # Test round_compact with varying sig_figs
    tests = [_random_number() for _ in range(100)]
    for t in tests:
        for sig_num in range(1, 20):
            print(f'{sig_num=}, {t=}')
            print('python .e:', f'{t:.{sig_num-1}e}')
            result = round_compact(t, sig_num)
            print('17', float(round_compact(t, 17)))
            assert t == float(round_compact(t, 17))
            print('   result:', result)
            print('---')

def _test_round_compact_in_interval():
    # Test round_compact_in_interval with varying numbers and intervals
    for _ in range(100):
        r = [_random_number(), _random_number(), _random_number()]
        a, x, b = sorted(r)

        x_round = round_compact_in_interval(x, a, b)
        x_round_f = float(x_round)
        assert (a <= x_round_f) and (x_round_f <= b)

        print(f'{[x, a, b] = }')
        print(f'{x_round = }')
        print('---')
        

def _test_formatted_float_encoder():
    # Test FormattedFloatEncoder
    data = [FormattedFloat(str(_random_number())) for _ in range(5)]
    data.append({ 'key1': FormattedFloat('1.000'), 'key2': FormattedFloat('3.14e-9') })
    print(f'{data = }')
    print()
    json_str = json.dumps(data, cls=FormattedFloatEncoder, indent=4)
    print('Before:\n', json_str)
    print('After:\n', FormattedFloatEncoder.apply_formatting(json_str))

    print('---')

    data = [8.1e6, FormattedFloat('8.1e6')]
    json_str = json.dumps(data, cls=FormattedFloatEncoder)
    json_str = FormattedFloatEncoder.apply_formatting(json_str)
    print(f'{json_str = }')

if __name__ == '__main__':
    _test_round_compact()
    _test_round_compact_in_interval()
    _test_formatted_float_encoder()