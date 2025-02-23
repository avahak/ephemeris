class FixedLengthReader:
    """
    Reads files with fixed length format.
    Usage example: 
    ```
    reader = FixedLengthReader('3,i1,i1,10,12i3,x1,2f18,2f14,f20')
    coord, alpha, *a_coeffs, s, k, a, b, c = reader.read(line)
    ```
    """
    def __init__(self, format_string):
        """
        Initialize the FixedLengthReader with a comma-separated format string.
        Each part can be an integer (skip to index) or of the form `mfl`, where:
        - `m` is multiplicity (optional),
        - `f` is field type (`i` for int, `f` for float, `s` for str, `x` for skip),
        - `l` is length (optional).
        """
        self.format_list = self._parse_format_string(format_string.lower())

    def _parse_field_spec(self, part):
        """
        Parse a field specification like '13i3' or 's' into (multiplicity, field_type, length).
        """
        multiplicity = ''
        i = 0
        while i < len(part) and part[i].isdigit():
            multiplicity += part[i]
            i += 1
        multiplicity = int(multiplicity) if multiplicity else 1

        field_type = ''
        while i < len(part) and part[i].isalpha():
            field_type += part[i]
            i += 1
        field_type = { 'x': None, 'i': int, 'f': float, 's': str }[field_type]

        length = ''
        while i < len(part) and part[i].isdigit():
            length += part[i]
            i += 1
        length = int(length) if length else None

        return multiplicity, field_type, length

    def _parse_format_string(self, format_string):
        """
        Parse the format string into a list of (start_index, field_type, length) entries.
        """
        format_list = []
        parts = format_string.split(sep=',')
        current_index = 0

        for part_index, part in enumerate(parts):
            # Handle explicit index (e.g., '78')
            if part.isdigit():
                current_index = int(part)
                continue

            multiplicity, field_type, length = self._parse_field_spec(part)

            if length is None and part_index + 1 < len(parts):
                try:
                    total_length = int(parts[part_index + 1]) - current_index
                except ValueError:
                    raise ValueError(f'Implicit multiplicity cannot be computed: end index not specified.')
                if total_length % multiplicity != 0:
                    raise ValueError(f'Implicit multiplicity cannot be computed: length {length} is not divisible by multiplicity {multiplicity}.')
                length = total_length // multiplicity

            # Add entries for each element in the multiplicity
            for _ in range(multiplicity):
                if field_type is not None:
                    format_list.append((current_index, field_type, length))
                if length is not None:
                    current_index += length

        return format_list

    def read(self, line, strip_strings=True):
        """
        Parse a line of text using the format list and return a list of parsed values.
        """
        result = []
        for start_index, field_type, length in self.format_list:
            substring = line[start_index : start_index+length] if length is not None else line[start_index:]
            if field_type == str:
                result.append(substring.strip() if strip_strings else substring)
            elif field_type == float:
                # Replace d with e (Fortran double)
                result.append(field_type(substring.strip().replace('d', 'e').replace('D', 'e')))
            else:
                result.append(field_type(substring.strip()))

        return result

def _test():
    line = ' 1932   14  0  0  1  0 -1  0  0  0  0  0  0  0 -0.00000104352     0.00000676414     0.00000684416 1.30699021227    5753.38488489680 '
    reader = FixedLengthReader('3,i1,i1,10,12i3,x1,2f18,2f14,f')
    v = reader.read(line)
    print(v)

if __name__ == '__main__':
    _test()