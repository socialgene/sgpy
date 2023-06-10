from math import floor, log10


def find_exp(number) -> int:
    number = float(number)
    if number < 0:
        raise ValueError("Function only meant for positive numbers")
    if number == 0:
        return 0
    base10 = log10(abs(number))
    return int(floor(base10))
