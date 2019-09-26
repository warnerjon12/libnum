#-*- coding:utf-8 -*-

import math
import random

from functools import reduce
from .compat import xrange, math_log2


def len_in_bits(n):
    """
    Return number of bits in binary representation of @n.
    """
    try:
        return n.bit_length() # new in Python 2.7
    except AttributeError:
        if n == 0:
            return 0
        return int(math_log2(n)) + 1


def randint_bits(size):
    low = 1 << (size - 1)
    hi = (1 << size) - 1
    return random.randint(low, hi)


def ceil(x, y):
    q, r = divmod(x, y)
    return q + int(r != 0)


def nroot(x, n):
    """
    Return truncated n'th root of x.
    """
    if n < 0:
        raise ValueError("can't extract negative root")

    if n == 0:
        raise ValueError("can't extract zero root")

    sign = 1
    if x < 0:
        sign = -1
        x = -x
        if n % 2 == 0:
            raise ValueError("can't extract even root of negative")

    high = 1
    while high ** n <= x:
        high <<= 1

    low = high >> 1
    while low < high:
        mid = (low + high) >> 1
        if low < mid and mid ** n < x:
            low = mid
        elif high > mid and mid ** n > x:
            high = mid
        else:
            return sign * mid
    return sign * (mid + 1)


try:
    _gcd = math.gcd
except AttributeError:

    def _gcd(a, b):
        """
        Return greatest common divisor using Euclid's Algorithm.
        """
        if a == 0: return abs(b)
        while b:
            a, b = b, a % b
        return abs(a)


def _lcm(a, b):
    """
    Return lowest common multiple.
    """
    if not a or not b:
        raise ZeroDivisionError("lcm arguments may not be zeros")
    return abs(a * b) // _gcd(a, b)


def gcd(*lst):
    """
    Return gcd of a variable number of arguments.
    """
    return abs(reduce(lambda a, b: _gcd(a, b), lst))


def lcm(*lst):
    """
    Return lcm of a variable number of arguments.
    """
    return reduce(lambda a, b: _lcm(a, b), lst)


def xgcd(a, b):
    """
    Extented Euclid GCD algorithm.
    Return (x, y, g) : a * x + b * y = gcd(a, b) = g.
    """
    if a == 0: return 0, 1, b

    px, ppx = 0, 1
    py, ppy = 1, 0

    while b:
        q = a // b
        a, b = b, a % b
        x = ppx - q * px
        y = ppy - q * py
        ppx, px = px, x
        ppy, py = py, y

    return ppx, ppy, a


def extract_power2(a):
    """
    Optimized version of extract_prime_power for p=2
    """
    s = 0
    while True:
        if a & 1:
            break
        a >>= 1
        s += 1
    return s, a


def extract_prime_power(a, p):
    """
    Return s, t such that  a = p**s * t,  t % p = 0
    """
    s = 0
    while True:
        q, r = divmod(a, p)
        if r != 0:
            break
        a = q
        s += 1
    return s, a


def solve_linear(a, b, c):
    """
    Solve a*x + b*y = c.
    Solution (x0 + b*n, y0 - a*n).
    Return None or (x0, y0).
    """
    x, y, g = xgcd(a, b)
    try:
        q, r = divmod(c, g)
    except ZeroDivisionError: # i.e. both a and b are zero
        return (0, 0) if c == 0 else None
    if r != 0:
        return None
    return x * q, y * q
