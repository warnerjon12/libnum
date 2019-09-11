#-*- coding:utf-8 -*-

import math
import random

from functools import reduce


def len_in_bits(n):
    """
    Return number of bits in binary representation of @n.
    """
    try:
        return n.bit_length() # new in Python 2.7
    except AttributeError:
        if n == 0:
            return 0
        return math.floor(math.log2(abs(n))) + 1


def randint_bits(size):
    low = 1 << (size - 1)
    hi = (1 << size) - 1
    return random.randint(low, hi)


def ceildiv(x, y):
    q, r = divmod(x, y)
    return q + int(r != 0)


def nroot(x, n):
    """
    Return the floor of the n'th root of x.
    """
    if n <= 0:
        raise ValueError("can't extract %s root" % ("negative"
                                                    if n < 0
                                                    else "zero"))

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


def _gcd(a, b):
    """
    Return the greatest common divisor of a and b using Euclid's Algorithm.
    """
    if a == 0: return abs(b)
    while b:
        a, b = b, a % b
    return abs(a)


def _lcm(a, b):
    """
    Return the least common multiple of a and b.
    """
    if not a or not b:
        raise ZeroDivisionError("lcm arguments may not be zeros")
    return abs(a * b) // _gcd(a, b)


def gcd(*lst):
    """
    Return gcd of a variable number of arguments.
    """
    return abs(reduce(_gcd, lst))


def lcm(*lst):
    """
    Return lcm of a variable number of arguments.
    """
    return abs(reduce(_lcm, lst))


def xgcd(a, b):
    """
    Extended Euclidian GCD algorithm.
    Return (x, y, g) : a * x + b * y = gcd(a, b) = g.
    """
    if a == 0: return 0, 1, b

    px, ppx = 0, 1
    py, ppy = 1, 0

    while b:
        a, q, b = b, *divmod(a, b)
        x = ppx - q * px
        y = ppy - q * py
        ppx, px = px, x
        ppy, py = py, y

    return ppx, ppy, a


def extract_power(a, p):
    """
    Return s, t such that  a = p**s * t,  t % p != 0
    """
    s = 0
    if p > 2:
        while a:
            q, r = divmod(a, p)
            if r != 0:
                break
            s += 1
            a = q
    elif p == 2:
        while a and a & 1 == 0:
            s += 1
            a >>= 1
    else:
        raise ValueError("Number %d is smaller than 2" % p)
    return s, a


def extract_prime_power(a, p):
    return extract_power(a, p)


def solve_linear(a, b, c):
    """
    Solve a*x + b*y = c.
    Solution (x0 + b*n, y0 + a*n).
    Return None or (x0, y0).
    """
    #TODO: do
    return None
