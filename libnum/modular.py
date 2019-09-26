#-*- coding:utf-8 -*-

import operator

from functools import reduce
from .compat import builtins, xrange
from .common import *
from .stuff import *
from .factorize import factorize


def has_invmod(a, m):
    """
    Check if a has an inverse modulo m.
    """
    if m < 2:
        raise ValueError("modulus must be greater than 1")

    # This statement along with the above check covers all edge cases
    return _gcd(a, m) == 1


def invmod(a, m):
    """
    Return a**-1 (mod m).
    a and m must be coprime integers.
    """
    if m < 2:
        raise ValueError("modulus must be greater than 1")

    x, y, g = xgcd(a, m)

    if g != 1:
        raise ValueError("%d has no inverse mod %d" % (a, m))
    else:
        return x % m


def pow(a, b, m=None):
    if b >= 0 or m is None:
        return builtins.pow(a, b, m)
    return builtins.pow(invmod(a, m), b, m)


def solve_crt(remainders, moduli):
    """
    Solve Chinese Remainder Theorem.
    All moduli must be pairwise coprime.
    Return x : all((x - remainders[i]) % moduli[i] == 0
                   for i in range(len(moduli))).
    """
    if len(moduli) != len(remainders):
        raise ValueError("Different number of remainders and moduli")

    if len(moduli) == 0:
        raise ValueError("No moduli were given")

    if len(moduli) == 1:
        return remainders[0] % moduli[0]

    x = 0
    N = reduce(operator.mul, moduli)
    for i, modulus in enumerate(moduli):
        if modulus == 1:
            continue

        Ni = N // modulus
        b = invmod(Ni, modulus)

        x += remainders[i] * Ni * b
    return x % N


def nCk_mod(n, k, m, factors=None):
    """
    Compute n choose k mod m.
    The factorization of m may be given or will be computed.
    """
    if factors is None:
        try:
            factors = m.items()
        except AttributeError:
            factors = factorize(m).items()
    rems = []
    mods = []
    for p, e in factors:
        rems.append(nCk_mod_prime_power(n, k, p, e))
        mods.append(p ** e)
    return solve_crt(rems, mods)


def factorial_mod(n, m, factors=None):
    """
    Compute n! mod m.
    The factorization of m may be given or will be computed.
    """
    if factors is None:
        try:
            factors = m.items()
        except AttributeError:
            factors = factorize(m).items()
    rems = []
    mods = []
    for p, e in factors:
        pe = p ** e
        if n >= pe or factorial_get_prime_pow(n, p) >= e:
            factmod = 0
        else:
            factmod = factorial(n) % pe
        rems.append(factmod)
        mods.append(pe)
    return solve_crt(rems, mods)


def nCk_mod_prime_power(n, k, p, e):
    """
    Compute n choose k mod p**e
    Algorithm by Andrew Granville:
        http://www.dms.umontreal.ca/~andrew/PDF/BinCoeff.pdf
    What can be optimized:
        - compute (n-k)*(n-k+1)*...*n / 1*2*...*k instead of n!, k!, r!
        - ...
    """

    def nCk_get_prime_pow(n, k, p):
        res = factorial_get_prime_pow(n, p)
        res -= factorial_get_prime_pow(k, p)
        res -= factorial_get_prime_pow(n - k, p)
        return res

    def nCk_get_non_prime_part(n, k, p, e):
        pe = p ** e
        r = n - k

        fact_pe = [1]
        acc = 1
        for x in xrange(1, pe):
            if x % p == 0:
                x = 1
            acc = (acc * x) % pe
            fact_pe.append(acc)

        top = bottom = 1
        is_negative = 0
        digits = 0

        while n != 0:
            if acc != 1:
                if digits >= e:
                    is_negative ^= n & 1
                    is_negative ^= r & 1
                    is_negative ^= k & 1

            top = (top * fact_pe[n % pe]) % pe
            bottom = (bottom * fact_pe[r % pe]) % pe
            bottom = (bottom * fact_pe[k % pe]) % pe

            n //= p
            r //= p
            k //= p

            digits += 1

        res = (top * invmod(bottom, pe)) % pe
        if p != 2 or e < 3:
            if is_negative:
                res = pe - res
        return res

    prime_part_pow = nCk_get_prime_pow(n, k, p)
    if prime_part_pow >= e:
        return 0

    modpow = e - prime_part_pow

    r = nCk_get_non_prime_part(n, k, p, modpow) % (p ** modpow)
    return ((p ** prime_part_pow) * r) % (p ** e)
