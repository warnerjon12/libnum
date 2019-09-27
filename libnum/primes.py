#-*- coding:utf-8 -*-

import random
import operator

from functools import reduce
from bisect import bisect
from .common import _gcd, len_in_bits, extract_prime_power, randint_bits
from .compat import xrange
from .strings import *

_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
_small_primes_product = 1
_primes_bits = [[] for i in range(11)]
_primes_mask = []


def _init():
    global _small_primes_product, _primes, _primes_bits, _primes_mask
    for p in primes(1024, setglobal=True):
        _primes_bits[len_in_bits(p)].append(p)
    _small_primes_product = reduce(operator.mul, _primes)
    _primes_mask = [(x in _primes) for x in xrange(_primes[-1] + 1)]
    primes(100000, setglobal=True)
    return


def primes(max, k=25, setglobal=False):
    """
    Return a list of primes not greater than max.
    Optional argument k (default 25) determines
    the number of rounds used in primality testing.
    Optional argument setglobal specifies whether or
    not to add new primes to the global variable _primes.
    """
    global _primes, _primes_mask
    try:
        max = int(max)
    except ValueError:
        raise TypeError("%s is not a number" % repr(max))

    if max < 2:
        return []

    if max <= _primes[-1]:
        return _primes[:bisect(_primes, max)]

    if setglobal:
        plist = _primes
    else:
        plist = _primes[:]

    pmod30 = (1, 7, 11, 13, 17, 19, 23, 29)

    m30 = _primes[-1] % 30
    p = _primes[-1] - m30
    i = bisect(pmod30, m30)
    if i == 8:
        p += 30
        i = 0
    while p < max:
        for j in range(i, 8):
            q = p + pmod30[j]
            if prime_test(q, k):
                plist.append(q)
        p += 30
        i = 0
    # Possibly went too far. Also causes the list to be copied.
    return plist[:bisect(plist, max)]


def generate_prime(size, k=25):
    """
    Generate a probable prime with bit-length size.
    Optional argument k (default 25) determines
    the number of rounds used in primality testing.
    """
    if size < 2:
        raise ValueError("No primes smaller than 2 bits!")

    if size <= 10:
        return random.choice(_primes_bits[size])

    while True:
        n = randint_bits(size) | 1  # only odd

        if _gcd(_small_primes_product, n) != 1:
            continue

        if prime_test(n, k):
            return n
    return


def generate_prime_from_string(s, size=None, k=25):
    """
    Generate a probable prime starting with s in its string representation.
    Optional argument size determines the bit-length of the returned value.
    Optional argument k (default 25) determines
    the number of rounds used in primality testing.
    """
    if not size:
        if len(s) > 512:
            size = len(s) * 8 + 32
        else:
            size = len(s) * 8 + 16

    if len(s) * 8 >= size:
        raise ValueError("given size is smaller than string length")

    if size % 8:
        raise ValueError("size must be 8*n")

    extend_len = size - len(s) * 8

    visible_part = s2n(s) << extend_len
    hi = 2 ** extend_len

    while True:
        n = visible_part | random.randint(1, hi) | 1  # only even

        if _gcd(_small_primes_product, n) != 1:
            continue

        if prime_test(n, k):
            return n
    return


def prime_test_fermat(p, k=25):
    """
    Test for primality based on Fermat's Little Theorem.
    Fails on Carmichael numbers.
    """
    if p < 2: return False
    if p <= 3: return True
    if p & 1 == 0: return False

    for j in xrange(k):
        a = random.randint(2, p - 1)
        if _gcd(a, p) != 1:
            return False

        result = pow(a, p - 1, p)
        if result != 1:
            return False
    return True


def prime_test_solovay_strassen(p, k=25):
    """
    Test for primality using Solovay-Strassen.
    Stronger than Fermat's test.
    """
    if p < 2: return False
    if p <= 3: return True
    if p & 1 == 0: return False

    for j in xrange(k):
        a = random.randint(2, p - 1)
        if _gcd(a, p) != 1:
            return False

        result = pow(a, (p - 1) // 2, p)
        if result not in (1, p - 1):
            return False

        if result != jacobi(a, p) % p:
            return False
    return True


def prime_test_miller_rabin(p, k=25):
    """
    Test for primality using Miller-Rabin.
    Stronger than Solovay-Strassen's test.
    """
    if p < 2: return False
    if p <= 3: return True
    if p & 1 == 0: return False

    # p - 1 = 2**s * m
    s, m = extract_prime_power(p - 1, 2)

    for j in range(k):
        a = random.randint(2, p - 2)
        if _gcd(a, p) != 1:
            return False

        b = pow(a, m, p)
        if b in (1, p - 1):
            continue

        for i in range(s):
            b = pow(b, 2, p)

            if b == 1:
                return False

            if b == p - 1:
                # is there one more squaring left to result in 1 ?
                if i < s - 1: break  # good
                else: return False   # bad
        else:
            # result is not 1
            return False
    return True


def jacobi(a, n):
    """
    Compute the Jacobi symbol.
    """
    s = 1
    while True:
        if n < 1:
            raise ValueError("Modulus is too small" +
                             "to compute the Jacobi symbol: " + str(n))
        if n & 1 == 0:
            raise ValueError("The Jacobi symbol is defined" +
                             "only for odd moduli")
        if n == 1:
            return 1

        a %= n
        if a in (0, 1): return s * a

        if a & 1 == 0:
            if n % 8 in (3, 5):
                s = -s
            a >>= 1
            continue

        if a % 4 == 3 and n % 4 == 3:
            s = -s

        a, n = n, a
    return


def legendre(a, p):
    return jacobi(a, p)


prime_test = prime_test_miller_rabin


_init()
