#-*- coding:utf-8 -*-

import random

from functools import reduce
from libnum.factorize import factorize
from libnum.factorize import is_power
from utcompat import *


def test_powers():
    numbers = [2, 3, 5, 6, 7, 10, 1993, 1995]
    pows = [2, 3, 4, 5, 6, 7, 8, 10, 30]
    for n in numbers:
        for k in pows:
            val = n ** k
            assertEqual(is_power(val), (n, k))


def test_not_powers():
    numbers = [2, 3, 5, 6, 7, 10, 1993, 1995]
    for n in numbers:
        assertFalse(is_power(n))


def test_samples():
    test_lists = [
        set([(3, 1), (5, 5), (19, 2), (1993, 1), (37, 1), (2, 4)]),
        set([(2, 100), (3, 50), (5, 20), (7, 15), (11, 10), (13, 5)]),
        set([(1993, 5)]),
        set([(2, 4000)]),
    ]

    for primes_list in test_lists:
        while primes_list:
            n = reduce(lambda a, b: a * (b[0] ** b[1]), primes_list, 1)
            primes_list_test = set(sorted(factorize(n).items()))
            assertEqual(primes_list, primes_list_test)
            primes_list.pop()


def test_zero():
    assertEqual(factorize(0), {0: 1})


def test_small():
    assertEqual(factorize(1), {})
    assertEqual(factorize(2), {2: 1})
    assertEqual(factorize(3), {3: 1})
    assertEqual(factorize(4), {2: 2})
    assertEqual(factorize(5), {5: 1})
    assertEqual(factorize(6), {2: 1, 3: 1})
    assertEqual(factorize(7), {7: 1})
    assertEqual(factorize(8), {2: 3})
    assertEqual(factorize(9), {3: 2})
    assertEqual(factorize(10), {2: 1, 5: 1})


pvals = [2, 3, 5, 7, 11, 13, 17, 23, 31, 41, 53, 67, 83, 107, 139, 179, 227,
         277, 353, 439, 557, 673, 829, 1021, 1249, 1523, 1861, 2243, 2689,
         3221, 3821, 4517, 5323, 6203, 7187, 8287, 9437, 10723, 12143, 13591,
         15131, 16703, 18371, 20161, 21961, 23767, 25693, 27689, 29573, 31627,
         33587, 35597, 37643, 39821, 41903, 43997, 46181, 48371, 50423, 52631,
         54787, 56929, 59083, 61363, 63599, 65827, 67993, 70379, 72533, 74771,
         77069, 79333, 81533, 83777, 86137, 88469, 90647, 92867, 95177, 97463,
         99793, 102061, 104473]


def test_prime_powers():
    print("")
    for i in range(len(pvals)):
        print("Testing prime power % 3d/%d" % (i, len(pvals)))
        p = pvals[i]
        for e in (list(range(1, 10)) + [10, 20, 50, 100]):
            assertEqual(factorize(p ** e), {p: e})


def test_composite_powers():
    count = 20
    print("")
    for i in range(count):
        print("Testing composite power % 3d/%d" % (i, count))
        p = list(random.sample(pvals, random.randint(2, 5)))
        k = [random.randint(1, 3) for _ in range(len(p))]
        n = reduce(lambda x, y: x * y, (p[i] ** k[i]
                                        for i in range(len(p))))
        for e in (1, 2, 5, 10, 20, 50):
            assertEqual(factorize(n ** e), {p[i]: k[i] * e
                                            for i in range(len(p))})


def test_small_negative():
    assertEqual(factorize(-1), {-1: 1})
    assertEqual(factorize(-2), {-1: 1, 2: 1})
    assertEqual(factorize(-3), {-1: 1, 3: 1})
    assertEqual(factorize(-4), {-1: 1, 2: 2})
    assertEqual(factorize(-5), {-1: 1, 5: 1})
    assertEqual(factorize(-6), {-1: 1, 2: 1, 3: 1})
    assertEqual(factorize(-7), {-1: 1, 7: 1})
    assertEqual(factorize(-8), {-1: 1, 2: 3})
    assertEqual(factorize(-9), {-1: 1, 3: 2})
    assertEqual(factorize(-10), {-1: 1, 2: 1, 5: 1})


def test_errors():
    assertRaises(TypeError, factorize, "1")
    assertRaises(TypeError, factorize, 10.3)
    assertRaises(TypeError, factorize, complex(10, 3))
    assertRaises(TypeError, factorize, (2, 3))
