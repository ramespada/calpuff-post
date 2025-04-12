#!/bin/python3

from struct   import unpack

def _decompress(xwork):
    xdat = []
    for value in xwork:
        if value > 0.0:
            xdat.append(value)
        else:
            xdat.extend([0.0]*int(-value))
    return xdat

def _skip_n_lines(f,n):
    for _ in range(n):
        f.seek(unpack("i", f.read(4))[0] + 4, 1)
