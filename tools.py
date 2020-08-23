# /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper Doi:10.1080/19393555.2020.1718248

This module is tools used in a
selective digital image cryptosystem
using DNA encoding and dul dual hyperchaos map.

1. dec_bin: convert decimal matrix into binary matrix
2. bin_dec: convert binary matrix into decimal matrix
3. dna_encoding:
convert binary matrix into DNA-encoded matrix
convert DNA-encoded matrix into binary matrix
"""

__author__ = "Yuming Chen"
__copyright__ = "Copyright 2020, Study Project in Lanzhou University , China"
__license__ = "GPL V3"
__maintainer__ = "Yuming Chen"
__email__ = ["chenym18@lzu.edu.cn"]
__status__ = "Experimental"

from copy import deepcopy
import numpy as np
import json

RULE_PATH = 'rules.json'

def dec_bin(matrix):
    width, height = matrix.shape
    copy = np.zeros(matrix.shape).astype(str)
    for x in range(width):
        for y in range(height):
            bin_ = '{:08b}'.format(matrix[x, y])
            copy[x, y] = bin_
    return copy


def bin_dec(matrix):
    width, height = matrix.shape
    copy = deepcopy(matrix).astype(np.uint8)
    for x in range(width):
        for y in range(height):
            dec = int(matrix[x, y], 2)
            copy[x, y] = dec
    return copy


def dna_encoding(m, encode=True):
    with open(RULE_PATH, 'r') as f:
        rules = json.load(f).get(str(encode))

    width, height = m.shape
    copy = deepcopy(m).astype(str)
    index = 0
    for x in range(width):
        for y in range(height):
            rule = rules.get(str(np.mod(index, 8)))
            dna_code = ''
            for _ in range(0, len(m[x, y]), 2):
                dna_code += rule.get(m[x, y][_:_ + 2])
            copy[x, y] = dna_code
            index += 1
    return copy


if __name__ == "__main__":
    a = {'00': 'A', '01': 'T', '10': 'G', '11': 'C'}
    b = {'A': '00', 'T': '01', 'G': '10', 'C': '11'}
    with open(RULE_PATH, 'r') as f:
        rules = json.load(f).get(str(False))
    for _ in rules:
        print('rule_{}'.format(_),end=' ')
        for key in rules[_]:
            print(str(a[rules[_][key]]) + '=' + str(key), end=' ')
        print('')
    with open(RULE_PATH, 'r') as f:
        rules = json.load(f).get(str(True))
    print('----------------------')
    for _ in rules:
        for key in rules[_]:
            print(str(rules[_][key]) + '=' + str(a[key]), end=' ')
        print('')
