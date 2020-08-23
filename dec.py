# /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper Doi:10.1080/19393555.2020.1718248

This module is decryption program for a
selective digital image cryptosystem
using DNA encoding and dul dual hyperchaos map.
"""

__author__ = "Yuming Chen"
__copyright__ = "Copyright 2020, Study Project in Lanzhou University , China"
__license__ = "GPL V3"
__maintainer__ = "Yuming Chen"
__email__ = ["chenym18@lzu.edu.cn"]
__status__ = "Experimental"

from tools import bin_dec, dec_bin, dna_encoding
from copy import deepcopy
from scipy.integrate import odeint
from chaotic import chen
import numpy as np
import warnings
import base64
import cv2
import ast
import os

warnings.filterwarnings('ignore')


class DCD(object):
    def __init__(self, path):
        assert os.path.exists(path), "Key doesn't exit"
        with open(path, 'rb') as f:
            self.__key = ast.literal_eval(base64.b64decode(f.read()).decode())
        self.__x0, self.__y0, self.__k, self.__a, self.__b, self.__c, self.__d, self.__l,self.__n = self.__key['params']
        self.__image = None
        self.__map = None

    @property
    def image(self):
        assert isinstance(self.__image, np.ndarray), "Load Image First"
        return self.__image

    def load_image(self, path):
        assert os.path.exists(path), "Image doesn't exit"
        self.__image = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        assert self.__image.shape > (0, 0), "Image doesn't exit"
        return

    @property
    def key(self):
        return self.__key

    @property
    def map(self):
        assert isinstance(self.__map, np.ndarray), "Encrypt First"
        return self.__map

    def decrypt(self):
        """
        The inverse of encryption
        """
        if not self.__map:
            self._gen_map()
        map_x, map_y = self.__map
        self.__key['array'] = np.array(self.__key['array'], dtype=np.uint8)
        self.__image = dec_bin(self.__image)
        self.__image = dna_encoding(self.__image, encode=True)
        self.__image = bin_dec(self.__image)
        c1 = self.__image ^ self.__key['array']
        c0, c1 = dec_bin(self.__key['array']), dec_bin(c1)
        c0_ = deepcopy(c0)
        for x0, x1 in zip(range(self.image.shape[0]), map_x):
            for y0, y1 in zip(range(self.image.shape[1]), map_y):
                c0[x1, y1] = c0_[x0, y0]
        b0, b1 = dna_encoding(c0, False), dna_encoding(c1, False)
        self.__image = bin_dec(b0) | bin_dec(b1)
        return self.image

    def _gen_map(self):
        """
        1. describe
        Generate the Taylor-Chirikov map and Chen's hyper chaotic map
        Use the index of sorted chaotic map to shuffle the pixel of the digital image
        1) The Taylor-Chirikov map provides initial values for system parameters of Chen's hyper chaotic map.
        2) The chaotic sequence of Chen's shuffles the pixels of the digital image.
        """
        x, y, z = self.__x0, self.__y0, self.__k + 0.4
        # Taylor-Chirikov map
        xn = x
        for _ in range(self.__n):
            xn = x
            x = (x + self.__k * np.math.sin(y)) / 2
            y = (y + np.mod(x, 2 * np.math.pi)) - 0.4
        # Chen's hyper chaotic map
        t = np.arange(x, self.__l * max(self.image.shape) - x, self.__l)
        # Use odeint to solve chen
        map_ = odeint(chen, (x,y,z), t, args=([self.__a,self.__b,self.__c,self.__d,xn],))
        # Get the index of sorted sequence which is used to suffle the digital image
        map_x, map_y = np.argsort(map_[:self.image.shape[0], 0]), np.argsort(map_[:self.image.shape[1], 1])
        self.__map = map_x,map_y
        return

if __name__ == "__main__":
    dcd = DCD('test\\dec\\key')
    dcd.load_image('test\\dec\\test_e.png')
    result = dcd.decrypt()
    cv2.imshow('after decrypt', dcd.image)
    cv2.waitKey()
    cv2.imwrite('test\\dec\\test_d.png',result)
