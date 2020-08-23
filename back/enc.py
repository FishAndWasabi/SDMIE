# /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module is encryption program for a
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
import numpy as np
import cv2
import warnings
warnings.filterwarnings('ignore')

class DCE(object):
    def __init__(self, x0: float, y0: float, k: float, a: float = 36,
                 b: float = 3, c: float = 28, d: float = 16, l: float = 0.2):
        self.__x0 = x0
        self.__y0 = y0
        self.__K = k
        self.__a = a
        self.__b = b
        self.__c = c
        self.__d = d
        self.__l = l
        self.__image = None
        self.__map = None
        self.__key = None

    @property
    def image(self):
        assert not isinstance(self.__image,bool), "Load Image First"
        return self.__image

    def load_image(self, path):
        self.__image = cv2.imread(path, cv2.IMREAD_GRAYSCALE)

    @property
    def key(self):
        assert not isinstance(self.__key,bool), "Encrypt First"
        return self.__key

    @property
    def map(self):
        assert not isinstance(self.__map,bool), "Encrypt First"
        return self.__map

    def encrypt(self):
        """
        step1
        """
        assert not isinstance(self.__image, bool), "Load Image First"
        if not self.__map:
            self._gen_map()
        m0, m1 = self._pixel_selection()
        b0, b1 = dec_bin(m0), dec_bin(m1)
        c0, c1 = dna_encoding(b0), dna_encoding(b1)
        map_x, map_y = self.__map
        self.__key = deepcopy(c0)
        for x0, x1 in zip(range(self.image.shape[0]), map_x):
            for y0, y1 in zip(range(self.image.shape[1]), map_y):
                self.__key[x0, y0] = c0[x1, y1]
        self.__key, c1 = bin_dec(self.__key), bin_dec(c1)
        self.__image = self.__key ^ c1
        self.__image = dec_bin(self.__image)
        self.__image = dna_encoding(self.__image, encode=False)
        self.__image = bin_dec(self.__image)
        return self.image

    def _pixel_selection(self):
        width, height = self.image.shape
        m0, m1 = np.zeros(self.image.shape).astype(np.uint8), np.zeros(self.image.shape).astype(np.uint8)
        for x in range(width):
            for y in range(height):
                a = self.image[x, y] / 3
                b = np.floor(a)
                if b < a:
                    m0[x, y] = self.image[x, y]
                else:
                    m1[x, y] = self.image[x, y]
        return m0, m1

    # def _gen_map(self):
    #     """
    #     """
    #     xn, yn, z, q = self.__x0, self.__y0, self.__k + 0.4, self.__x0 + self.__l
    #     map_x, map_y = np.zeros(max(self.image.shape), dtype=float), \
    #                    np.zeros(max(self.image.shape), dtype=float)
    #     for _ in range(max(self.image.shape)):
    #         x, y, z, xn, yn, q = self._equation(xn, yn, z, q)
    #         map_x[_], map_y[_] = x, y
    #     map_x, map_y = np.argsort(map_x[:self.image.shape[0]]), np.argsort(map_y[:self.image.shape[1]])
    #     self.__map = map_x, map_y
    #     return
    #
    # def _equation(self, xn, yn, zn, q):
    #     """
    #     Equation to generate the Taylor-Chirikov map and Chen's hyper chaotic map
    #     The taylor chirikov map provides initial values
    #     """
    #     # Taylor Chirikov map
    #     xn1 = (xn + self.__k * np.math.sin(yn)) / 2
    #     yn1 = (yn + np.mod(xn1, 2 * np.math.pi)) - 0.4
    #     # Chen's hyper chaotic map
    #     x = self.__a * (yn1 - xn1)
    #     y = (-xn1 * zn) + self.__d * xn1 + self.__c * yn1 - q
    #     z = xn1 * yn1 - self.__b * zn
    #     q = xn1 + self.__l
    #     return x, y, z, xn1, yn1, q

    def _gen_map(self):
        """
        Generate the Taylor-Chirikov map and Chen's hyper chaotic map
        1) The Taylor-Chirikov map provides initial values for system parameters of Chen's hyper chaotic map.
        2) The chaotic sequence of Chen's shuffles the pixels of the digital image.
        """
        # Taylor-Chirikov map
        x, y, z, q = self.__x0, self.__y0, self.__K + 0.4, self.__x0 + self.__l
        x = (x + self.__K * np.math.sin(y)) / 2
        y = (y + np.mod(x, 2 * np.math.pi)) - 0.4
        # Chen's hyper chaotic map
        map_x, map_y = np.zeros(max(self.image.shape), dtype=float), \
                       np.zeros(max(self.image.shape), dtype=float)
        for _ in range(max(self.image.shape)):
            x = self.__a * (y - x)
            y = (-x * z) + self.__d * x + self.__c * y - q
            z = x * y - self.__b * z
            q = x + self.__l
            map_x[_], map_y[_] = x, y
        map_x, map_y = np.argsort(map_x[:self.image.shape[0]]), np.argsort(map_y[:self.image.shape[1]])
        self.__map = map_x, map_y
        return

if __name__ == "__main__":
    dce = DCE(x0=1, y0=0.5, k=0.8)
    dce.load_image('test.png')
    result = dce.encrypt()
    cv2.imshow('test_e',dce.key)
    cv2.waitKey()
