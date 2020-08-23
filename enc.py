# /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper Doi:10.1080/19393555.2020.1718248

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
from scipy.integrate import odeint
from chaotic import chen
import numpy as np
import warnings
import base64
import cv2
import os
import json
warnings.filterwarnings('ignore')


class DCE(object):
    def __init__(self, x0: float, y0: float, k: float, a: float = 36,
                 b: float = 3, c: float = 28, d: float = 16, l: float = 0.2):
        self.__x0 = x0
        self.__y0 = y0
        self.__k = k
        self.__a = a
        self.__b = b
        self.__c = c
        self.__d = d
        self.__l = l
        self.__image = None
        self.__map = None
        self.__key = {}

    @property
    def image(self):
        assert isinstance(self.__image, np.ndarray), "Load Image First"
        return self.__image

    def load_image(self, path):
        assert os.path.exists(path) , "Image doesn't exit"
        self.__image = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        assert self.__image.shape > (0, 0), "Image doesn't exit"
        return

    @property
    def key(self):
        assert self.__key, "Encrypt First"
        self.__key['params'] = [self.__x0, self.__y0, self.__k, self.__a, self.__b, self.__c, self.__d, self.__l]
        return base64.b64encode(json.dumps(self.__key).encode())

    @property
    def map(self):
        assert isinstance(self.__map, np.ndarray), "Encrypt First"
        return self.__map

    def encrypt(self):
        """
        step1
        Divide into two matrices
        step2
        Convert decimal image into binary images
        step3
        Convert  8-bit binary image into DNA-encoded matrix using DNA base encoding rules
        step4
        Generate the chaotic sequences x,y using dual hyperchaos map
        step5
        Jumble the pixels of image by the index value of sorted chaotic sequences
        step6
        Using DNA operation XOR to fusion the matrices
        step7
        Convert into DNA-encoded matrix into 8-bit binary image using DNA base decoding rules
        step8
        Convert binary image into cipher image
        """
        assert isinstance(self.__image, np.ndarray), "Load Image First"
        if not self.__map:
            self._gen_map()
        # step1
        m0, m1 = self._pixel_selection()
        # step2
        b0, b1 = dec_bin(m0), dec_bin(m1)
        # step3
        c0, c1 = dna_encoding(b0), dna_encoding(b1)
        # step4
        map_x, map_y = self.__map
        # step5
        self.__key['array'] = deepcopy(c0)
        for x0, x1 in zip(range(self.image.shape[0]), map_x):
            for y0, y1 in zip(range(self.image.shape[1]), map_y):
                self.__key['array'][x0, y0] = c0[x1, y1]
        # step6
        self.__key['array'], c1 = bin_dec(self.__key['array']), bin_dec(c1)
        self.__image = self.__key['array'] ^ c1
        self.__image = dec_bin(self.__image)
        # step7
        self.__image = dna_encoding(self.__image, encode=False)
        # step8
        self.__image = bin_dec(self.__image)
        self.__key['array'] = self.__key['array'].tolist()
        return self.image

    def _pixel_selection(self):
        """
        1. describe
        The pixels are selected from the original digitized medical image
        :return:
        m0: selected pixels of image.
        m1: remaining pixels of image
        """
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
        for _ in range(500):
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
    dce = DCE(x0=1, y0=0.5, k=0.8)
    dce.load_image('test\\enc\\test.png')
    result = dce.encrypt()
    with open('test\\dec\\key', 'wb') as f:
        f.write(dce.key)
    cv2.imwrite('test\\dec\\test_e.png',result)
