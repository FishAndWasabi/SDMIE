# /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Paper Doi:10.1080/19393555.2020.1718248

This module is used to generate chen's hyperchaotic map
"""

import numpy as np
from scipy.integrate import odeint


def chen(points, t, params):
    """
    p：位置矢量
    sets：其他参数
    """
    a, b, c, d, xn = params
    x, y, z = points
    return np.array([a*(y - x), -x * z + d * x + c * y - xn, x * y - b * z])

def chen_(points, t, params):
    """
    p：初始值
    sets：其他参数
    """
    a, b, c = params
    x, y, z = points
    return np.array([a*(y - x), -x * z + (c-a) * x + c * y , x * y - b * z])

if __name__ == "__main__":
    x, y, k,z = 1, 0.5, 0.8,0.8 + 0.4
    xn = x
    for _ in range(500):
        xn = x
        x = (x + k * np.math.sin(y)) / 2
        y = (y + np.mod(x, 2 * np.math.pi)) - 0.4
    t1 = np.arange(0, 0.2*30000, 0.0001)
    t2 = np.arange(0, 30, 0.0001)
    P1 = odeint(chen, (x, y, z), t2, args=([36., 3., 28.,16,2],))
    P2 = odeint(chen_, (x, y, z), t2, args=([36., 3., 28.],))

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(P1[:, 0], P1[:, 1], P1[:, 2])
    plt.show()
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(P2[:, 0], P2[:, 1], P2[:, 2])
    plt.show()