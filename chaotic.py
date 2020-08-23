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

if __name__ == "__main__":
    t = np.arange(0, 0.2*1000, 0.2)
    # 调用odeint对chen进行求解
    P1 = odeint(chen, (1., 0.5, 1.2), t, args=([36., 3., 28.,16,2],))
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(P1[:, 0], P1[:, 1], P1[:, 2])
    plt.show()