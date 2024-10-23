# Trailer path planning library
# author: Atsushi Sakai(@Atsushi_twi)

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

# Vehicle parameters
WB = 3.7  # [m] wheel base: rear to front steer
LT = 8.0  # [m] rear to trailer wheel
W = 2.6   # [m] width of vehicle
LF = 4.5  # [m] distance from rear to vehicle front end of vehicle
LB = 1.0  # [m] distance from rear to vehicle back end of vehicle
LTF = 1.0  # [m] distance from rear to vehicle front end of trailer
LTB = 9.0  # [m] distance from rear to vehicle back end of trailer
MAX_STEER = 0.6  # [rad] maximum steering angle 
TR = 0.5  # Tyre radius [m] for plot
TW = 1.0  # Tyre width [m] for plot

# For collision check
WBUBBLE_DIST = 3.5  # Distance from rear to center of whole bubble
WBUBBLE_R = 10.0  # Whole bubble radius
B = 4.45  # Distance from rear to vehicle back end
C = 11.54  # Distance from rear to vehicle front end
I = 8.55  # Width of vehicle
VRX = np.array([C, C, -B, -B, C])
VRY = np.array([-I / 2.0, I / 2.0, I / 2.0, -I / 2.0, -I / 2.0])


def check_collision(x, y, yaw, kdtree, ox, oy, wbd, wbr, vrx, vry):
    for ix, iy, iyaw in zip(x, y, yaw):
        cx = ix + wbd * np.cos(iyaw)
        cy = iy + wbd * np.sin(iyaw)

        # Whole bubble check
        ids = kdtree.query_ball_point([cx, cy], wbr)
        if len(ids) == 0:
            continue

        # if not rect_check(ix, iy, iyaw, ox[ids], oy[ids], vrx, vry):
        #     return False  # Collision

    return True  # OK


def rect_check(ix, iy, iyaw, ox, oy, vrx, vry):
    c = np.cos(-iyaw)
    s = np.sin(-iyaw)

    for iox, ioy in zip(ox, oy):
        tx = iox - ix
        ty = ioy - iy
        lx = (c * tx - s * ty)
        ly = (s * tx + c * ty)

        sumangle = 0.0
        for i in range(len(vrx) - 1):
            x1 = vrx[i] - lx
            y1 = vry[i] - ly
            x2 = vrx[i + 1] - lx
            y2 = vry[i + 1] - ly
            d1 = np.hypot(x1, y1)
            d2 = np.hypot(x2, y2)
            theta1 = np.arctan2(y1, x1)
            tty = (-np.sin(theta1) * x2 + np.cos(theta1) * y2)
            tmp = (x1 * x2 + y1 * y2) / (d1 * d2)

            tmp = np.clip(tmp, -1.0, 1.0)

            if tty >= 0.0:
                sumangle += np.arccos(tmp)
            else:
                sumangle -= np.arccos(tmp)

        if abs(sumangle) >= np.pi:
            return False  # Collision

    return True  # OK


def calc_trailer_yaw_from_xyyaw(x, y, yaw, init_tyaw, steps):
    """
    Calculate trailer yaw from x, y, yaw lists
    """
    tyaw = np.zeros(len(x))
    tyaw[0] = init_tyaw

    for i in range(1, len(x)):
        tyaw[i] += tyaw[i - 1] + steps[i - 1] / LT * np.sin(yaw[i - 1] - tyaw[i - 1])

    return tyaw


def trailer_motion_model(x, y, yaw0, yaw1, D, d, L, delta):
    """
    Motion model for trailer 
    see: http://planning.cs.uiuc.edu/node661.html#77556
    """
    x += D * np.cos(yaw0)
    y += D * np.sin(yaw0)
    yaw0 += D / L * np.tan(delta)
    yaw1 += D / d * np.sin(yaw0 - yaw1)

    return x, y, yaw0, yaw1


def check_trailer_collision(ox, oy, x, y, yaw0, yaw1, kdtree=None):
    """
    Collision check function for trailer
    """
    if kdtree is None:
        kdtree = KDTree(np.vstack((ox, oy)).T)

    vrxt = np.array([LTF, LTF, -LTB, -LTB, LTF])
    vryt = np.array([-W / 2.0, W / 2.0, W / 2.0, -W / 2.0, -W / 2.0])

    # Bubble parameters
    DT = (LTF + LTB) / 2.0 - LTB
    DTR = (LTF + LTB) / 2.0 + 0.3

    # Check trailer
    if not check_collision(x, y, yaw1, kdtree, ox, oy, DT, DTR, vrxt, vryt):
        return False

    vrxf = np.array([LF, LF, -LB, -LB, LF])
    vryf = np.array([-W / 2.0, W / 2.0, W / 2.0, -W / 2.0, -W / 2.0])

    # Bubble parameters
    DF = (LF + LB) / 2.0 - LB
    DFR = (LF + LB) / 2.0 + 0.3

    # Check front trailer
    if not check_collision(x, y, yaw0, kdtree, ox, oy, DF, DFR, vrxf, vryf):
        return False

    return True  # OK


def plot_trailer(x, y, yaw, yaw1, steer):
    truckcolor = "-k"

    LENGTH = LB + LF
    LENGTHt = LTB + LTF

    truckOutLine = np.array([[-LB, (LENGTH - LB), (LENGTH - LB), -LB, -LB],
                             [W / 2, W / 2, -W / 2, -W / 2, W / 2]])

    trailerOutLine = np.array([[-LTB, (LENGTHt - LTB), (LENGTHt - LTB), -LTB, -LTB],
                               [W / 2, W / 2, -W / 2, -W / 2, W / 2]])

    rr_wheel = np.array([[TR, -TR, -TR, TR, TR],
                         [-W / 12.0 + TW, -W / 12.0 + TW, W / 12.0 + TW, W / 12.0 + TW, -W / 12.0 + TW]])

    rl_wheel = np.array([[TR, -TR, -TR, TR, TR],
                         [-W / 12.0 - TW, -W / 12.0 - TW, W / 12.0 - TW, W / 12.0 - TW, -W / 12.0 - TW]])

    fr_wheel = rr_wheel.copy()
    fl_wheel = rl_wheel.copy()

    tr_wheel = rr_wheel.copy()
    tl_wheel = rl_wheel.copy()

    # Rotation matrices
    Rot1 = np.array([[np.cos(yaw), np.sin(yaw)],
                     [-np.sin(yaw), np.cos(yaw)]])
    Rot2 = np.array([[np.cos(steer), np.sin(steer)],
                     [-np.sin(steer), np.cos(steer)]])
    Rot3 = np.array([[np.cos(yaw1), np.sin(yaw1)],
                     [-np.sin(yaw1), np.cos(yaw1)]])

    # Front wheels transformation
    fr_wheel = Rot2 @ fr_wheel
    fl_wheel = Rot2 @ fl_wheel
    fr_wheel[0, :] += WB
    fl_wheel[0, :] += WB
    fr_wheel = Rot1 @ fr_wheel
    fl_wheel = Rot1 @ fl_wheel

    # Trailer wheels transformation
    tr_wheel[0, :] -= LT
    tl_wheel[0, :] -= LT
    tr_wheel = Rot3 @ tr_wheel
    tl_wheel = Rot3 @ tl_wheel

    # Transform truck and trailer outlines
    truckOutLine = Rot1 @ truckOutLine
    trailerOutLine = Rot3 @ trailerOutLine

    # Plotting
    plt.plot(trailerOutLine[0, :] + x, trailerOutLine[1, :] + y, truckcolor)
    plt.plot(tr_wheel[0, :] + x, tr_wheel[1, :] + y, truckcolor)
    plt.plot(tl_wheel[0, :] + x, tl_wheel[1, :] + y, truckcolor)

    plt.plot(truckOutLine[0, :] + x, truckOutLine[1, :] + y, truckcolor)
    plt.plot(fr_wheel[0, :] + x, fr_wheel[1, :] + y, truckcolor)
    plt.plot(fl_wheel[0, :] + x, fl_wheel[1, :] + y, truckcolor)
    plt.plot(rr_wheel[0, :] + x, rr_wheel[1, :] + y, truckcolor)
    plt.plot(rl_wheel[0, :] + x, rl_wheel[1, :] + y, truckcolor)


if __name__ == "__main__":
    # Test example
    ox = [i for i in range(50)]
    oy = [i for i in range(50)]
    
    # Create kd-tree for collision checking
    kdtree = KDTree(np.vstack((ox, oy)).T)

    # Test for collision check function
    x = np.array([0.0])
    y = np.array([0.0])
    yaw = np.array([0.0])
    yaw1 = np.array([0.0])

    is_collision = check_trailer_collision(ox, oy, x, y, yaw, yaw1, kdtree)
    print(f"Is collision: {is_collision}")

    # Test plot
    plt.figure()
    plot_trailer(0.0, 0.0, 0.0, 0.0, 0.0)
    plt.axis("equal")
    plt.show()
