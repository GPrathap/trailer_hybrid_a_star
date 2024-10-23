import numpy as np
from queue import PriorityQueue
from scipy.spatial import KDTree
import math

XY_GRID_RESOLUTION = 2.0  # [m]
YAW_GRID_RESOLUTION = np.deg2rad(15.0)  # [rad]
GOAL_TYAW_TH = np.deg2rad(5.0)  # [rad]
MOTION_RESOLUTION = 0.1  # [m] path interpolate resolution
N_STEER = 20.0  # number of steer command
EXTEND_AREA = 5.0  # [m] map extend length
SKIP_COLLISION_CHECK = 4  # skip number for collision check

SB_COST = 100.0  # switch back penalty cost
BACK_COST = 5.0  # backward penalty cost
STEER_CHANGE_COST = 5.0  # steer angle change penalty cost
STEER_COST = 1.0  # steer angle change penalty cost
JACKKNIF_COST = 200.0  # Jackknif cost
H_COST = 5.0  # Heuristic cost

# Assuming `trailerlib` is some utility library that provides constants WB, LT, MAX_STEER, etc.
WB = trailerlib.WB  # [m] Wheel base
LT = trailerlib.LT  # [m] length of trailer
MAX_STEER = trailerlib.MAX_STEER  # [rad] maximum steering angle


class Node:
    def __init__(self, xind, yind, yawind, direction, x, y, yaw, yaw1, directions, steer, cost, pind):
        self.xind = xind  # x index
        self.yind = yind  # y index
        self.yawind = yawind  # yaw index
        self.direction = direction  # moving direction forward: true, backward: false
        self.x = x  # x position [m]
        self.y = y  # y position [m]
        self.yaw = yaw  # yaw angle [rad]
        self.yaw1 = yaw1  # trailer yaw angle [rad]
        self.directions = directions  # directions of each point, forward: true, backward: false
        self.steer = steer  # steer input
        self.cost = cost  # cost
        self.pind = pind  # parent index


class Config:
    def __init__(self, minx, miny, minyaw, minyawt, maxx, maxy, maxyaw, maxyawt, xw, yw, yaww, yawtw, xyreso, yawreso):
        self.minx = minx
        self.miny = miny
        self.minyaw = minyaw
        self.minyawt = minyawt
        self.maxx = maxx
        self.maxy = maxy
        self.maxyaw = maxyaw
        self.maxyawt = maxyawt
        self.xw = xw
        self.yw = yw
        self.yaww = yaww
        self.yawtw = yawtw
        self.xyreso = xyreso
        self.yawreso = yawreso


class Path:
    def __init__(self, x, y, yaw, yaw1, direction, cost):
        self.x = x  # x position [m]
        self.y = y  # y position [m]
        self.yaw = yaw  # yaw angle [rad]
        self.yaw1 = yaw1  # trailer angle [rad]
        self.direction = direction  # direction forward: true, back: false
        self.cost = cost  # cost


def is_same_grid(node1, node2):
    if node1.xind != node2.xind:
        return False
    if node1.yind != node2.yind:
        return False
    if node1.yawind != node2.yawind:
        return False
    return True


def calc_index(node, c):
    ind = (node.yawind - c.minyaw) * c.xw * c.yw + (node.yind - c.miny) * c.xw + (node.xind - c.minx)

    # 4D grid
    yaw1ind = round(node.yaw1[-1] / c.yawreso)
    ind += (yaw1ind - c.minyawt) * c.xw * c.yw * c.yaww

    if ind <= 0:
        print("Error(calc_index):", ind)
    return ind


def get_final_path(closed, ngoal, nstart, c):
    rx = np.array(list(reversed(ngoal.x)))
    ry = np.array(list(reversed(ngoal.y)))
    ryaw = np.array(list(reversed(ngoal.yaw)))
    ryaw1 = np.array(list(reversed(ngoal.yaw1)))
    direction = np.array(list(reversed(ngoal.directions)))
    nid = ngoal.pind
    finalcost = ngoal.cost

    while True:
        n = closed[nid]
        rx = np.concatenate((rx, np.array(list(reversed(n.x)))))
        ry = np.concatenate((ry, np.array(list(reversed(n.y)))))
        ryaw = np.concatenate((ryaw, np.array(list(reversed(n.yaw)))))
        ryaw1 = np.concatenate((ryaw1, np.array(list(reversed(n.yaw1)))))
        direction = np.concatenate((direction, np.array(list(reversed(n.directions)))))
        nid = n.pind
        if is_same_grid(n, nstart):
            break

    rx = np.array(list(reversed(rx)))
    ry = np.array(list(reversed(ry)))
    ryaw = np.array(list(reversed(ryaw)))
    ryaw1 = np.array(list(reversed(ryaw1)))
    direction = np.array(list(reversed(direction)))

    return Path(rx, ry, ryaw, ryaw1, direction, finalcost)


def calc_holonomic_with_obstacle_heuristic(gnode, ox, oy, xyreso):
    h_dp = grid_a_star.calc_dist_policy(gnode.x[-1], gnode.y[-1], ox, oy, xyreso, 1.0)
    return h_dp


def calc_config(ox, oy, xyreso, yawreso):
    min_x_m = min(ox) - EXTEND_AREA
    min_y_m = min(oy) - EXTEND_AREA
    max_x_m = max(ox) + EXTEND_AREA
    max_y_m = max(oy) + EXTEND_AREA

    ox = np.append(ox, [min_x_m, max_x_m])
    oy = np.append(oy, [min_y_m, max_y_m])

    minx = round(min_x_m / xyreso)
    miny = round(min_y_m / xyreso)
    maxx = round(max_x_m / xyreso)
    maxy = round(max_y_m / xyreso)

    xw = round(maxx - minx)
    yw = round(maxy - miny)

    minyaw = round(-math.pi / yawreso) - 1
    maxyaw = round(math.pi / yawreso)
    yaww = round(maxyaw - minyaw)

    minyawt = minyaw
    maxyawt = maxyaw
    yawtw = yaww

    config = Config(minx, miny, minyaw, minyawt, maxx, maxy, maxyaw, maxyawt, xw, yw, yaww, yawtw, xyreso, yawreso)

    return config


# Assuming the trailerlib, grid_a_star, and rs_path functions are imported and work in Python similar to Julia
