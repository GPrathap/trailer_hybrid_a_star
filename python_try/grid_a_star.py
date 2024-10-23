import matplotlib.pyplot as plt
import numpy as np
import math
from heapq import heappop, heappush

class Node:
    def __init__(self, x, y, cost, pind):
        self.x = x  # x index
        self.y = y  # y index
        self.cost = cost  # cost
        self.pind = pind  # parent index

    def __lt__(self, other):
        return self.cost < other.cost

def calc_dist_policy(gx, gy, ox, oy, reso, vr):
    """
    gx: goal x position [m]
    gy: goal y position [m]
    ox: x position list of Obstacles [m]
    oy: y position list of Obstacles [m]
    reso: grid resolution [m]
    vr: vehicle radius[m]
    """
    ngoal = Node(round(gx / reso), round(gy / reso), 0.0, -1)

    ox = [iox / reso for iox in ox]
    oy = [ioy / reso for ioy in oy]

    obmap, minx, miny, maxx, maxy, xw, yw = calc_obstacle_map(ox, oy, reso, vr)

    open_set, closed_set = {}, {}
    open_set[calc_index(ngoal, xw, minx, miny)] = ngoal

    motion = get_motion_model()
    pq = []
    heappush(pq, (ngoal.cost, calc_index(ngoal, xw, minx, miny)))

    while pq:
        c_id = heappop(pq)[1]
        current = open_set[c_id]

        del open_set[c_id]
        closed_set[c_id] = current

        for i in range(len(motion)):
            node = Node(current.x + motion[i][0], current.y + motion[i][1], current.cost + motion[i][2], c_id)

            if not verify_node(node, minx, miny, xw, yw, obmap):
                continue

            node_ind = calc_index(node, xw, minx, miny)

            if node_ind in closed_set:
                continue

            if node_ind in open_set:
                if open_set[node_ind].cost > node.cost:
                    open_set[node_ind].cost = node.cost
                    open_set[node_ind].pind = c_id
            else:
                open_set[node_ind] = node
                heappush(pq, (node.cost, node_ind))

    pmap = calc_policy_map(closed_set, xw, yw, minx, miny)

    return pmap

def calc_policy_map(closed_set, xw, yw, minx, miny):
    pmap = np.full((xw, yw), np.inf)

    for n in closed_set.values():
        pmap[n.x - minx, n.y - miny] = n.cost

    return pmap

def calc_astar_path(sx, sy, gx, gy, ox, oy, reso, vr):
    """
    sx: start x position [m]
    sy: start y position [m]
    gx: goal x position [m]
    gy: goal y position [m]
    ox: x position list of Obstacles [m]
    oy: y position list of Obstacles [m]
    reso: grid resolution [m]
    """
    nstart = Node(round(sx / reso), round(sy / reso), 0.0, -1)
    ngoal = Node(round(gx / reso), round(gy / reso), 0.0, -1)

    ox = [iox / reso for iox in ox]
    oy = [ioy / reso for ioy in oy]

    obmap, minx, miny, maxx, maxy, xw, yw = calc_obstacle_map(ox, oy, reso, vr)

    open_set, closed_set = {}, {}
    open_set[calc_index(nstart, xw, minx, miny)] = nstart

    motion = get_motion_model()
    pq = []
    heappush(pq, (calc_cost(nstart, ngoal), calc_index(nstart, xw, minx, miny)))

    while pq:
        c_id = heappop(pq)[1]
        current = open_set[c_id]

        if current.x == ngoal.x and current.y == ngoal.y:
            closed_set[c_id] = current
            break

        del open_set[c_id]
        closed_set[c_id] = current

        for i in range(len(motion)):
            node = Node(current.x + motion[i][0], current.y + motion[i][1], current.cost + motion[i][2], c_id)

            if not verify_node(node, minx, miny, xw, yw, obmap):
                continue

            node_ind = calc_index(node, xw, minx, miny)

            if node_ind in closed_set:
                continue

            if node_ind in open_set:
                if open_set[node_ind].cost > node.cost:
                    open_set[node_ind].cost = node.cost
                    open_set[node_ind].pind = c_id
            else:
                open_set[node_ind] = node
                heappush(pq, (calc_cost(node, ngoal), node_ind))

    rx, ry = get_final_path(closed_set, ngoal, nstart, xw, minx, miny, reso)

    return rx, ry

def verify_node(node, minx, miny, xw, yw, obmap):
    if (node.x - minx) >= xw or (node.x - minx) <= 0:
        return False
    if (node.y - miny) >= yw or (node.y - miny) <= 0:
        return False
    if obmap[int(node.x - minx), int(node.y - miny)]:
        return False
    return True

def calc_cost(n, ngoal):
    return n.cost + h(n.x - ngoal.x, n.y - ngoal.y)

def get_motion_model():
    # dx, dy, cost
    return np.array([[1, 0, 1],
                     [0, 1, 1],
                     [-1, 0, 1],
                     [0, -1, 1],
                     [-1, -1, math.sqrt(2)],
                     [-1, 1, math.sqrt(2)],
                     [1, -1, math.sqrt(2)],
                     [1, 1, math.sqrt(2)]])

def calc_index(node, xwidth, xmin, ymin):
    return (node.y - ymin) * xwidth + (node.x - xmin)

def calc_obstacle_map(ox, oy, reso, vr):
    minx = int(min(ox))
    miny = int(min(oy))
    maxx = int(max(ox))
    maxy = int(max(oy))

    xwidth = maxx - minx
    ywidth = maxy - miny

    obmap = np.zeros((xwidth, ywidth), dtype=bool)

    for ix in range(xwidth):
        x = ix + minx
        for iy in range(ywidth):
            y = iy + miny
            if np.min(np.hypot(np.array(ox) - x, np.array(oy) - y)) <= vr / reso:
                obmap[ix, iy] = True

    return obmap, minx, miny, maxx, maxy, xwidth, ywidth

def get_final_path(closed_set, ngoal, nstart, xw, minx, miny, reso):
    rx, ry = [ngoal.x], [ngoal.y]
    nid = calc_index(ngoal, xw, minx, miny)
    
    while True:
        n = closed_set[nid]
        rx.append(n.x)
        ry.append(n.y)
        nid = n.pind
        if rx[-1] == nstart.x and ry[-1] == nstart.y:
            break

    rx = np.array(rx) * reso
    ry = np.array(ry) * reso

    return rx[::-1], ry[::-1]

def h(x, y):
    return math.hypot(x, y)

def main():
    print("A* path planning")

    sx = 10.0  # [m]
    sy = 10.0  # [m]
    gx = 50.0  # [m]
    gy = 50.0  # [m]

    ox, oy = [], []
    for i in range(61):
        ox.append(i)
        oy.append(0.0)
    for i in range(61):
        ox.append(60.0)
        oy.append(i)
    for i in range(61):
        ox.append(i)
        oy.append(60.0)
    for i in range(61):
        ox.append(0.0)
        oy.append(i)
    for i in range(41):
        ox.append(20.0)
        oy.append(i)
    for i in range(41):
        ox.append(40.0)
        oy.append(60.0 - i)

    VEHICLE_RADIUS = 5.0  # [m]
    GRID_RESOLUTION = 1.0  # [m]

    rx, ry = calc_astar_path(sx, sy, gx, gy, ox, oy, GRID_RESOLUTION, VEHICLE_RADIUS)

    plt.plot(ox, oy, ".k")
    plt.plot(sx, sy, "og")
    plt.plot(gx, gy, "xb")
    plt.plot(rx, ry, "-r")
    plt.grid(True)
    plt.axis("equal")
    plt.show()

if __name__ == '__main__':
    main()
