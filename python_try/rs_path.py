import numpy as np
import matplotlib.pyplot as plt

STEP_SIZE = 0.1

class Path:
    def __init__(self):
        self.lengths = []  # lengths of each part of the path (+: forward, -: backward)
        self.ctypes = []  # type of each part of the path
        self.L = 0.0  # total path length
        self.x = []  # final x positions [m]
        self.y = []  # final y positions [m]
        self.yaw = []  # final yaw angles [rad]
        self.directions = []  # forward: 1, backward: -1

def pi_2_pi(iangle):
    while iangle > np.pi:
        iangle -= 2.0 * np.pi
    while iangle < -np.pi:
        iangle += 2.0 * np.pi
    return iangle

def calc_shortest_path(sx, sy, syaw, gx, gy, gyaw, maxc, step_size=STEP_SIZE):
    paths = calc_paths(sx, sy, syaw, gx, gy, gyaw, maxc, step_size=step_size)
    
    minL = float('inf')
    best_path_index = -1
    for i, path in enumerate(paths):
        if path.L <= minL:
            minL = path.L
            best_path_index = i
    
    return paths[best_path_index]

def calc_shortest_path_length(sx, sy, syaw, gx, gy, gyaw, maxc, step_size=STEP_SIZE):
    q0 = [sx, sy, syaw]
    q1 = [gx, gy, gyaw]
    paths = generate_path(q0, q1, maxc)

    minL = float('inf')
    for path in paths:
        L = path.L / maxc
        if L <= minL:
            minL = L

    return minL

def calc_paths(sx, sy, syaw, gx, gy, gyaw, maxc, step_size=STEP_SIZE):
    q0 = [sx, sy, syaw]
    q1 = [gx, gy, gyaw]
    
    paths = generate_path(q0, q1, maxc)
    
    for path in paths:
        x, y, yaw, directions = generate_local_course(path.L, path.lengths, path.ctypes, maxc, step_size * maxc)

        # Convert to global coordinates
        path.x = [np.cos(-q0[2]) * ix + np.sin(-q0[2]) * iy + q0[0] for ix, iy in zip(x, y)]
        path.y = [-np.sin(-q0[2]) * ix + np.cos(-q0[2]) * iy + q0[1] for ix, iy in zip(x, y)]
        path.yaw = [pi_2_pi(iyaw + q0[2]) for iyaw in yaw]
        path.directions = directions
        path.lengths = [l / maxc for l in path.lengths]
        path.L = path.L / maxc

    return paths

def get_label(path):
    label = ""
    for m, l in zip(path.ctypes, path.lengths):
        label += m
        if l > 0.0:
            label += "+"
        else:
            label += "-"
    return label

def polar(x, y):
    r = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    return r, theta

def mod2pi(x):
    v = x % (2.0 * np.pi)
    if v < -np.pi:
        v += 2.0 * np.pi
    elif v > np.pi:
        v -= 2.0 * np.pi
    return v

def LSL(x, y, phi):
    u, t = polar(x - np.sin(phi), y - 1.0 + np.cos(phi))
    if t >= 0.0:
        v = mod2pi(phi - t)
        if v >= 0.0:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0

def LSR(x, y, phi):
    u1, t1 = polar(x + np.sin(phi), y - 1.0 - np.cos(phi))
    u1 = u1 ** 2
    if u1 >= 4.0:
        u = np.sqrt(u1 - 4.0)
        theta = np.arctan2(2.0, u)
        t = mod2pi(t1 + theta)
        v = mod2pi(t - phi)
        if t >= 0.0 and v >= 0.0:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0

def LRL(x, y, phi):
    u1, t1 = polar(x - np.sin(phi), y - 1.0 + np.cos(phi))
    if u1 <= 4.0:
        u = -2.0 * np.arcsin(0.25 * u1)
        t = mod2pi(t1 + 0.5 * u + np.pi)
        v = mod2pi(phi - t + u)
        if t >= 0.0 and u <= 0.0:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0

def set_path(paths, lengths, ctypes):
    path = Path()
    path.ctypes = ctypes
    path.lengths = lengths

    # Check for identical paths
    for tpath in paths:
        if tpath.ctypes == path.ctypes:
            if np.sum(np.abs(np.array(tpath.lengths) - np.array(path.lengths))) <= 0.01:
                return paths  # Do not insert path

    path.L = np.sum(np.abs(lengths))

    if path.L >= 0.01:
        paths.append(path)

    return paths

def SCS(x, y, phi, paths):
    flag, t, u, v = SLS(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["S", "L", "S"])

    flag, t, u, v = SLS(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["S", "R", "S"])

    return paths

def SLS(x, y, phi):
    phi = mod2pi(phi)
    if y > 0.0 and phi > 0.0 and phi < np.pi * 0.99:
        xd = -y / np.tan(phi) + x
        t = xd - np.tan(phi / 2.0)
        u = phi
        v = np.sqrt((x - xd) ** 2 + y ** 2) - np.tan(phi / 2.0)
        return True, t, u, v
    elif y < 0.0 and phi > 0.0 and phi < np.pi * 0.99:
        xd = -y / np.tan(phi) + x
        t = xd - np.tan(phi / 2.0)
        u = phi
        v = -np.sqrt((x - xd) ** 2 + y ** 2) - np.tan(phi / 2.0)
        return True, t, u, v
    return False, 0.0, 0.0, 0.0

def CSC(x, y, phi, paths):
    flag, t, u, v = LSL(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["L", "S", "L"])

    flag, t, u, v = LSL(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["L", "S", "L"])

    flag, t, u, v = LSL(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["R", "S", "R"])

    flag, t, u, v = LSL(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["R", "S", "R"])

    flag, t, u, v = LSR(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["L", "S", "R"])

    flag, t, u, v = LSR(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["L", "S", "R"])

    flag, t, u, v = LSR(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["R", "S", "L"])

    flag, t, u, v = LSR(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["R", "S", "L"])

    return paths

def CCC(x, y, phi, paths):
    flag, t, u, v = LRL(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["L", "R", "L"])

    flag, t, u, v = LRL(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["L", "R", "L"])

    flag, t, u, v = LRL(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, v], ["R", "L", "R"])

    flag, t, u, v = LRL(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, -u, -v], ["R", "L", "R"])

    return paths

def generate_path(q0, q1, maxc):
    dx = q1[0] - q0[0]
    dy = q1[1] - q0[1]
    dth = q1[2] - q0[2]
    c = np.cos(q0[2])
    s = np.sin(q0[2])
    x = c * dx + s * dy
    y = -s * dx + c * dy

    paths = []
    paths = CSC(x / maxc, y / maxc, dth, paths)
    paths = CCC(x / maxc, y / maxc, dth, paths)
    paths = SCS(x / maxc, y / maxc, dth, paths)

    return paths

def generate_local_course(total_length, lengths, modes, maxc, step_size):
    npoint = int(abs(total_length) / step_size) + len(lengths) + 4
    px = [0.0] * npoint
    py = [0.0] * npoint
    pyaw = [0.0] * npoint
    directions = [0] * npoint
    ind = 1
    directions[0] = 1 if lengths[0] > 0.0 else -1

    for (m, length, i) in zip(modes, lengths, range(len(modes))):
        pd = 1 if length > 0.0 else -1
        d = step_size * pd
        if pd == 1:
            n = int(length / step_size)
        else:
            n = int(abs(length) / step_size)
        if length == 0:
            n = 0

        x, y, yaw = 0.0, 0.0, 0.0
        for i in range(n):
            x, y, yaw = step_xyd(x, y, yaw, d, m, maxc)
            px[ind] = x
            py[ind] = y
            pyaw[ind] = yaw
            directions[ind] = pd
            ind += 1

        pd = -pd
        x, y, yaw = step_xyd(x, y, yaw, length - d * n, m, maxc)
        px[ind] = x
        py[ind] = y
        pyaw[ind] = yaw
        directions[ind] = pd
        ind += 1

    while px[ind - 1] == px[ind]:
        ind -= 1

    return px[:ind], py[:ind], pyaw[:ind], directions[:ind]

def step_xyd(x, y, yaw, d, mode, maxc):
    if mode == "S":
        x += d * np.cos(yaw)
        y += d * np.sin(yaw)
    elif mode == "L":
        l = d / maxc
        x += np.sin(l) * np.cos(yaw)
        y += np.sin(l) * np.sin(yaw)
        yaw += l
    elif mode == "R":
        l = d / maxc
        x += -np.sin(l) * np.cos(yaw)
        y += -np.sin(l) * np.sin(yaw)
        yaw -= l

    return x, y, yaw


import numpy as np

def mod2pi(x):
    return x % (2 * np.pi)

def calc_tauOmega(u, v, xi, eta, phi):
    delta = mod2pi(u - v)
    A = np.sin(u) - np.sin(delta)
    B = np.cos(u) - np.cos(delta) - 1.0

    t1 = np.arctan2(eta * A - xi * B, xi * A + eta * B)
    t2 = 2.0 * (np.cos(delta) - np.cos(v) - np.cos(u)) + 3.0

    if t2 < 0:
        tau = mod2pi(t1 + np.pi)
    else:
        tau = mod2pi(t1)

    omega = mod2pi(tau - u + v - phi)
    return tau, omega

def LRLRn(x, y, phi):
    xi = x + np.sin(phi)
    eta = y - 1.0 - np.cos(phi)
    rho = 0.25 * (2.0 + np.sqrt(xi * xi + eta * eta))

    if rho <= 1.0:
        u = np.arccos(rho)
        t, v = calc_tauOmega(u, -u, xi, eta, phi)
        if t >= 0.0 and v <= 0.0:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0

def LRLRp(x, y, phi):
    xi = x + np.sin(phi)
    eta = y - 1.0 - np.cos(phi)
    rho = (20.0 - xi * xi - eta * eta) / 16.0

    if 0.0 <= rho <= 1.0:
        u = -np.arccos(rho)
        if u >= -0.5 * np.pi:
            t, v = calc_tauOmega(u, u, xi, eta, phi)
            if t >= 0.0 and v >= 0.0:
                return True, t, u, v
    return False, 0.0, 0.0, 0.0

def set_path(paths, lengths, modes):
    # This is a placeholder for the actual implementation of path setting
    paths.append({'lengths': lengths, 'modes': modes})
    return paths

def CCCC(x, y, phi, paths):
    flag, t, u, v = LRLRn(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, -u, v], ["L", "R", "L", "R"])

    flag, t, u, v = LRLRn(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, u, -v], ["L", "R", "L", "R"])

    flag, t, u, v = LRLRn(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, -u, v], ["R", "L", "R", "L"])

    flag, t, u, v = LRLRn(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, -u, u, -v], ["R", "L", "R", "L"])

    flag, t, u, v = LRLRp(x, y, phi)
    if flag:
        paths = set_path(paths, [t, u, u, v], ["L", "R", "L", "R"])

    flag, t, u, v = LRLRp(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, -u, -u, -v], ["L", "R", "L", "R"])

    flag, t, u, v = LRLRp(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, u, u, v], ["R", "L", "R", "L"])

    flag, t, u, v = LRLRp(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, -u, -u, -v], ["R", "L", "R", "L"])

    return paths

def polar(x, y):
    rho = np.hypot(x, y)
    theta = np.arctan2(y, x)
    return rho, theta

def LRSR(x, y, phi):
    xi = x + np.sin(phi)
    eta = y - 1.0 - np.cos(phi)
    rho, theta = polar(-eta, xi)

    if rho >= 2.0:
        t = theta
        u = 2.0 - rho
        v = mod2pi(t + 0.5 * np.pi - phi)
        if t >= 0.0 and u <= 0.0 and v <= 0.0:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0

def LRSL(x, y, phi):
    xi = x - np.sin(phi)
    eta = y - 1.0 + np.cos(phi)
    rho, theta = polar(xi, eta)

    if rho >= 2.0:
        r = np.sqrt(rho * rho - 4.0)
        u = 2.0 - r
        t = mod2pi(theta + np.arctan2(r, -2.0))
        v = mod2pi(phi - 0.5 * np.pi - t)
        if t >= 0.0 and u <= 0.0 and v <= 0.0:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0

def CCSC(x, y, phi, paths):
    flag, t, u, v = LRSL(x, y, phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * np.pi, u, v], ["L", "R", "S", "L"])

    flag, t, u, v = LRSL(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * np.pi, -u, -v], ["L", "R", "S", "L"])

    flag, t, u, v = LRSL(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * np.pi, u, v], ["R", "L", "S", "R"])

    flag, t, u, v = LRSL(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * np.pi, -u, -v], ["R", "L", "S", "R"])

    flag, t, u, v = LRSR(x, y, phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * np.pi, u, v], ["L", "R", "S", "R"])

    flag, t, u, v = LRSR(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * np.pi, -u, -v], ["L", "R", "S", "R"])

    flag, t, u, v = LRSR(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * np.pi, u, v], ["R", "L", "S", "L"])

    flag, t, u, v = LRSR(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * np.pi, -u, -v], ["R", "L", "S", "L"])

    return paths

def CCSCC(x: float, y: float, phi: float, paths: list):
    flag, t, u, v = LRSLR(x, y, phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * np.pi, u, -0.5 * np.pi, v], ["L", "R", "S", "L", "R"])

    flag, t, u, v = LRSLR(-x, y, -phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * np.pi, -u, 0.5 * np.pi, -v], ["L", "R", "S", "L", "R"])

    flag, t, u, v = LRSLR(x, -y, -phi)
    if flag:
        paths = set_path(paths, [t, -0.5 * np.pi, u, -0.5 * np.pi, v], ["R", "L", "S", "R", "L"])

    flag, t, u, v = LRSLR(-x, -y, phi)
    if flag:
        paths = set_path(paths, [-t, 0.5 * np.pi, -u, 0.5 * np.pi, -v], ["R", "L", "S", "R", "L"])

    return paths


def generate_local_course(L, lengths, mode, maxc, step_size):
    npoint = int(L / step_size) + len(lengths) + 3

    px = np.zeros(npoint)
    py = np.zeros(npoint)
    pyaw = np.zeros(npoint)
    directions = np.zeros(npoint, dtype=int)
    ind = 1

    if lengths[0] > 0.0:
        directions[0] = 1
    else:
        directions[0] = -1

    if lengths[0] > 0.0:
        d = step_size
    else:
        d = -step_size

    pd = d
    ll = 0.0

    for m, l, i in zip(mode, lengths, range(len(mode))):

        if l > 0.0:
            d = step_size
        else:
            d = -step_size

        ox, oy, oyaw = px[ind], py[ind], pyaw[ind]

        ind -= 1
        if i >= 1 and (lengths[i-1] * lengths[i]) > 0:
            pd = -d - ll
        else:
            pd = d - ll

        while abs(pd) <= abs(l):
            ind += 1
            px, py, pyaw, directions = interpolate(ind, pd, m, maxc, ox, oy, oyaw, px, py, pyaw, directions)
            pd += d

        ll = l - pd - d

        ind += 1
        px, py, pyaw, directions = interpolate(ind, l, m, maxc, ox, oy, oyaw, px, py, pyaw, directions)

    # Remove unused data
    while px[-1] == 0.0:
        px = np.delete(px, -1)
        py = np.delete(py, -1)
        pyaw = np.delete(pyaw, -1)
        directions = np.delete(directions, -1)

    return px, py, pyaw, directions


def interpolate(ind, l, m, maxc, ox, oy, oyaw, px, py, pyaw, directions):
    if m == "S":
        px[ind] = ox + l / maxc * math.cos(oyaw)
        py[ind] = oy + l / maxc * math.sin(oyaw)
        pyaw[ind] = oyaw
    else:  # curve
        ldx = math.sin(l) / maxc
        if m == "L":  # left turn
            ldy = (1.0 - math.cos(l)) / maxc
        elif m == "R":  # right turn
            ldy = (1.0 - math.cos(l)) / -maxc

        gdx = math.cos(-oyaw) * ldx + math.sin(-oyaw) * ldy
        gdy = -math.sin(-oyaw) * ldx + math.cos(-oyaw) * ldy
        px[ind] = ox + gdx
        py[ind] = oy + gdy

    if m == "L":  # left turn
        pyaw[ind] = oyaw + l
    elif m == "R":  # right turn
        pyaw[ind] = oyaw - l

    if l > 0.0:
        directions[ind] = 1
    else:
        directions[ind] = -1

    return px, py, pyaw, directions


def generate_path(q0, q1, maxc):
    dx = q1[0] - q0[0]
    dy = q1[1] - q0[1]
    dth = q1[2] - q0[2]
    c = math.cos(q0[2])
    s = math.sin(q0[2])
    x = (c * dx + s * dy) * maxc
    y = (-s * dx + c * dy) * maxc

    paths = []
    paths = SCS(x, y, dth, paths)
    paths = CSC(x, y, dth, paths)
    paths = CCC(x, y, dth, paths)
    paths = CCCC(x, y, dth, paths)
    paths = CCSC(x, y, dth, paths)
    paths = CCSCC(x, y, dth, paths)

    return paths


def calc_curvature(x, y, yaw, directions):
    c = []
    ds = []

    for i in range(1, len(x) - 1):
        dxn = x[i] - x[i - 1]
        dxp = x[i + 1] - x[i]
        dyn = y[i] - y[i - 1]
        dyp = y[i + 1] - y[i]
        dn = math.sqrt(dxn ** 2 + dyn ** 2)
        dp = math.sqrt(dxp ** 2 + dyp ** 2)
        dx = 1.0 / (dn + dp) * (dp / dn * dxn + dn / dp * dxp)
        ddx = 2.0 / (dn + dp) * (dxp / dp - dxn / dn)
        dy = 1.0 / (dn + dp) * (dp / dn * dyn + dn / dp * dyp)
        ddy = 2.0 / (dn + dp) * (dyp / dp - dyn / dn)
        curvature = (ddy * dx - ddx * dy) / (dx ** 2 + dy ** 2)
        d = (dn + dp) / 2.0

        if math.isnan(curvature):
            curvature = 0.0

        if directions[i] <= 0.0:
            curvature = -curvature

        if len(c) == 0:
            ds.append(d)
            c.append(curvature)

        ds.append(d)
        c.append(curvature)

    ds.append(ds[-1])
    c.append(c[-1])

    return c, ds


def check_path(start_x, start_y, start_yaw, end_x, end_y, end_yaw, max_curvature):
    paths = calc_paths(start_x, start_y, start_yaw, end_x, end_y, end_yaw, max_curvature)

    assert len(paths) >= 1

    for path in paths:
        assert abs(path.x[0] - start_x) <= 0.01
        assert abs(path.y[0] - start_y) <= 0.01
        assert abs(path.yaw[0] - start_yaw) <= 0.01
        assert abs(path.x[-1] - end_x) <= 0.01
        assert abs(path.y[-1] - end_y) <= 0.01
        assert abs(path.yaw[-1] - end_yaw) <= 0.01

        d = [math.sqrt(dx**2 + dy**2) for dx, dy in zip(np.diff(path.x[:-1]), np.diff(path.y[:-1]))]

        for i in range(len(d)):
            assert abs(d[i] - STEP_SIZE) <= 0.001


def test():
    print("Test1")
    start_x, start_y, start_yaw = 0.0, 0.0, math.radians(10.0)
    end_x, end_y, end_yaw = 7.0, -8.0, math.radians(50.0)
    max_curvature = 2.0

    check_path(start_x, start_y, start_yaw, end_x, end_y, end_yaw, max_curvature)

    start_x, start_y, start_yaw = 0.0, 0.0, math.radians(10.0)
    end_x, end_y, end_yaw = 7.0, -8.0, math.radians(-50.0)
    max_curvature = 2.0

    check_path(start_x, start_y, start_yaw, end_x, end_y, end_yaw, max_curvature)

    start_x, start_y, start_yaw = 0.0, 10.0, math.radians(-10.0)
    end_x, end_y, end_yaw = -7.0, -8.0, math.radians(-50.0)
    max_curvature = 2.0

    check_path(start_x, start_y, start_yaw, end_x, end_y, end_yaw, max_curvature)

    start_x, start_y, start_yaw = 0.0, 10.0, math.radians(-10.0)
    end_x, end_y, end_yaw = -7.0, -8.0, math.radians(150.0)
    max_curvature = 1.0

    check_path(start_x, start_y, start_yaw, end_x, end_y, end_yaw, max_curvature)

    start_x, start_y, start_yaw = 0.0, 10.0, math.radians(-10.0)
    end_x, end_y, end_yaw = 7.0, 8.0, math.radians(50.0)
    max_curvature = 1.0

    check_path(start_x, start_y, start_yaw, end_x, end_y, end_yaw, max_curvature)

    print("Done")


if __name__ == "__main__":
    test()

# # Test the code with example start and goal points
# def plot_path(path):
#     plt.plot(path.x, path.y, label=get_label(path))
#     plt.grid(True)
#     plt.axis("equal")
#     plt.show()

# # Example usage
# sx, sy, syaw = 0, 0, np.deg2rad(45)  # start x, y, yaw
# gx, gy, gyaw = 5, 5, np.deg2rad(-45)  # goal x, y, yaw
# max_curvature = 0.1

# path = calc_shortest_path(sx, sy, syaw, gx, gy, gyaw, max_curvature)
# plot_path(path)
