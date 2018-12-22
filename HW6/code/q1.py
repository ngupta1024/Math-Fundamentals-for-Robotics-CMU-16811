import numpy as np
import matplotlib.pyplot as plt
import sys


def decide2eliminate(p1, p2, p3):
    v1 = p1 - p2
    v2 = p3 - p2

    cross_prod = np.cross(v1, v2)

    if cross_prod >= 0:
        return True
    else:
        return False


def generate_points(n, max_coord):
    p = np.random.randint(0, max_coord, (n, 2))

    center = np.sum(p, axis=0)/p.shape[0]

    return p, center


def sort_points(center, p):
    sorted_p = p - center

    sorted_p = sorted(sorted_p, key=lambda x: (np.arctan2(x[1], x[0]), x[0]**2+x[1]**2))

    return sorted_p


def process(p):
    l, ind, r = 0, 1, 2
    eliminate_set = set()
    keep_set = set()
    num = len(p)
    start = ind
    while len(eliminate_set | keep_set) != num or ind != start:
        if ind in eliminate_set:
            ind = (ind+1) % num
            continue

        eliminate = decide2eliminate(p[l], p[ind], p[r])

        if eliminate:
            if ind in keep_set:
                keep_set.remove(ind)
            eliminate_set.add(ind)
            if ind == start:
                start = min(keep_set) if keep_set else l
            ind = l
            l = (l-1)%num
            while l in eliminate_set:
                l = (l-1)%num
        else:
            if ind in eliminate_set:
                eliminate_set.remove(ind)
            keep_set.add(ind)
            l = ind
            ind = (ind+1) % num
            r = (r+1)%num
            while r in eliminate_set:
                r = (r+1)%num

    while len(eliminate_set | keep_set) == num and ind == start:
        if not decide2eliminate(p[l], p[ind], p[r]):
            break
        while decide2eliminate(p[l], p[ind], p[r]):
            if ind in keep_set:
                keep_set.remove(ind)
            eliminate_set.add(ind)
            ind = l
            l = (l - 1) % num
            while l in eliminate_set:
                l = (l - 1) % num
        start = min(keep_set)

        l = ind
        ind = start
        r = (r+1)%num
        while r in eliminate_set:
            r = (r + 1) % num
        if not decide2eliminate(p[l], p[ind], p[r]):
            break

    keep_points = [p[i] for i in range(num) if i in keep_set]

    return keep_points


def plot_polygon(p, origin, keep_points):
    x = [pt[0] for pt in p]
    y = [pt[1] for pt in p]
    plt.figure()
    plt.scatter(x, y, marker='.')
    plt.scatter(origin[0], origin[1], marker='.', color='r')

    kx = [pt[0] for pt in keep_points]
    ky = [pt[1] for pt in keep_points]
    for i in range(len(keep_points)-1):
        xt, yt = kx[i:i+2], ky[i:i+2]
        plt.plot(xt, yt, color='r')

    xt, yt = [kx[-1], kx[0]], [ky[-1], ky[0]]
    plt.plot(xt, yt, color='r')

    plt.show()


if __name__ == '__main__':
    method = ''
    if len(sys.argv) == 1:
        method = 'random'
    elif len(sys.argv) == 2:
        method = 'specified'
    else:
        print('Error in command. ')

    n = 100
    max_coord = 10000

    p, center = None, None
    if method == 'random':
        p, center = generate_points(n, max_coord)
    else:
        pt_file = sys.argv[1]
        p = []
        with open(pt_file, 'r') as f:
            for line in f.readlines():
                tokens = line.strip('\n').replace(' ', '').split(',')
                p.append([float(tokens[0]), float(tokens[1])])

        if len(p) == 0:
            print('With method "specified", a list of points should be provided.')
            exit(1)
        p = np.array(p)
        center = np.sum(p, axis=0)/p.shape[0]

    sorted_p = sort_points(center, p)

    keep_p = process(sorted_p)
    keep_p = np.array(keep_p)+center

    plot_polygon(p, center, keep_p)
