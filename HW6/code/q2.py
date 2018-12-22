import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


def plot_poly(polys, start, end, graph, point_dict, pass_points):
    plt.figure()
    plt.scatter(start[0], start[1], marker='o')
    plt.scatter(end[0], end[1], marker='o')
    for key, neighs in graph.items():
        for n, _ in neighs:
            xt, yt = [point_dict[key][0], point_dict[n][0]], [point_dict[key][1], point_dict[n][1]]
            plt.plot(xt, yt, color='g', linestyle=':')

    for poly in polys:
        x = [pt[0] for pt in poly]
        y = [pt[1] for pt in poly]
        plt.scatter(x, y, marker='.')
        for i in range(poly.shape[0]-1):
            xt, yt = poly[i:i+2, 0], poly[i:i+2, 1]
            plt.plot(xt, yt, color='r')

        xt, yt = [poly[-1, 0], poly[0, 0]], [poly[-1, 1], poly[0, 1]]
        plt.plot(xt, yt, color='r')

    for i in range(len(pass_points)-1):
        xt = [point_dict[pass_points[i]][0], point_dict[pass_points[i+1]][0]]
        yt = [point_dict[pass_points[i]][1], point_dict[pass_points[i+1]][1]]
        plt.plot(xt, yt, color='b', linewidth=2)

    plt.show()


def check_intersect(p1, p2, p3, p4):
    p1x, p1y = p1
    p2x, p2y = p2
    p3x, p3y = p3
    p4x, p4y = p4

    l1 = np.array([p2y-p1y, p1x-p2x, -p1x*(p2y-p1y)+p1y*(p2x-p1x)])
    l2 = np.array([p4y-p3y, p3x-p4x, -p3x*(p4y-p3y)+p3y*(p4x-p3x)])

    l1 /= np.sqrt(np.sum(l1[:2]**2))
    l2 /= np.sqrt(np.sum(l2[:2]**2))

    div = l1[0]*l2[1]-l1[1]*l2[0]

    if div == 0:
        # parallel
        l1 = -l1 if l1[0]*l2[0] < 0 else l1
        if np.sum(np.abs(l2-l1)) < 1e-5:
            return False
        else:
            return False
    else:
        x = -(l1[2]*l2[1]-l1[1]*l2[2])/div
        y = -(l1[0]*l2[2]-l1[2]*l2[0])/div

        pp = np.array([x, y])
        if (np.sum(np.abs(pp-p1)) < 1e-5 or np.sum(np.abs(pp-p2)) < 1e-5 or
            np.sum(np.abs(pp-p3))<1e-5 or np.sum(np.abs(pp-p4))<1e-5):
            return False
        elif ((min(p1x, p2x) < x < max(p1x, p2x) or abs(x-p1x) < 1e-5 or abs(x-p2x)<1e-5) and
              (min(p3x, p4x) < x < max(p3x, p4x) or abs(x-p3x) < 1e-5 or abs(x-p4x)<1e-5) and
              (min(p1y, p2y) < y < max(p1y, p2y) or abs(y-p1y) < 1e-5 or abs(y-p2y)<1e-5) and
              (min(p3y, p4y) < y < max(p3y, p4y) or abs(y-p3y) < 1e-5 or abs(y-p4y)<1e-5)):
            return True
        else:
            return False


def points_in_polygon(p, polygons):
    x, y = p
    status = False
    for poly in polygons:
        i, j = 0, len(poly)-1
        status = False

        for i in range(len(poly)):
            if (poly[i][0] == x and poly[i][1] == y) or (poly[j][0] == x and poly[j][1] == y):
                status = False
                break
            if (((poly[i][1] < y and poly[j][1] >= y) or (poly[j][1] < y and poly[i][1] >= y)) and
                    (poly[i][0] <= x or poly[j][0] <= x)):
                status ^= (poly[i][0]+(y-poly[i][1])/(poly[j][1]-poly[i][1]) * (poly[j][0]-poly[i][0]) < x)
            j = i

        if status:
            break

    return status


def construct_graph(start, end, polygons, show_path=True):
    points = []
    lines = []
    graph = defaultdict(list)
    same_polygon = []
    count = 1
    for pk, poly in enumerate(polygons):
        pset = set()
        cur_points = [poly[i] for i in range(poly.shape[0]) if not points_in_polygon(poly[i], polygons)]
        points.extend(cur_points)
        lines.extend([[poly[i], poly[i+1]] for i in range(poly.shape[0]-1)])
        lines.append([poly[-1], poly[0]])
        for i in range(poly.shape[0]):
            if points_in_polygon(poly[i], polygons):
                continue
            pset.add(count)
            if i == poly.shape[0]-1:
                if not points_in_polygon(poly[0], polygons):
                    dist = np.sqrt(np.sum((poly[i] - poly[0]) ** 2))
                    graph[count-len(cur_points)+1].append((count, dist))
                    graph[count].append((count-len(cur_points)+1, dist))
            else:
                if not points_in_polygon(poly[i+1], polygons):
                    dist = np.sqrt(np.sum((poly[i] - poly[i + 1]) ** 2))
                    c = (poly[i] + poly[i+1])/2
                    if not points_in_polygon(c, polygons[:pk]+polygons[pk+1:]):
                        graph[count].append((count+1, dist))
                        graph[count+1].append((count, dist))
            count += 1

        same_polygon.append(pset)
    point_dict = {i+1: points[i] for i in range(len(points))}
    point_dict[0] = start
    point_dict[len(points)+1] = end

    for i in range(1, len(points)+1):
        for j in range(i+1, len(points)+1):
            flag = 0
            for item in same_polygon:
                if i in item and j in item:
                    flag = 1
                    break
            if flag:
                continue
            c = (point_dict[i] + point_dict[j]) / 2
            if points_in_polygon(c, polygons):
                continue
            intersect = False
            for l in lines:
                if check_intersect(point_dict[i], point_dict[j], l[0], l[1]):
                    intersect = True
                    break
            if not intersect:
                dist = np.sqrt(np.sum((point_dict[i]-point_dict[j])**2))
                graph[i].append((j, dist))
                graph[j].append((i, dist))

    for i, p in enumerate([start, end]):
        for j in range(1, len(points)+1):
            intersect = False
            for l in lines:
                if check_intersect(p, point_dict[j], l[0], l[1]):
                    intersect = True
                    break
            if not intersect:
                c = (point_dict[j] + p) / 2
                if points_in_polygon(c, polygons):
                    continue
                dist = np.sqrt(np.sum((p - point_dict[j]) ** 2))
                ind = 0 if i == 0 else len(points)+1
                graph[ind].append((j, dist))
                graph[j].append((ind, dist))

    pass_inds, path_length = get_min_dist(graph)
    pass_points = [point_dict[i] for i in pass_inds][::-1]
    print(pass_inds)
    print(pass_points)

    if show_path:
        plot_poly(polygons, start, end, graph, point_dict, pass_inds)

    return pass_points, path_length, graph, point_dict


def get_min_dist(graph):
    dist = [float('inf') for _ in range(len(graph))]
    cur = [0]
    dist[0] = 0
    while cur:
        next = []
        for item in cur:
            for neigh, d in graph[item]:
                nd = dist[item] + d
                if nd < dist[neigh]:
                    dist[neigh] = nd
                    next.append(neigh)
        cur = next

    path_length = dist[-1]
    pass_points = [len(dist)-1]
    target = len(dist)-1
    while pass_points[-1] != 0:
        for neigh, d in graph[target]:
            if abs(dist[target]-d - dist[neigh]) < 1e-3:
                pass_points.append(neigh)
                target = neigh
                path_length -= d
                break

    return pass_points, dist[-1]


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Error in command. '
              'python q2.py [POLYGON_TXT]')
        exit(1)

    polygon_txt = sys.argv[1]

    polygons = []
    start, end = None, None
    with open(polygon_txt, 'r') as f:
        for line in f.readlines():
            if line.startswith('start'):
                x, y = line.strip('\n').replace(' ', '').split(':')[-1][1:-1].split(',')
                start = np.array([float(x), float(y)])
                continue
            elif line.startswith('end'):
                x, y = line.strip('\n').replace(' ', '').split(':')[-1][1:-1].split(',')
                end = np.array([float(x), float(y)])
                continue
            polygon = []
            tokens = line.strip('\n').replace(' ', '')[1:-1].split('),(')
            for p in tokens:
                x, y = p.split(',')
                polygon.append([float(x), float(y)])
            xs = [a[0] for a in polygon]
            ys = [a[1] for a in polygon]
            center = np.array([np.sum(xs)/len(xs), np.sum(ys)/len(ys)])
            polygon = np.array(polygon) - center
            polygon = sorted(polygon, key=lambda a: (np.arctan2(a[1], a[0]), a[0]**2+a[1]**2))
            polygon += center

            polygons.append(polygon)

    if points_in_polygon(start, polygons) or points_in_polygon(end, polygons):
        print('There is no path from the start point to the goal point.')

    else:
        pass_points, path_length, _, _ = construct_graph(start, end, polygons)
        print("Minimum path: {}".format(pass_points))
        print("Minimum distance: {}".format(path_length))
