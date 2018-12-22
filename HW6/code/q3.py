import sys
import numpy as np
import matplotlib.pyplot as plt
import q1
import q2


def parse_input(input_txt):
    start, end = None, None
    start_points, end_points = None, None
    polygons = []
    polygons_points = []
    with open(input_txt, 'r') as f:
        for line in f.readlines():
            if line.startswith('start') or line.startswith('end'):
                points = []
                tokens = line.strip('\n').replace(' ', '').split(':')[-1][1:-1].split('),(')
                for pt in tokens:
                    x, y = pt.split(',')
                    points.append([float(x), float(y)])
                points = np.array(points)
                center = np.sum(points, axis=0) / points.shape[0]
                sorted_p = q1.sort_points(center, points)
                if line.startswith('start'):
                    start_points = points
                    start = q1.process(sorted_p)+center
                else:
                    end_points = points
                    end = q1.process(sorted_p)+center
                continue
            else:
                points = []
                tokens = line.strip('\n').replace(' ', '')[1:-1].split('),(')
                for pt in tokens:
                    x, y = pt.split(',')
                    points.append([float(x), float(y)])
                points = np.array(points)
                center = np.sum(points, axis=0)/points.shape[0]
                sorted_p = q1.sort_points(center, points)
                sorted_p = q1.process(sorted_p)
                polygons.append(sorted_p+center)
                polygons_points.append(points)

                # q1.plot_polygon(points, center, sorted_p+center)

    return start, start_points, end, end_points, polygons, polygons_points


def form_obstacle(robot, polygons):
    joined_polys = []

    robot_edges = []
    sort_pt = sorted([(i, robot[i]) for i in range(robot.shape[0])], key=lambda x: np.sum(x[1]**2))
    origin_ind = sort_pt[0][0]
    robot = -robot
    for i in range(robot.shape[0]-1):
        robot_edges.append((robot[i+1]-robot[i], None))
    robot_edges.append((robot[0]-robot[-1], None))
    robot_edges = robot_edges[origin_ind:]+robot_edges[:origin_ind]
    robot = -robot

    for poly in polygons:
        poly_edges = []
        for i in range(poly.shape[0]-1):
            poly_edges.append((poly[i+1]-poly[i], poly[i]))
        poly_edges.append((poly[0]-poly[-1], poly[-1]))

        # joined with robot edges
        poly_edges.extend(robot_edges)
        poly_edges = sorted(poly_edges, key=lambda x: np.arctan2(x[0][1], x[0][0]))

        # reconstruct edges
        begin_ind = 0
        while not (poly_edges[begin_ind][1] is None and poly_edges[begin_ind+1][1] is not None):
            begin_ind += 1
        pt = np.copy(poly_edges[begin_ind+1][1])
        start_edge = poly_edges[begin_ind][0]
        for item in robot_edges:
            if item[0][0] == start_edge[0] and item[0][1] == start_edge[1]:
                break
            pt += item[0]
        poly_edges = poly_edges[begin_ind:]+poly_edges[:begin_ind]

        new_poly_points = []
        for item in poly_edges:
            pt += item[0]
            new_poly_points.append(np.copy(pt))

        joined_polys.append(np.array(new_poly_points))

    return joined_polys, origin_ind


def plot_poly(polys, joined_polys, poly_points, start, start_points, end, end_points, graph, point_dict, pass_points):
    plt.figure()
    for i in range(len(poly_points)):
        plt.scatter(poly_points[i][:, 0], poly_points[i][:, 1], marker='.')
    plt.scatter(start_points[:, 0], start_points[:, 1])
    plt.scatter(end_points[:, 0], end_points[:, 1])

    for key, neighs in graph.items():
        for n, _ in neighs:
            xt, yt = [point_dict[key][0], point_dict[n][0]], [point_dict[key][1], point_dict[n][1]]
            plt.plot(xt, yt, color='g', linestyle=':')

    # plot start
    x = [pt[0] for pt in start]
    y = [pt[1] for pt in start]
    plt.scatter(x, y, marker='.')
    for i in range(start.shape[0]-1):
        xt, yt = start[i:i+2, 0], start[i:i+2, 1]
        plt.plot(xt, yt, color='y')
    xt, yt = [start[-1, 0], start[0, 0]], [start[-1, 1], start[0, 1]]
    plt.plot(xt, yt, color='y')

    # plot goal
    x = [pt[0] for pt in end]
    y = [pt[1] for pt in end]
    plt.scatter(x, y, marker='.')
    for i in range(end.shape[0] - 1):
        xt, yt = end[i:i + 2, 0], end[i:i + 2, 1]
        plt.plot(xt, yt, color='y')
    xt, yt = [end[-1, 0], end[0, 0]], [end[-1, 1], end[0, 1]]
    plt.plot(xt, yt, color='y')

    for poly in polys:
        for i in range(poly.shape[0]-1):
            xt, yt = poly[i:i+2, 0], poly[i:i+2, 1]
            plt.plot(xt, yt, color='r', linestyle=':')
        xt, yt = [poly[-1, 0], poly[0, 0]], [poly[-1, 1], poly[0, 1]]
        plt.plot(xt, yt, color='r', linestyle=':')
    for poly in joined_polys:
        for i in range(poly.shape[0]-1):
            xt, yt = poly[i:i+2, 0], poly[i:i+2, 1]
            plt.plot(xt, yt, color='r')
        xt, yt = [poly[-1, 0], poly[0, 0]], [poly[-1, 1], poly[0, 1]]
        plt.plot(xt, yt, color='r')

    for i in range(len(pass_points)-1):
        xt = [pass_points[i][0], pass_points[i+1][0]]
        yt = [pass_points[i][1], pass_points[i+1][1]]
        plt.plot(xt, yt, color='b', linewidth=2)

    plt.show()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Error in command. '
              'python q3.py [ROBOT_TXT]')
        exit(1)

    robot_txt = sys.argv[1]

    start, start_points, end, end_points, polygons, polygons_points = parse_input(robot_txt)

    joined_polys, origin_ind = form_obstacle(start, polygons)

    if q2.points_in_polygon(start[origin_ind], joined_polys) or q2.points_in_polygon(end[origin_ind], joined_polys):
        print('There is no path from the start point to the goal point.')
    else:
        pass_points, path_length, graph, point_dict = q2.construct_graph(start[origin_ind], end[origin_ind], joined_polys, show_path=False)
        plot_poly(polygons, joined_polys, polygons_points, start, start_points, end, end_points, graph, point_dict, pass_points)
        print("Minimum path: {}".format(pass_points))
        print("Minimum distance: {}".format(path_length))
