# -*- encoding=utf-8 -*-
#
# fast method for searching network-constrained point knn
import numpy as np

from numpy import array, tile


def init_interval_neighbors_on_same_road(points_on_road, step_num):
    interval_neighbors, interval_dists, interval_counts = dict(), dict(), dict()
    for pid in points_on_road:
        interval_neighbors[pid] = [list() for i in range(step_num)]
        interval_dists[pid] = [list() for i in range(step_num)]
        interval_counts[pid] = 0
    return interval_neighbors, interval_dists, interval_counts


def create_interval_neighbors_on_same_road(points_on_road, match_dists, step_size, step_num):
    interval_neighbors, interval_dists, interval_counts = init_interval_neighbors_on_same_road(points_on_road, step_num)
    num = len(points_on_road)
    if num == 1:
        return interval_neighbors, interval_dists, interval_counts
    for i in range(num - 1):
        p1 = points_on_road[i]
        for j in range(i + 1, num):
            p2 = points_on_road[j]
            dist = round(abs(match_dists[p1][0] - match_dists[p2][0]), 3)
            index = int(dist / step_size)
            if index < step_num:
                interval_neighbors[p1][index].append(p2)
                interval_dists[p1][index].append(dist)

                interval_neighbors[p2][index].append(p1)
                interval_dists[p2][index].append(dist)
    return interval_neighbors, interval_dists, interval_counts


def create_candidate_neighbor_dist_matrix(road_neighbors, road_contains, match_dists):
    node_node_dist_matrix = list()
    node_point_dist_matrix = list()
    pids = list()
    for rid in road_neighbors.keys():
        pids.extend(road_contains[rid])
        for pid in road_contains[rid]:
            node_node_dist_matrix.append(road_neighbors[rid])
            node_point_dist_matrix.append([match_dists[pid][0], match_dists[pid][1],
                                           match_dists[pid][0], match_dists[pid][1]])
    return array(pids), array(node_node_dist_matrix) + array(node_point_dist_matrix)


def get_given_interval_road_neighbors(road_neighbors, start_dist, end_dist):
    interval_road_neighbors = dict()
    rids = list()
    for rid in road_neighbors.keys():
        min_dist = min(road_neighbors[rid])
        if start_dist <= min_dist < end_dist:
            rids.append(rid)
            interval_road_neighbors[rid] = road_neighbors[rid]
    return rids, interval_road_neighbors


def calc_network_dist_by_matrix(node_point_dist, pids, other_dist_matrix, max_size):
    node_point_matrix = tile(node_point_dist, (len(pids), 1))
    dists_matrix = node_point_matrix + other_dist_matrix
    network_dist = dists_matrix.min(axis=1)
    flag = network_dist < max_size
    return pids[flag], network_dist[flag]


def assign_neighbors_2_intervals(interval_neighbors, interval_dists, neighbors, dists, step_size):
    indexes = (dists / step_size).astype(int)
    for i in range(len(indexes)):
        index, n, d = indexes[i], neighbors[i], dists[i]
        interval_neighbors[index].append(n)
        interval_dists[index].append(d)


def create_interval_neighbors(rpids, interval_index, step_size, road_neighbor, road_contains, match_dists,
                              interval_neighbors, interval_dists, interval_counts, max_size):
    start_dist, end_dist = interval_index * step_size, (interval_index + 1) * step_size
    rids, interval_road_neighbors = get_given_interval_road_neighbors(road_neighbor, start_dist, end_dist)
    pids, dist_matrix = create_candidate_neighbor_dist_matrix(interval_road_neighbors, road_contains, match_dists)
    if len(pids) > 0:
        for rpid in rpids:
            node_point_dist = array([match_dists[rpid][0], match_dists[rpid][0],
                                     match_dists[rpid][1], match_dists[rpid][1]])
            neighbors, network_dist = calc_network_dist_by_matrix(node_point_dist, pids, dist_matrix, max_size)
            assign_neighbors_2_intervals(interval_neighbors[rpid], interval_dists[rpid], neighbors, network_dist,
                                         step_size)
            interval_counts[rpid] += len(interval_neighbors[rpid][interval_index])
    for rid in rids:
        del road_neighbor[rid]


def assign_enough(current_interval, interval_counts, unenough, k):
    enough = [(i, current_interval) for i in unenough if interval_counts[i] >= k]
    for i in enough:
        unenough.remove(i[0])
    return enough


def get_knn_from_interval_neighbors(interval_neighbors, interval_dists, enough_intervals, k):
    knn_dist = dict()
    for pid in interval_neighbors.keys():
        neighbors, dists = list(), list()
        if enough_intervals[pid] < 0:
            num = len(interval_neighbors[pid])
        else:
            num = enough_intervals[pid] + 1
        for i in range(num):
            neighbors.extend(interval_neighbors[pid][i])
            dists.extend(interval_dists[pid][i])
        neighbors, dists = np.array(neighbors), np.array(dists)
        sorted_flag = dists.argsort()
        num = min(len(neighbors), k)
        knn = neighbors[sorted_flag]
        dists = dists[sorted_flag]
        knn_dist[pid] = {'indexes': knn[: num], 'dists': dists[: num]}
    return knn_dist


def create_points_knn_by_road_step_neighbors(rid, match_dists, road_contains, road_neighbor,
                                             step_size, step_num, k, max_size):
    interval_neighbors, interval_dists, interval_counts = create_interval_neighbors_on_same_road(road_contains[rid],
                                                                                                 match_dists, step_size,
                                                                                                 step_num)
    enough_intervals, unenough = {i: -1 for i in road_contains[rid]}, [i for i in road_contains[rid]]
    for i in range(step_num - 1):
        create_interval_neighbors(unenough, i, step_size, road_neighbor, road_contains, match_dists,
                                  interval_neighbors, interval_dists, interval_counts, max_size)
        enough = assign_enough(i, interval_counts, unenough, k)
        for pid, interval in enough:
            enough_intervals[pid] = interval
        if len(unenough) == 0:
            break
    knn_dists = get_knn_from_interval_neighbors(interval_neighbors, interval_dists, enough_intervals, k)
    return knn_dists


class NetworkPointKNN:
    def __init__(self, edge_neighbors, max_dist):
        self._max_dist = max_dist
        self._edge_neighbors = edge_neighbors
        self._buffer_size = None
        self._step_size = None
        self._road_neighbors = None
        self._road_contains = None
        self._rids = None
        self._match_dists = None
        self._steps = None

    def _create_point_info(self, rids, match_dists):
        road_contains = [list() for i in range(len(self._edge_neighbors))]
        for i in range(len(rids)):
            road_contains[rids[i]].append(i)
        self._road_contains = road_contains
        self._rids = rids
        self._match_dists = match_dists

    def _create_points_knn(self, k):
        match_dists = self._match_dists
        road_neighbors = self._road_neighbors
        road_contains = self._road_contains
        step_size, step_num, max_size = self._step_size, self._steps, self._buffer_size
        indexes, dists = [list() for i in range(len(self._rids))], [list() for i in range(len(self._rids))]
        for rid in range(len(self._edge_neighbors)):
            if len(road_contains[rid]) > 0 and len(self._road_neighbors[rid]) > 0:
                knn_dists = create_points_knn_by_road_step_neighbors(rid, match_dists, road_contains,
                                                                     road_neighbors[rid], step_size,
                                                                     step_num, k, max_size)
                for pid in knn_dists.keys():
                    indexes[pid] = knn_dists[pid]['indexes']
                    dists[pid] = knn_dists[pid]['dists']
        return indexes, dists

    def _configure(self, rids, match_dists, buffer_size, step_size):
        self._buffer_size, self._step_size = buffer_size, step_size
        self._steps = int(buffer_size / step_size) + 1
        self._create_point_info(rids, match_dists)

    def _create_valid_road_neighbors(self):
        road_neighbors = list()
        edges_buffers = self._edge_neighbors
        road_contains = self._road_contains
        for i in range(len(edges_buffers)):
            edges_buffer = edges_buffers[i]
            del edges_buffer[i]
            road_neighbor = dict()
            for rid in edges_buffer.keys():
                if len(road_contains[rid]) > 0:
                    road_neighbor[rid] = edges_buffer[rid]['dist']
            road_neighbors.append(road_neighbor)
        self._road_neighbors = road_neighbors

    def run(self, rids, match_dists, buffer_size, step_size, k):
        self._configure(rids, match_dists, buffer_size, step_size)
        self._create_valid_road_neighbors()
        indexes, dists = self._create_points_knn(k)
        return indexes, dists
