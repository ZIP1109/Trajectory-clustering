# -*- encoding=utf-8 -*-
#
# searching the other type flow neighbors
from os import path

import pandas as pd

from Tool.load import load_road_neighbors


def load_range_file(range_file):
    cores_range = dict()
    rf = open(range_file)
    for i in rf:
        i = i.rstrip('\n')
        i = i.split(';')
        cid = int(i[0])
        ns = [int(j) for j in i[1].split(',')]
        ds = [float(j) for j in i[2].split(',')]
        cores_range[cid] = {'ns': ns, 'range': ds}
    return cores_range


def load_one_type_data(data_dir, data_file, range_file):
    data_file = path.join(data_dir, data_file)
    range_file = path.join(data_dir, range_file)

    data = pd.read_csv(data_file)
    cores_range = load_range_file(range_file)
    return data, cores_range


def stat_road_contains(road_number, data):
    road_contains = {i: {'o': list(), 'd': list()} for i in range(road_number)}
    for i in range(len(data)):
        o_rid, d_rid = data[['o_rid', 'd_rid']].iloc[i]
        road_contains[o_rid]['o'].append(i)
        road_contains[d_rid]['d'].append(i)
    return road_contains


def filter_road_neighbors(dist, road_neighbor):
    new_road_neighbor = {rid: road_neighbor[rid] for rid in road_neighbor.keys() if
                         min(road_neighbor[rid]['dist']) < dist}
    return new_road_neighbor


def _find_road_neighbor_contains(road_neighbors, point_type,road_contains):
    pid = list()
    for rid in road_neighbors.keys():
        pid.extend(road_contains[rid][point_type])
    return pid


def _find_coarse_neighbors(od_rid, core_range, road_contains, roads_neighbors):
    o_rid, d_rid = od_rid
    fd, od, dd = core_range
    o_road_neighbors = filter_road_neighbors(od, roads_neighbors[o_rid])
    d_road_neighbors = filter_road_neighbors(dd, roads_neighbors[d_rid])
    o_neighbors = _find_road_neighbor_contains(o_road_neighbors, 'o', road_contains)
    d_neighbors = _find_road_neighbor_contains(d_road_neighbors, 'd', road_contains)
    return list(set(o_neighbors).intersection(set(d_neighbors)))


def _calc_shortest_path_dist(rid, sd, ed, rid1, sd1, ed1, road_neighbors):
    if rid == rid1:
        return abs(sd - sd1)
    d1 = sd + road_neighbors[rid1]['ss'] + sd1
    d2 = sd + road_neighbors[rid1]['se'] + ed1
    d3 = ed + road_neighbors[rid1]['es'] + sd1
    d4 = ed + road_neighbors[rid1]['es'] + ed1
    return min([d1, d2, d3, d4])


def _filter_flow_neighbors_by_dist(flow_info, core_range, near_flows, dataII, roads_neighbors):
    o_rid, o_sd, o_ed = flow_info[['o_rid', 'o_sd', 'o_ed']]
    d_rid, d_sd, d_ed = flow_info[['d_rid', 'd_sd', 'd_ed']]
    dist_neighbors = dict()
    for fid in near_flows:
        oo_rid, oo_sd, oo_ed = dataII.iloc[fid][['o_rid', 'o_sd', 'o_ed']]
        o_dist = _calc_shortest_path_dist(o_rid, o_sd, o_ed, oo_rid, oo_sd, oo_ed, roads_neighbors[o_rid])
        if o_dist > core_range[1]:
            continue
        dd_rid, dd_sd, dd_ed = dataII.iloc[fid][['d_rid', 'd_sd', 'd_ed']]
        d_dist = _calc_shortest_path_dist(d_rid, d_sd, d_ed, dd_rid, dd_sd, dd_ed, roads_neighbors[d_rid])
        if d_dist > core_range[2]:
            continue
        f_dist = round((o_dist + d_dist) / 2, 3)
        if f_dist < core_range[0]:
            dist_neighbors[fid] = (f_dist, round(o_dist, 3), round(d_dist, 3))
    return dist_neighbors


def _search_oof_4_one(flow_info, core_range, road_contains, dataII, roads_neighbors):
    near_flows = _find_coarse_neighbors(flow_info[['o_rid', 'd_rid']], core_range, road_contains, roads_neighbors)
    dist_neighbors = _filter_flow_neighbors_by_dist(flow_info, core_range, near_flows, dataII, roads_neighbors)
    return dist_neighbors


def search_cross_type_neighbors(rangeI, dataI, dataII, roads_neighbors):
    road_contains = stat_road_contains(len(roads_neighbors), dataII)
    core_neighbors = dict()
    for cid in rangeI.keys():
        flow_info = dataI.iloc[cid]
        core_neighbors[cid] = _search_oof_4_one(flow_info, rangeI[cid]['range'], road_contains, dataII, roads_neighbors)
    return core_neighbors


def other_neighbors_2_file(other_neighbors, fn):
    wf = open(fn, 'w')
    for cid in other_neighbors.keys():
        s = ','.join([str(i) for i in other_neighbors[cid].keys()])
        wf.write('{0};{1}\n'.format(cid, s))
    wf.close()


def run():
    root_dir = 'E:/Temp/'

    road_dir = path.join(root_dir, 'road network')
    data_dir = path.join(root_dir, 'synthetic data/70/Type')

    road_neighbor_file = path.join(road_dir, 'edges_neighbors.txt')

    roads_neighbors = load_road_neighbors(road_neighbor_file)
    dataI, rangeI = load_one_type_data(data_dir + 'I', 'data.csv', 'cores_range.txt')
    dataII, rangeII = load_one_type_data(data_dir + 'II', 'data.csv', 'cores_range.txt')

    other_neighborsI = search_cross_type_neighbors(rangeI, dataI, dataII, roads_neighbors)
    other_neighborsII = search_cross_type_neighbors(rangeII, dataII, dataI, roads_neighbors)

    fn1 = path.join(data_dir + 'I', 'other_neighbors.txt')
    other_neighbors_2_file(other_neighborsI, fn1)
    fn2 = path.join(data_dir + 'II', 'other_neighbors.txt')
    other_neighbors_2_file(other_neighborsII, fn2)


if __name__ == '__main__':
    run()