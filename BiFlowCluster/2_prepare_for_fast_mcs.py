# -*- encoding=utf-8 -*-
#
# preparing for the fast Monte Carlo simulation method
from os import path

import pandas as pd

from Tool.load import load_road_neighbors
from Method.network_point_knn import NetworkPointKNN


def calc_keep_num(data_size, k):
    n = (1.5 * k * data_size) ** 0.5
    return int(n)


def knn_2_file(knn, dists, fn):
    wf = open(fn, 'w')
    for ns, dist in zip(knn, dists):
        dist = [round(i, 3) for i in dist]
        s = ['{0},{1}'.format(i, j) for i, j in zip(ns, dist)]
        wf.write(';'.join(s) + '\n')
    wf.close()


def run():
    root_dir = 'E:/Temp/'

    road_dir = path.join(root_dir, 'road network')
    data_dir = path.join(root_dir, 'synthetic data/70/TypeII')

    road_neighbor_file = path.join(road_dir, 'edges_neighbors.txt')
    data_file = path.join(data_dir, 'data.csv')

    road_neighbors = load_road_neighbors(road_neighbor_file)
    data = pd.read_csv(data_file)

    k = 15
    keep_num = calc_keep_num(len(data), k)
    point_type = ['o', 'd'][0]
    match_ids = data['{0}_rid'.format(point_type)].to_numpy()
    match_dists = data[['{0}_sd'.format(point_type), '{0}_ed'.format(point_type)]].to_numpy()

    npk = NetworkPointKNN(road_neighbors, 250)
    indexes, dists = npk.run(match_ids, match_dists, buffer_size=250, step_size=50, k=keep_num)

    keep_file = path.join(data_dir, 'keep_{0}_neighbors.txt'.format(point_type))
    knn_2_file(indexes, dists, keep_file)


if __name__ == '__main__':
    run()