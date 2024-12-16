# -*- encoding=utf-8 -*-
#
# constructing the bivariate flow clusters by extending the concept of density connectivity.
from os import path

import pandas as pd
import numpy as np
import networkx as nx


def load_same_neighbors(nf):
    same_neighbors = dict()
    rf = open(nf)
    for i in rf:
        i = i.rstrip('\n')
        i = i.split(';')
        cid = int(i[0])
        ns = [int(j) for j in i[1].split(',')]
        same_neighbors[cid] = ns
    return same_neighbors


def load_other_neighbors(nf):
    other_neighbors = dict()
    rf = open(nf)
    for i in rf:
        i = i.rstrip('\n')
        i = i.split(';')
        if len(i[1]) == 0:
            other_neighbors[int(i[0])] = list()
        else:
            other_neighbors[int(i[0])] = [int(j) for j in i[1].split(',')]
    return other_neighbors


def load_one_type_data(data_dir, data_file, same_neighbor_file, other_neighbor_file):
    data_file = path.join(data_dir, data_file)
    same_neighbor_file = path.join(data_dir, same_neighbor_file)
    other_neighbor_file = path.join(data_dir, other_neighbor_file)

    data = pd.read_csv(data_file)
    same_neighbors = load_same_neighbors(same_neighbor_file)
    other_neighbors = load_other_neighbors(other_neighbor_file)
    return data, same_neighbors, other_neighbors


def set_flow_type(flow_num, neighbors):
    flow_type = np.zeros(flow_num, dtype=int)
    flow_type[list(neighbors.keys())] = 1
    for cid in neighbors.keys():
        for n in neighbors[cid]:
            if not flow_type[n] == 1:
                flow_type[n] = 2
    return flow_type


def _create_core_flag(sizeI, sizeII, neighborsI, neighborsII):
    core_flag = np.zeros(sizeI + sizeII, dtype=int)
    core_flag[list(neighborsI.keys())] = 1
    for i in neighborsII.keys():
        core_flag[i + sizeI] = 1
    return core_flag


def _overlap_neighbors(core_neighbors):
    for cid in core_neighbors:
        for on in core_neighbors[cid]['o_core_n']:
            if cid not in core_neighbors[on]['o_core_n']:
                core_neighbors[on]['o_core_n'].append(cid)


def create_core_neighbors(type_I_size, type_II_size, same_neighborsI, other_neighborsI,
                          same_neighborsII, other_neighborsII):
    core_flag = _create_core_flag(type_I_size, type_II_size, same_neighborsI, same_neighborsII)
    core_neighbors = dict()
    for cid in same_neighborsI.keys():
        core_neighbor = {'s_core_n': list(), 's_non_core_n': list(),
                         'o_core_n': list(), 'o_non_core_n': list()}
        for sn in same_neighborsI[cid]:
            if core_flag[sn] == 1:
                core_neighbor['s_core_n'].append(sn)
            else:
                core_neighbor['s_non_core_n'].append(sn)
                core_flag[sn] = 2
        for on in other_neighborsI[cid]:
            on += type_I_size
            if core_flag[on] == 1:
                core_neighbor['o_core_n'].append(on)
            else:
                core_neighbor['o_non_core_n'].append(on)
        core_neighbors[cid] = core_neighbor

    for cid in same_neighborsII.keys():
        core_neighbor = {'s_core_n': list(), 's_non_core_n': list(),
                         'o_core_n': list(), 'o_non_core_n': list()}
        for sn in same_neighborsII[cid]:
            sn += type_I_size
            if core_flag[sn] == 1:
                core_neighbor['s_core_n'].append(sn)
            else:
                core_neighbor['s_non_core_n'].append(sn)
                core_flag[sn] = 2
        for on in other_neighborsII[cid]:
            if core_flag[on] == 1:
                core_neighbor['o_core_n'].append(on)
            else:
                core_neighbor['o_non_core_n'].append(on)
        core_neighbors[cid + type_I_size] = core_neighbor
    _overlap_neighbors(core_neighbors)
    return core_flag, core_neighbors


def identify_core_types(core_neighbors, size):
    core_types = dict()
    for cid in core_neighbors.keys():
        if len(core_neighbors[cid]['o_core_n']) > 0:
            core_types[cid] = 'HH'
        else:
            if cid >= size:
                core_types[cid] = 'LH'
            else:
                core_types[cid] = 'HL'
    return core_types


def _get_clusters_from_graph(graph, min_num=15):
    clusters = list()
    for c in nx.connected_components(graph):
        nodes = graph.subgraph(c).nodes()
        if len(nodes) >= min_num:
            clusters.append([i for i in nodes])
    clusters = sorted(clusters, key=lambda a: len(a), reverse=True)
    return clusters


def _create_core_graphs(core_neighbors, size):
    core_types = identify_core_types(core_neighbors, size)
    h_h_g, h_l_g, l_h_g = nx.Graph(), nx.Graph(), nx.Graph()
    for cid in core_types.keys():
        if core_types[cid] == 'HH':
            for on in core_neighbors[cid]['o_core_n']:
                h_h_g.add_edge(cid, on)
            for sn in core_neighbors[cid]['s_core_n']:
                if core_types[sn] == 'HH':
                    h_h_g.add_edge(cid, sn)
        elif core_types[cid] == 'HL':
            for sn in core_neighbors[cid]['s_core_n']:
                if core_types[sn] == 'HL':
                    h_l_g.add_edge(cid, sn)
        else:
            for sn in core_neighbors[cid]['s_core_n']:
                if core_types[sn] == 'LH':
                    l_h_g.add_edge(cid, sn)
    return h_h_g, h_l_g, l_h_g


def _add_no_cores_2_cluster(clusters, core_neighbors, core_flag):
    cluster_flag = np.zeros(len(core_flag), dtype=bool)
    for ct in clusters.keys():
        for cluster in clusters[ct]:
            to_add = list()
            if ct == 'HH':
                for cid in cluster:
                    for sn in core_neighbors[cid]['s_non_core_n']:
                        if not cluster_flag[sn]:
                            cluster_flag[sn] = True
                            to_add.append(sn)
            else:
                for cid in cluster:
                    for on in core_neighbors[cid]['o_non_core_n']:
                        if core_flag[on] == 0 and not cluster_flag[on]:
                            to_add.append(on)
                            cluster_flag[on] = True
            cluster.extend(to_add)
    return clusters


def create_bi_clusters(core_neighbors, size, core_flag):
    clusters = dict()
    h_h_g, h_l_g, l_h_g = _create_core_graphs(core_neighbors, size)
    clusters['HH'] = _get_clusters_from_graph(h_h_g, min_num=10)
    clusters['HL'] = _get_clusters_from_graph(h_l_g, min_num=10)
    clusters['LH'] = _get_clusters_from_graph(l_h_g, min_num=10)
    _add_no_cores_2_cluster(clusters, core_neighbors, core_flag)
    return clusters


def parse_cluster(cluster, size):
    parsed_cluster = {'I': list(), 'II': list()}
    for cid in cluster:
        if cid >= size:
            parsed_cluster['II'].append(cid - size)
        else:
            parsed_cluster['I'].append(cid)
    return parsed_cluster


def _cluster_2_file(cf, parsed_cluster):
    wf = open(cf, 'w')
    for i in ['I', 'II']:
        if len(parsed_cluster[i]) > 0:
            members = sorted(parsed_cluster[i])
        else:
            members = list()
        s = [str(j) for j in members]
        wf.write('{0}:{1}\n'.format(i, ','.join(s)))
    wf.close()


def clusters_2_file(result_dir, clusters, size):
    for ct in clusters.keys():
        for i, cluster in enumerate(clusters[ct]):
            cf = path.join(result_dir, '{0}_{1}.txt'.format(ct, i + 1))
            parsed_cluster = parse_cluster(cluster, size)
            _cluster_2_file(cf, parsed_cluster)


def run():
    root_dir = 'E:/Temp/'

    data_dir = path.join(root_dir, 'synthetic data/70/Type')

    dataI, same_neighborsI, other_neighborsI = load_one_type_data(data_dir + 'I', 'data.csv',
                                                                  'cores_range.txt', 'other_neighbors.txt')
    dataII, same_neighborsII, other_neighborsII = load_one_type_data(data_dir + 'II', 'data.csv',
                                                                     'cores_range.txt', 'other_neighbors.txt')
    core_flag, core_neighbors = create_core_neighbors(len(dataI), len(dataII), same_neighborsI,
                                                      other_neighborsI, same_neighborsII, other_neighborsII)
    clusters = create_bi_clusters(core_neighbors, len(dataI), core_flag)
    clusters_2_file(data_dir.replace('Type', ''), clusters, len(dataI))


if __name__ == '__main__':
    run()

