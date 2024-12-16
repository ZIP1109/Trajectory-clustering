# -*- encoding=utf-8 -*-
#
# identifying core flows using the fast Monte Carlo simulation method
from os import path

import numpy as np
from numpy import array


def parse_line(line, data_type):
    line = line.rstrip('\n')
    line = line.split(';')
    point_dict = dict()
    num = len(line)
    if num == 1 and len(line[0]) == 0:
        return point_dict
    for i in range(num):
        if data_type == 'point':
            pid, dist = line[i].split(',')
            nid = int(pid)
            point_dict[nid] = float(dist)
        elif data_type == 'flow':
            pid, dist, od, ed = line[i].split(',')
            nid = int(pid)
            point_dict[nid] = (float(dist), float(od), float(ed))
        else:
            raise ValueError('data_type must be "point" or "flow"')
    return point_dict


def load_neighbors_with_dist(nf, data_type):
    neighbors = list()
    f = open(nf)
    for i in f:
        one_neighbor_dict = parse_line(i, data_type)
        neighbors.append(one_neighbor_dict)
    return neighbors


def calc_densities(flow_knn, k):
    neighbor_sets = [set(neighbor.keys()) for neighbor in flow_knn]
    neighbors_size = list()
    denominator = k * (k - 1)
    densities = list()
    for i in range(len(flow_knn)):
        neighbors_size.append(len(neighbor_sets[i]))
        density = 0
        for n in neighbor_sets[i]:
            common = neighbor_sets[i].intersection(neighbor_sets[n])
            density += len(common)
        density = density / denominator
        densities.append(density)
    return densities, neighbors_size


def flows_2_dict(flows):
    flows_start_dict = dict()
    flows_end_dict = dict()
    for i in range(len(flows)):
        s, e = flows[i]
        flows_start_dict[s] = e
        flows_end_dict[e] = {'start': s, 'id': i}
    return flows_start_dict, flows_end_dict


def create_simulate_flow_info(indexes):
    np.random.shuffle(indexes)
    flows = list()
    flows_start_dict, flows_end_dict = dict(), dict()
    simulated_flows = dict()
    for i in range(len(indexes)):
        s, e = i, indexes[i]
        flows.append((s, e))
        flows_start_dict[s] = e
        flows_end_dict[e] = {'start': s, 'id': i}
    simulated_flows['flows'] = flows
    simulated_flows['flows_start_dict'] = flows_start_dict
    simulated_flows['flows_end_dict'] = flows_end_dict
    return simulated_flows


def get_common_flow(o_dict, d_dict, flows_start_dict, flows_end_dict):
    end_ids = set([flows_start_dict[i] for i in o_dict.keys()])
    d_neighbors = set(d_dict.keys())

    common_ends = end_ids.intersection(d_neighbors)

    flows, o_dists, d_dists = list(), list(), list()
    for end in common_ends:
        start = flows_end_dict[end]['start']
        flow_id = flows_end_dict[end]['id']
        flows.append(flow_id)
        o_dists.append(o_dict[start])
        d_dists.append(d_dict[end])
    dists = (array(o_dists) + array(d_dists)) / 2
    flows = array(flows)
    sorted_flag = dists.argsort()
    return flows[sorted_flag], dists[sorted_flag]


def create_simulated_knn(simulated_flows, o_neighbors, d_neighbors, k):
    knn, dists = list(), list()
    flows, flows_start_dict, flows_end_dict = simulated_flows['flows'], simulated_flows['flows_start_dict'], \
                                              simulated_flows['flows_end_dict']
    for i in range(len(flows)):
        start, end = flows[i]
        one_flows, one_dists = get_common_flow(o_neighbors[start], d_neighbors[end],
                                               flows_start_dict, flows_end_dict)
        num = min(len(one_flows), k)
        knn.append(one_flows[: num])
        dists.append(one_dists[: num])
    return knn, dists


def create_observed_simulated_knn(simulated_flows, o_neighbors, d_neighbors, k):
    knn, dists = list(), list()
    flows, flows_start_dict, flows_end_dict = simulated_flows['flows'], simulated_flows['flows_start_dict'], \
                                              simulated_flows['flows_end_dict']
    for i in range(len(flows)):
        one_flows, one_dists = get_common_flow(o_neighbors[i], d_neighbors[i],
                                               flows_start_dict, flows_end_dict)
        num = min(len(one_flows), k)
        knn.append(one_flows[: num])
        dists.append(one_dists[: num])
    return knn, dists


def get_shared_num(set1, set2):
    return len(set1.intersection(set2))


def correction_share_count(knn, simulated_flows, simulated_dists, o_neighbor, d_neighbor, k):
    count = 0
    for n in knn:
        if len(simulated_dists[n]) == k:
            start, end = simulated_flows[n]
            dist = round((o_neighbor[start] + d_neighbor[end]) / 2, 3)
            if dist < simulated_dists[n][-1]:
                count += 1
    return count


def calc_random_densities(observed_simulated_knn, simulated_flows, simulated_knn, simulated_dists,
                          o_neighbors, d_neighbors, k):
    densities = list()
    for i in range(len(observed_simulated_knn)):
        total_share_count = 0
        set_knn = set(observed_simulated_knn[i])
        for n in set_knn:
            total_share_count += get_shared_num(set_knn, set(simulated_knn[n]))
        cc = correction_share_count(set_knn, simulated_flows, simulated_dists, o_neighbors[i], d_neighbors[i], k)
        density = (total_share_count - cc) / (k * (k - 1))
        densities.append(density)
    return np.array(densities)


def correct_p_values(flow_knn, p_values, k):
    for i in range(len(p_values)):
        if len(flow_knn[i].keys()) < k:
            p_values[i] = 1.0
    return p_values


def fast_mcs(flow_knn, o_neighbors, d_neighbors, k):
    densities, neighbors_size = calc_densities(flow_knn, k)
    size = len(densities)
    indexes = np.arange(size)
    p_values = np.ones(size)
    for i in range(99):
        simulated_flows = create_simulate_flow_info(indexes)
        simulated_knn, simulated_dists = create_simulated_knn(simulated_flows, o_neighbors, d_neighbors, k)
        observed_simulated_knn, _ = create_observed_simulated_knn(simulated_flows, o_neighbors, d_neighbors, k)
        random_densities = calc_random_densities(observed_simulated_knn, simulated_flows['flows'], simulated_knn,
                                                 simulated_dists, o_neighbors, d_neighbors, k)
        p_values[random_densities >= densities] += 1
    p_values = p_values / 100
    return correct_p_values(flow_knn, p_values, k)


def _calc_core_range(index, knn_dict, knn_s, min_num):
    core_range = {'sn': list()}
    max_os, max_ds, max_fs = 0, 0, 0
    for n in knn_dict.keys():
        common = knn_s[index].intersection(knn_s[n])
        if len(common) > min_num:
            core_range['sn'].append(n)
            if max_fs < knn_dict[n][0]:
                max_fs = knn_dict[n][0]
            if max_os < knn_dict[n][1]:
                max_os = knn_dict[n][1]
            if max_ds < knn_dict[n][2]:
                max_ds = knn_dict[n][2]
    core_range['range'] = (max_fs, max_os, max_ds)
    return core_range


def get_cores_range(flow_knn, k, p_values, sig_level):
    cores_range = dict()
    min_num = k / 2
    knn_sets = [set(neighbor.keys()) for neighbor in flow_knn]
    for i in range(len(flow_knn)):
        if p_values[i] <= sig_level:
            cores_range[i] = _calc_core_range(i, flow_knn[i], knn_sets, min_num)
    return cores_range


def cores_range_2_file(cores_range, fn):
    wf = open(fn, 'w')
    for core_index in cores_range.keys():
        neighbor_str = ','.join([str(i) for i in cores_range[core_index]['sn']])
        range_str = ','.join([str(i) for i in cores_range[core_index]['range']])
        wf.write('{0};{1};{2}\n'.format(core_index, neighbor_str, range_str))
    wf.close()


def run():
    root_dir = 'E:/Temp/'

    data_dir = path.join(root_dir, 'synthetic data/70/TypeI')

    flow_knn_file = path.join(data_dir, 'knn_dist.txt')
    o_neighbors_file = path.join(data_dir, 'keep_o_neighbors.txt')
    d_neighbors_file = path.join(data_dir, 'keep_d_neighbors.txt')
    range_file = path.join(data_dir, 'cores_range.txt')

    k, sig_level = 15, 0.05
    flow_knn = load_neighbors_with_dist(flow_knn_file, 'flow')
    o_neighbors = load_neighbors_with_dist(o_neighbors_file, 'point')
    d_neighbors = load_neighbors_with_dist(d_neighbors_file, 'point')

    p_values = fast_mcs(flow_knn, o_neighbors, d_neighbors, k)
    cores_range = get_cores_range(flow_knn, k, p_values, sig_level)
    cores_range_2_file(cores_range, range_file)


if __name__ == '__main__':
    run()