# -*- coding=utf-8 -*-

import numpy as np
import random
from os import path


class flow(object):
    def __init__(self, origin, destination):
        self.origin = origin
        self.destination = destination
        self.BF_value = 0
        self.p_value = 0
        self.neighbor = list()
        self.sign = 0


# 读取流数据
def read_flows(ff, flows, sign):
    file = open(ff, 'r')

    if sign == 1:
        while True:
            str_line = file.readline()
            if not str_line:
                break

            str_lines = str_line.split(',')
            indices = '{}-{}'.format(int(str_lines[0]), int(str_lines[1]))
            if indices not in flows.keys():
                flows[indices] = [1, 0]
            else:
                flows[indices][0] += 1

    else:
        while True:
            str_line = file.readline()
            if not str_line:
                break

            # str_lines = str_line.split('\t')
            str_lines = str_line.split(',')
            indices = '{}-{}'.format(int(str_lines[0]), int(str_lines[1]))
            if indices not in flows.keys():
                flows[indices] = [0, 1]
            else:
                flows[indices][1] += 1
    file.close()


def read_cluster_region(cf):
    file = open(cf, 'r')
    cluster_region = list()

    while True:
        str_line = file.readline()
        if not str_line:
            break
        str_lines = str_line.split(',')
        new_region = list()
        for i in range(len(str_lines) - 1):
            new_region.append(int(str_lines[i]))
        cluster_region.append(new_region)

    return cluster_region


def read_cluster_flow(cluster_file):
    cluster_flow = list()
    file = open(cluster_file, 'r')

    for f in file:
        i = f.rstrip('\n')
        i = i.split(';')
        new_cluster = list()
        for flows in i:
            fl = flows.split(',')
            new_cluster.append(fl[0] + '-' + fl[1])
        cluster_flow.append(new_cluster)
    return cluster_flow


# def read_cluster_flow(cluster_region):
#     cluster_flow = [list() for i in range(int(len(cluster_region)/2))]
#     for i in range(int(len(cluster_region)/2)):
#         for l in cluster_region[2 * i]:
#             for n in cluster_region[2 * i + 1]:
#                 cluster_flow[i].append('{}-{}'.format(l, n))
#     return cluster_flow


# 读取流数据
def read_neighbor(nf):
    file = open(nf, 'r')
    spatial_neighbors = list()

    while True:
        str_line = file.readline()
        if not str_line:
            break
        str_lines = str_line.split(',')

        neighbor = list()
        if str_line != '\n':
            for i in range(len(str_lines)):
                neighbor.append(int(str_lines[i]))
        spatial_neighbors.append(neighbor)

    file.close()

    return spatial_neighbors


def read_simulation_flow(ff, value_type):
    file = open(ff, 'r')
    attribute_value = dict()

    if value_type == 1:
        while True:
            str_line = file.readline()
            if not str_line:
                break
            str_lines = str_line.split(',')

            indices = '{}-{}'.format(int(str_lines[0]), int(str_lines[1]))
            attribute_value[indices] = [float(str_lines[2]), float(str_lines[3])]

    else:
        while True:
            str_line = file.readline()
            if not str_line:
                break
            str_lines = str_line.split(',')

            indices = '{}-{}'.format(int(str_lines[0]), int(str_lines[1]))
            attribute_value[indices] = [float(str_lines[3]), float(str_lines[2])]

    file.close()

    return attribute_value


# 统计流数据并标准化属性值
def normal_flow_attribute_value(attribute_value):
    origin_value1 = list()
    origin_value2 = list()

    for l in attribute_value.keys():
        origin_value1.append(attribute_value[l][0])
        origin_value2.append(attribute_value[l][1])

    origin_value1 = np.array(origin_value1)
    origin_value2 = np.array(origin_value2)
    average1 = np.mean(origin_value1)
    std1 = np.std(origin_value1)
    average2 = np.mean(origin_value2)
    std2 = np.std(origin_value2)

    nor_value1 = (origin_value1 - average1) / std1
    nor_value2 = (origin_value2 - average2) / std2

    cir_num = 0
    for l in attribute_value.keys():
        attribute_value[l] = [nor_value1[cir_num], nor_value2[cir_num]]
        cir_num += 1


def search_flow_neighbors(origin, destination, attribute_value, besides_value, neighbor, new_flow):
    for l in neighbor[destination - 1]:
        index1 = '{}-{}'.format(origin, l)
        if index1 in attribute_value.keys():
            new_flow.neighbor.append(index1)
            besides_value.append(attribute_value[index1][1])

    for l in neighbor[origin - 1]:
        index1 = '{}-{}'.format(l, destination)
        if index1 in attribute_value.keys():
            new_flow.neighbor.append(index1)
            besides_value.append(attribute_value[index1][1])


def distinguish_sign(value_1, value_2):
    if value_1 > 0 and value_2 > 0:
        new_sign = 1
    elif value_2 < 0 < value_1:
        new_sign = 2
    elif value_1 < 0 < value_2:
        new_sign = 3
    elif value_1 < 0 and value_2 < 0:
        new_sign = 4

    return new_sign


def calculate_bf_value(attribute_value, neighbor, Sign):
    flows = list()
    flow_indices = dict()
    index_id = 0

    for l in attribute_value.keys():
        besides_value = list()
        besides_value.append(attribute_value[l][1])
        loc = l.split('-')
        new_flow = flow(int(loc[0]), int(loc[1]))
        search_flow_neighbors(new_flow.origin, new_flow.destination, attribute_value, besides_value, neighbor, new_flow)
        new_flow.BF_value = attribute_value[l][0] * np.sum(besides_value)
        new_sign = distinguish_sign(attribute_value[l][0], np.sum(besides_value))
        if new_sign == Sign:
            flows.append(new_flow)
            flow_indices[l] = index_id
            index_id += 1

    return flows, flow_indices


# 假设检验
def Monte_carlo_simulation(flows, attribute_value, rep_num, sig_level, Sign):
    sort_num = [i for i in range(len(attribute_value))]
    keys = [l for l in attribute_value.keys()]
    cluster_num = len(flows)
    p_value = [0 for i in range(cluster_num)]
    final_true_flows = list()

    for i in range(rep_num):
        # print(sort_num)
        random.shuffle(sort_num)
        for j in range(cluster_num):
            if p_value[j] <= sig_level * 2000:
                new_value = list()
                origin_index = '{}-{}'.format(flows[j].origin, flows[j].destination)
                location = keys.index(origin_index)
                new_value.append(attribute_value[keys[sort_num[location]]][1])
                for e in flows[j].neighbor:
                    # index1 = '{}-{}'.format(e[0], e[1])
                    location = keys.index(e)
                    # print(sort_num[location])
                    new_value.append(attribute_value[keys[sort_num[location]]][1])
                final_sum = np.sum(new_value)
                core_value = attribute_value['{}-{}'.format(flows[j].origin, flows[j].destination)][0]
                new_sign = distinguish_sign(core_value, final_sum)
                if new_sign == Sign:
                    new_stat = core_value * final_sum
                    if abs(new_stat) >= abs(flows[j].BF_value):
                        p_value[j] += 1
        # print('simulation')
        # print(time.process_time())

    for i in range(cluster_num):
        flows[i].p_value = (p_value[i]) / (rep_num + 1)
    final_flows = sorted(flows, key=lambda x: x.p_value)

    for i in range(cluster_num, 0, -1):
        if final_flows[i - 1].p_value <= sig_level:
            final_true_flows = final_flows[:i]
            break

    return final_true_flows


def delete_overlap_region(clusters):
    origin_region = list()
    destination_region = list()
    cluster_flow = list()

    for cluster in clusters:
        new_flow = list()
        new_origin = list()
        new_destination = list()
        for f in cluster:
            loc = f.split('-')
            new_flow.append('{} {}'.format(int(loc[0]), int(loc[1])))
            if loc[0] not in new_origin:
                new_origin.append(loc[0])
            if loc[1] not in new_destination:
                new_destination.append(loc[1])
        origin_region.append(new_origin)
        destination_region.append(new_destination)
        cluster_flow.append(new_flow)

    return origin_region, destination_region, cluster_flow


def iteration_neighbor(tar_flow, cluster, flow_indices, flows, final_flows):
    for nei in tar_flow.neighbor:
        if nei in flow_indices.keys():
            if flows[flow_indices[nei]] in final_flows:
                if flows[flow_indices[nei]].sign == 0:
                    cluster.append(
                        '{}-{}'.format(flows[flow_indices[nei]].origin, flows[flow_indices[nei]].destination))
                    flows[flow_indices[nei]].sign = 1
                    iteration_neighbor(flows[flow_indices[nei]], cluster, flow_indices, flows, final_flows)


def construct_cluster(flows, final_flows, flow_indices):
    cluster = list()

    for l in final_flows:
        if l.sign == 0:
            l.sign = 1
            new_cluster = list()
            new_cluster.append('{}-{}'.format(l.origin, l.destination))
            iteration_neighbor(l, new_cluster, flow_indices, flows, final_flows)
            cluster.append(new_cluster)

    return cluster


def construct_true_cluster(flows, final_flows, flow_indices):
    cluster = list()
    cluster_num = 0

    for l in final_flows:
        if l.sign == 0:
            candidate_num = search_neighbor_cluster(l, flow_indices, flows, final_flows)
            if candidate_num != 0:
                l.sign = candidate_num
                cluster[candidate_num - 1].append('{}-{}'.format(l.origin, l.destination))
            else:
                cluster_num += 1
                l.sign = cluster_num
                new_cluster = list()
                new_cluster.append('{}-{}'.format(l.origin, l.destination))
                cluster.append(new_cluster)
    final_cluster = sorted(cluster, key=lambda x: len(x), reverse=True)

    return final_cluster


def search_neighbor_cluster(tar_flow, flow_indices, flows, final_flows):
    for nei in tar_flow.neighbor:
        if nei in flow_indices.keys():
            if flows[flow_indices[nei]] in final_flows:
                if flows[flow_indices[nei]].sign != 0:
                    return flows[flow_indices[nei]].sign
    return 0


# 存储显著类
# def save_cluster(final_true_flow, cf):
#     file = open(cf, 'w')
#
#     for l in final_true_flow:
#         file.writelines('{},'.format(l.origin))
#
#         file.writelines('{}'.format(l.destination))
#
#         file.writelines(';{};'.format(l.BF_value))
#
#         neighbor_flow = ['{} {}'.format(j[0], j[1]) for j in l.neighbor]
#         s_flow = ','.join(neighbor_flow)
#         file.writelines(s_flow + ';')
#
#         file.writelines('{}\n'.format(l.p_value))
#
#     file.close()


def save_cluster(origin_region, destination_region, cluster_flow, cf):
    file = open(cf, 'w')

    for i in range(len(origin_region)):
        origins = ['{}'.format(j) for j in origin_region[i]]
        s_origin = ','.join(origins)
        file.writelines(s_origin + ';')

        destinations = ['{}'.format(j) for j in destination_region[i]]
        s_destination = ','.join(destinations)
        file.writelines(s_destination + ';')

        flows = ','.join(cluster_flow[i])
        file.writelines(flows + '\n')


def run():
    rep_num = 999
    sig_level = 0.01
    print("Result for data sets contain a single cluster")
    cluster_name = ['H-H Cluster:', 'H-L Cluster:', 'L-H Cluster:', 'L-L Cluster:']
    for i in range(4):
        data_dir = 'C:/Users/12449/Desktop/BiflowAMOEBA/IJGIS-code&data(2021)/Synthetic data/D{}'.format(i + 1)
        Sign = i + 1
        nf = path.join(data_dir, 'neighbors.txt')
        cf = path.join(data_dir, 'Cluster_flow.txt')
        ff = path.join(data_dir, 'flow.txt')

        cluster_flow = read_cluster_flow(cf)
        neighbor = read_neighbor(nf)
        TP = 0
        final_num = 0
        cluster = cluster_flow[0]
        for j in range(2):
            if j == 0:
                Sign = Sign
            elif j == 1 and Sign != 1:
                Sign = 3 - Sign + 2
            attribute_value = read_simulation_flow(ff, j + 1)
            normal_flow_attribute_value(attribute_value)

            flows, flow_indices = calculate_bf_value(attribute_value, neighbor, Sign)
            final_flows = Monte_carlo_simulation(flows, attribute_value, rep_num, sig_level, Sign)

            final_cluster = construct_true_cluster(flows, final_flows, flow_indices)

            for l in final_cluster:
                final_num += len(l)
                for n in l:
                    if n in cluster:
                        TP += 1
        if TP != 0:
            Precision = TP / final_num
            Recall = TP / (2 * len(cluster))
            F1 = 2 * Precision * Recall / (Precision + Recall)
        else:
            Precision = 0
            Recall = 0
            F1 = 0
        print(cluster_name[i])
        print(Precision)
        print(Recall)
        print(F1)

    print("----------------------------------------------")
    print("Result for data set contains multiple clusters")

    data_dir = 'C:/Users/12449/Desktop/BiflowAMOEBA/IJGIS-code&data(2021)/Synthetic data/D5'
    nf = path.join(data_dir, 'neighbors.txt')
    cf = path.join(data_dir, 'Cluster_flow.txt')
    ff = path.join(data_dir, 'flow.txt')

    cluster_flow = read_cluster_flow(cf)
    neighbor = read_neighbor(nf)
    for i in range(4):
        Sign = i + 1
        TP = 0
        final_num = 0
        cluster = cluster_flow[i]
        for j in range(2):
            if j == 0:
                Sign = Sign
            elif j == 1 and Sign != 1:
                Sign = 3 - Sign + 2
            attribute_value = read_simulation_flow(ff, j+1)
            normal_flow_attribute_value(attribute_value)

            flows, flow_indices = calculate_bf_value(attribute_value, neighbor, Sign)
            final_flows = Monte_carlo_simulation(flows, attribute_value, rep_num, sig_level, Sign)
            print(len(final_flows))

            final_cluster = construct_true_cluster(flows, final_flows, flow_indices)

            for l in final_cluster:
                final_num += len(l)
                for n in l:
                    if n in cluster:
                        TP += 1
        if TP != 0:
            Precision = TP / final_num
            Recall = TP / (2 * len(cluster))
            F1 = 2 * Precision * Recall / (Precision + Recall)
        else:
            Precision = 0
            Recall = 0
            F1 = 0
        print(cluster_name[i])
        print(Precision)
        print(Recall)
        print(F1)


if __name__ == '__main__':
    run()
