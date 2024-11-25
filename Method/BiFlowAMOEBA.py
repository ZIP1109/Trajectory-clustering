# -*- coding=utf-8 -*-

import numpy as np
from math import sqrt
import random
from os import path


# 冷热点类
class Hot_Spot_region(object):
    def __init__(self, origin, destination):
        self.member_o = [origin]
        self.member_d = [destination]
        self.member_flow = list()
        self.member_value1 = list()
        self.member_value2 = list()
        self.stat = 0
        self.next_stat = 0
        self.neighbor_flow = dict()
        self.neigh_flow_flag = dict()
        self.sign = 1
        self.p_value = 0


# 扩展邻居流类
class Flow(object):
    def __init__(self):
        self.flow = list()
        self.Zx = list()
        self.Zy = list()
        self.value = 0


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
        if str_lines[0] != '\n':
            for i in range(len(str_lines)):
                neighbor.append(int(str_lines[i]))
        spatial_neighbors.append(neighbor)

    return spatial_neighbors


def read_cluster_flow(cluster_file):
    cluster_flow = list()
    file = open(cluster_file, 'r')

    for f in file:
        i = f.rstrip('\n')
        i = i.split(';')
        new_cluster = list()
        for flow in i:
            fl = flow.split(',')
            new_cluster.append(fl[0] + '-' + fl[1])
        cluster_flow.append(new_cluster)
    return cluster_flow


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
    print(len(flows.keys()))
    file.close()


def read_simulation_flow(ff, flows):
    file = open(ff, 'r')

    while True:
        str_line = file.readline()
        if not str_line:
            break
        str_lines = str_line.split(',')

        indices = '{}-{}'.format(int(str_lines[0]), int(str_lines[1]))
        flows[indices] = [float(str_lines[2]), float(str_lines[3])]

    file.close()


# 统计流数据并标准化属性值
def normal_flow_attribute_value(attribute_value):
    unique_mode = list()
    origin_value1 = list()
    origin_value2 = list()

    for l in attribute_value.keys():
        str_lines = l.split('-')
        unique_mode.append([int(str_lines[0]), int(str_lines[1])])
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

    return unique_mode


# 寻找种子流
def find_seed_flows(G_sign, unique_mode, attribute_value):
    seed_flow = list()
    seed_flow_value = dict()

    if G_sign == 1:
        for l in unique_mode:
            index1 = '{}-{}'.format(l[0], l[1])
            value = attribute_value[index1]
            if value[0] > 0 and value[1] > 0:
                seed_flow.append(l)
                seed_flow_value[index1] = [value[0], value[1]]

    elif G_sign == 2:
        for l in unique_mode:
            index1 = '{}-{}'.format(l[0], l[1])
            value = attribute_value[index1]
            if value[0] > 0 > value[1]:
                seed_flow.append(l)
                seed_flow_value[index1] = [value[0], value[1]]

    elif G_sign == 3:
        for l in unique_mode:
            index1 = '{}-{}'.format(l[0], l[1])
            value = attribute_value[index1]
            if value[1] > 0 > value[0]:
                seed_flow.append(l)
                seed_flow_value[index1] = [value[0], value[1]]

    elif G_sign == 4:
        for l in unique_mode:
            index1 = '{}-{}'.format(l[0], l[1])
            value = attribute_value[index1]
            if value[0] < 0 and value[1] < 0:
                seed_flow.append(l)
                seed_flow_value[index1] = [value[0], value[1]]

    return seed_flow, seed_flow_value


# 标识流
def cal_flow_flag(seed_flow, region_flag, G_sign, attribution_value):

    if G_sign == 1:
        for l in seed_flow:
            index1 = '{}-{}'.format(l[0], l[1])
            if attribution_value[index1][0] > 0 or attribution_value[index1][1] > 0:
                region_flag[l[0] - 1, l[1] - 1] = 1

    elif G_sign == 2:
        for l in seed_flow:
            index1 = '{}-{}'.format(l[0], l[1])
            if attribution_value[index1][0] > 0 or attribution_value[index1][1] < 0:
                region_flag[l[0] - 1, l[1] - 1] = 1

    elif G_sign == 3:
        for l in seed_flow:
            index1 = '{}-{}'.format(l[0], l[1])
            if attribution_value[index1][0] < 0 or attribution_value[index1][1] > 0:
                region_flag[l[0] - 1, l[1] - 1] = 1

    else:
        for l in seed_flow:
            index1 = '{}-{}'.format(l[0], l[1])
            if attribution_value[index1][0] < 0 or attribution_value[index1][1] < 0:
                region_flag[l[0] - 1, l[1] - 1] = 1


def initial_flow_neighbor(unique_mode, neighbor, region_flag):
    flow_neighbor = dict()

    for l in unique_mode:
        if region_flag[l[0] - 1, l[1] - 1] == 1:
            neighbor_o = neighbor[l[0] - 1]
            neighbor_d = neighbor[l[1] - 1]
            neighbor_o.append(l[0])
            neighbor_d.append(l[1])
            new_neighbor = list()
            for i in neighbor_o:
                for j in neighbor_d:
                    if region_flag[i - 1, j - 1] == 1:
                        new_neighbor.append('{}-{}'.format(i, j))
            index1 = '{}-{}'.format(l[0], l[1])
            flow_neighbor[index1] = new_neighbor

    return flow_neighbor


# 寻找流邻居
def initial_cluster_neighbor(cluster_flag, attribution_value, index_matrix,
                             type_num, G_sign, all_cluster, neighbor_flow):
    new_cluster = all_cluster[cluster_flag - 1]

    index1 = '{}-{}'.format(new_cluster.member_o[0], new_cluster.member_d[0])
    for nei in neighbor_flow[index1]:
        new_flow = nei.split('-')
        if nei not in new_cluster.neigh_flow_flag.keys():
            new_cluster.neigh_flow_flag[nei] = 1
            neighbor_flag = index_matrix[int(new_flow[0]) - 1, int(new_flow[1]) - 1]
            if neighbor_flag != 0 and neighbor_flag > cluster_flag:
                neighbor_cluster = all_cluster[neighbor_flag - 1]
                neighbor_cluster.neigh_flow_flag[index1] = 1
                k_value = len(new_cluster.member_flow) + len(neighbor_cluster.member_flow)
                new_Zx = (sum(new_cluster.member_value1) + sum(neighbor_cluster.member_value1)) / k_value
                new_Zy = (sum(new_cluster.member_value2) + sum(neighbor_cluster.member_value2)) / k_value
                new_sign = distinguish_sign(new_Zx, new_Zy)

                if new_sign == G_sign:
                    new_value = sqrt((type_num - 1) / (type_num / k_value - 1)) * new_Zx * new_Zy
                    old_value = sqrt((type_num - 1) / (type_num / len(neighbor_cluster.member_value1) - 1)) \
                                * np.mean(neighbor_cluster.member_value1) * \
                                np.mean(neighbor_cluster.member_value2)

                    if abs(new_value) > abs(new_cluster.stat) and abs(new_value) > abs(old_value):
                        new_cluster.neighbor_flow[neighbor_flag] = new_value
                        neighbor_cluster.neighbor_flow[cluster_flag] = new_value
                    else:
                        new_cluster.neighbor_flow[neighbor_flag] = 0
                        neighbor_cluster.neighbor_flow[cluster_flag] = 0
                else:
                    new_cluster.neighbor_flow[neighbor_flag] = 0
                    neighbor_cluster.neighbor_flow[cluster_flag] = 0
            elif neighbor_flag == 0:
                k_value = len(new_cluster.member_flow) + 1
                new_Zx = (sum(new_cluster.member_value1) + attribution_value[nei][0]) / k_value
                new_Zy = (sum(new_cluster.member_value2) + attribution_value[nei][1]) / k_value
                new_sign = distinguish_sign(new_Zx, new_Zy)

                if new_sign == G_sign:
                    new_value = sqrt((type_num - 1) / (type_num / k_value - 1)) * new_Zx * new_Zy

                    if abs(new_value) > abs(new_cluster.stat):
                        new_cluster.neighbor_flow[nei] = new_value
                    else:
                        new_cluster.neighbor_flow[nei] = 0
                else:
                    new_cluster.neighbor_flow[nei] = 0

    neighbor_flow.pop(index1)

    if len(new_cluster.neighbor_flow) != 0:
        new_cluster.neighbor_flow = dict(sorted(new_cluster.neighbor_flow.items(),
                                                key=lambda x: abs(x[1]), reverse=True))

        new_cluster.next_stat = list(new_cluster.neighbor_flow.values())[0]
    else:
        new_cluster.next_stat = 0


def distinguish_sign(value_1, value_2):
    if value_1 > 0 and value_2 > 0:
        new_sign = 1

    elif value_2 < 0 < value_1:
        new_sign = 2

    elif value_1 < 0 < value_2:
        new_sign = 3

    else:
        new_sign = 4

    return new_sign


def add_cluster_neighbor(add_flag, cluster_flag, new_cluster,
                         attribution_value, all_cluster, G_sign, type_num, change_cluster):
    for l in all_cluster[add_flag - 1].neighbor_flow.keys():
        if type(l) != str:
            neighbor_cluster = all_cluster[l - 1]
            flag_key = '{}-{}'.format(neighbor_cluster.member_flow[0][0], neighbor_cluster.member_flow[0][1])
            if flag_key not in new_cluster.neigh_flow_flag.keys():
                k_value = len(new_cluster.member_flow) + len(neighbor_cluster.member_flow)
                new_Zx = (sum(new_cluster.member_value1) + sum(neighbor_cluster.member_value1)) / k_value
                new_Zy = (sum(new_cluster.member_value2) + sum(neighbor_cluster.member_value2)) / k_value
                new_sign = distinguish_sign(new_Zx, new_Zy)
                for e in new_cluster.member_flow:
                    flag_key = '{}-{}'.format(e[0], e[1])
                    neighbor_cluster.neigh_flow_flag[flag_key] = 1
                for e in neighbor_cluster.member_flow:
                    flag_key = '{}-{}'.format(e[0], e[1])
                    new_cluster.neigh_flow_flag[flag_key] = 1
                change_cluster[l] = 0

                if new_sign == G_sign:
                    new_value = sqrt((type_num - 1) / (type_num / k_value - 1)) * new_Zx * new_Zy
                    old_value = sqrt((type_num - 1) / (type_num / len(neighbor_cluster.member_value1) - 1)) \
                                * np.mean(neighbor_cluster.member_value1) * \
                                np.mean(neighbor_cluster.member_value2)

                    if abs(new_value) > abs(new_cluster.stat) and abs(new_value) > abs(old_value):
                        new_cluster.neighbor_flow[l] = new_value
                        neighbor_cluster.neighbor_flow[cluster_flag] = new_value
                    else:
                        new_cluster.neighbor_flow[l] = 0
                        neighbor_cluster.neighbor_flow[cluster_flag] = 0
                else:
                    new_cluster.neighbor_flow[l] = 0
                    neighbor_cluster.neighbor_flow[cluster_flag] = 0
        else:
            new_flow = l.split('-')
            if l not in new_cluster.neigh_flow_flag.keys():
                k_value = len(new_cluster.member_flow) + 1
                new_Zx = (sum(new_cluster.member_value1) + attribution_value[l][0]) / k_value
                new_Zy = (sum(new_cluster.member_value2) + attribution_value[l][1]) / k_value
                new_sign = distinguish_sign(new_Zx, new_Zy)

                new_cluster.neigh_flow_flag[l] = 1

                if new_sign == G_sign:
                    new_value = sqrt((type_num - 1) / (type_num / k_value - 1)) * new_Zx * new_Zy

                    if abs(new_value) > abs(new_cluster.stat):
                        new_cluster.neighbor_flow[l] = new_value
                    else:
                        new_cluster.neighbor_flow[l] = 0
                else:
                    new_cluster.neighbor_flow[l] = 0


def add_flow_neighbor(add_flag, cluster_flag, new_cluster, neighbor_flow, index_matrix,
                        attribution_value, all_cluster, G_sign, type_num, change_cluster):
    for nei in neighbor_flow[add_flag]:
        new_flow = nei.split('-')
        if nei not in new_cluster.neigh_flow_flag.keys():
            new_cluster.neigh_flow_flag[nei] = 1
            neighbor_flag = index_matrix[int(new_flow[0]) - 1, int(new_flow[1]) - 1]
            if neighbor_flag != 0 and all_cluster[neighbor_flag - 1].sign == 1:
                neighbor_cluster = all_cluster[neighbor_flag - 1]
                k_value = len(new_cluster.member_flow) + len(neighbor_cluster.member_flow)
                new_Zx = (sum(new_cluster.member_value1) + sum(neighbor_cluster.member_value1)) / k_value
                new_Zy = (sum(new_cluster.member_value2) + sum(neighbor_cluster.member_value2)) / k_value
                new_sign = distinguish_sign(new_Zx, new_Zy)
                for e in new_cluster.member_flow:
                    flag_key = '{}-{}'.format(e[0], e[1])
                    neighbor_cluster.neigh_flow_flag[flag_key] = 1
                for e in neighbor_cluster.member_flow:
                    flag_key = '{}-{}'.format(e[0], e[1])
                    new_cluster.neigh_flow_flag[flag_key] = 1
                change_cluster[neighbor_flag] = 0

                if new_sign == G_sign:
                    new_value = sqrt((type_num - 1) / (type_num / k_value - 1)) * new_Zx * new_Zy
                    old_value = sqrt((type_num - 1) / (type_num / len(neighbor_cluster.member_value1) - 1)) \
                                * np.mean(neighbor_cluster.member_value1) * \
                                np.mean(neighbor_cluster.member_value2)

                    if abs(new_value) > abs(new_cluster.stat) and abs(new_value) > abs(old_value):
                        new_cluster.neighbor_flow[neighbor_flag] = new_value
                        neighbor_cluster.neighbor_flow[cluster_flag] = new_value
                    else:
                        new_cluster.neighbor_flow[neighbor_flag] = 0
                        neighbor_cluster.neighbor_flow[cluster_flag] = 0
                else:
                    new_cluster.neighbor_flow[neighbor_flag] = 0
                    neighbor_cluster.neighbor_flow[cluster_flag] = 0
                if add_flag in neighbor_cluster.neighbor_flow.keys():
                    neighbor_cluster.neighbor_flow.pop(add_flag)
            elif neighbor_flag == 0:
                k_value = len(new_cluster.member_flow) + 1
                new_Zx = (sum(new_cluster.member_value1) + attribution_value[nei][0]) / k_value
                new_Zy = (sum(new_cluster.member_value2) + attribution_value[nei][1]) / k_value
                new_sign = distinguish_sign(new_Zx, new_Zy)

                if new_sign == G_sign:
                    new_value = sqrt((type_num - 1) / (type_num / k_value - 1)) * new_Zx * new_Zy

                    if abs(new_value) > abs(new_cluster.stat):
                        new_cluster.neighbor_flow[nei] = new_value
                    else:
                        new_cluster.neighbor_flow[nei] = 0
                else:
                    new_cluster.neighbor_flow[nei] = 0


def update_old_neighbor(new_cluster, cluster_flag, add_flag, type_num,
                        G_sign, all_cluster, change_cluster, attribution_value):
    for l in new_cluster.neighbor_flow.keys():
        if type(l) != str:
            neighbor_cluster = all_cluster[l - 1]
            if neighbor_cluster.sign == 1:
                k_value = len(new_cluster.member_flow) + len(neighbor_cluster.member_flow)
                new_Zx = (sum(new_cluster.member_value1) + sum(neighbor_cluster.member_value1)) / k_value
                new_Zy = (sum(new_cluster.member_value2) + sum(neighbor_cluster.member_value2)) / k_value
                new_sign = distinguish_sign(new_Zx, new_Zy)
                if type(add_flag) != str:
                    if add_flag not in neighbor_cluster.neighbor_flow.keys():
                        for e in all_cluster[add_flag - 1].member_flow:
                            flag_key = '{}-{}'.format(e[0], e[1])
                            neighbor_cluster.neigh_flow_flag[flag_key] = 1
                else:
                    if add_flag not in neighbor_cluster.neighbor_flow.keys():
                        # new_flow = add_flag.split('-')
                        neighbor_cluster.neigh_flow_flag[add_flag] = 1

                if new_sign == G_sign:
                    new_value = sqrt((type_num - 1) / (type_num / k_value - 1)) * new_Zx * new_Zy
                    old_value = sqrt((type_num - 1) / (type_num / len(neighbor_cluster.member_value1) - 1)) \
                                * np.mean(neighbor_cluster.member_value1) * \
                                np.mean(neighbor_cluster.member_value2)

                    if abs(new_value) > abs(new_cluster.stat) and abs(new_value) > abs(old_value):
                        new_cluster.neighbor_flow[l] = new_value
                        neighbor_cluster.neighbor_flow[cluster_flag] = new_value
                    else:
                        new_cluster.neighbor_flow[l] = 0
                        neighbor_cluster.neighbor_flow[cluster_flag] = 0
                else:
                    new_cluster.neighbor_flow[l] = 0
                    neighbor_cluster.neighbor_flow[cluster_flag] = 0

                change_cluster[l] = 0

        else:
            k_value = len(new_cluster.member_flow) + 1
            new_Zx = (sum(new_cluster.member_value1) + attribution_value[l][0]) / k_value
            new_Zy = (sum(new_cluster.member_value2) + attribution_value[l][1]) / k_value
            new_sign = distinguish_sign(new_Zx, new_Zy)

            if new_sign == G_sign:
                new_value = sqrt((type_num - 1) / (type_num / k_value - 1)) * new_Zx * new_Zy

                if abs(new_value) > abs(new_cluster.stat):
                    new_cluster.neighbor_flow[l] = new_value
                else:
                    new_cluster.neighbor_flow[l] = 0
            else:
                new_cluster.neighbor_flow[l] = 0


# 冷热点探测
def get_desired_cluster(seed_flow, region_num, attribute_value, type_num, G_sign, neighbor_flow):
    current_cluster = list()
    index_matrix = np.zeros((region_num, region_num), dtype=np.int)

    current_index = 0
    for l in seed_flow:
        new_cluster = Hot_Spot_region(l[0], l[1])

        index_matrix[l[0] - 1, l[1] - 1] = current_index + 1
        current_index += 1

        new_cluster.member_flow.append(l)

        new_Zx = attribute_value['{}-{}'.format(l[0], l[1])][0]
        new_Zy = attribute_value['{}-{}'.format(l[0], l[1])][1]
        new_cluster.member_value1.append(new_Zx)
        new_cluster.member_value2.append(new_Zy)
        new_cluster.stat = new_Zx * new_Zy

        flag_key = '{}-{}'.format(l[0], l[1])
        new_cluster.neigh_flow_flag[flag_key] = 1

        current_cluster.append(new_cluster)
    # print(time.process_time())

    sorted_dict = dict()
    sort_index = 1

    for i in range(len(current_cluster)):
        initial_cluster_neighbor(i + 1, attribute_value, index_matrix, type_num, G_sign, current_cluster, neighbor_flow)
        sorted_dict[sort_index] = abs(current_cluster[i].next_stat)
        sort_index += 1

    sorted_list = list(sorted(sorted_dict.items(), key=lambda x: x[1], reverse=True))
    # print(time.process_time())
    finish_cluster = 0
    change_cluster = dict()

    while sorted_list[0][1] != 0:
        sign = 0
        # print(finish_cluster)
        for p in sorted_list:
            l = current_cluster[p[0] - 1]
            if l.sign == 1:
                for fir_neigh in l.neighbor_flow.keys():
                    if type(fir_neigh) != str:
                        neighbor_cluster = current_cluster[fir_neigh - 1]

                        l.member_value1 = l.member_value1 + neighbor_cluster.member_value1
                        l.member_value2 = l.member_value2 + neighbor_cluster.member_value2
                        l.stat = l.next_stat

                        neighbor_cluster.sign = 0
                        neighbor_cluster.next_stat = 0
                        sorted_dict.pop(fir_neigh)
                        finish_cluster += 1
                        if fir_neigh in change_cluster:
                            change_cluster.pop(fir_neigh)

                        for j in neighbor_cluster.member_flow:
                            index_matrix[j[0] - 1, j[1] - 1] = \
                                index_matrix[l.member_flow[0][0] - 1, l.member_flow[0][1] - 1]
                            l.member_flow.append(j)

                        update_old_neighbor(l, p[0], fir_neigh, type_num, G_sign,
                                            current_cluster, change_cluster, attribute_value)
                        add_cluster_neighbor(fir_neigh, p[0], l, attribute_value, current_cluster,
                                             G_sign, type_num, change_cluster)
                        for e in neighbor_cluster.neighbor_flow.keys():
                            if type(e) != str:
                                current_cluster[e - 1].neighbor_flow.pop(fir_neigh)

                        if len(l.neighbor_flow) != 0:
                            l.neighbor_flow = dict(sorted(l.neighbor_flow.items(),
                                                          key=lambda x: abs(x[1]), reverse=True))
                            l.next_stat = list(l.neighbor_flow.values())[0]
                        else:
                            l.next_stat = 0
                        sorted_dict[p[0]] = abs(l.next_stat)

                    else:
                        l.member_value1.append(attribute_value[fir_neigh][0])
                        l.member_value2.append(attribute_value[fir_neigh][1])
                        l.stat = l.next_stat
                        l.neighbor_flow.pop(fir_neigh)

                        new_flow = fir_neigh.split('-')
                        index_matrix[int(new_flow[0]) - 1, int(new_flow[1]) - 1] = \
                                index_matrix[l.member_flow[0][0] - 1, l.member_flow[0][1] - 1]
                        l.member_flow.append([int(new_flow[0]), int(new_flow[1])])

                        update_old_neighbor(l, p[0], fir_neigh, type_num, G_sign,
                                            current_cluster, change_cluster, attribute_value)
                        add_flow_neighbor(fir_neigh, p[0], l, neighbor_flow, index_matrix,
                                          attribute_value, current_cluster, G_sign, type_num, change_cluster)
                        # for e in neighbor_flow[fir_neigh]:
                        #     new_flow = e.split('-')
                        #     neighbor_flag = index_matrix[int(new_flow[0]) - 1, int(new_flow[1]) - 1]
                        #     if neighbor_flag != 0:
                        #         current_cluster[neighbor_flag - 1].neighbor_flow.pop(fir_neigh)
                        # neighbor_flow.pop(fir_neigh)

                        if len(l.neighbor_flow) != 0:
                            l.neighbor_flow = dict(sorted(l.neighbor_flow.items(), key=lambda x: abs(x[1]),
                                                          reverse=True))
                            l.next_stat = list(l.neighbor_flow.values())[0]
                        else:
                            l.next_stat = 0
                        sorted_dict[p[0]] = abs(l.next_stat)
                    break
                if l.next_stat != 0:
                    oo = 0
                else:
                    for change_index in change_cluster:
                        current_cluster[change_index - 1].neighbor_flow = dict(sorted(
                            current_cluster[change_index - 1].neighbor_flow.items(), key=lambda x: abs(x[1]),
                            reverse=True))
                        if len(current_cluster[change_index - 1].neighbor_flow) != 0:
                            current_cluster[change_index - 1].next_stat = list(
                                current_cluster[change_index - 1].neighbor_flow.values())[0]
                        else:
                            current_cluster[change_index - 1].next_stat = 0

                        sorted_dict[change_index] = abs(current_cluster[change_index - 1].next_stat)
                    sorted_list = list(sorted(sorted_dict.items(), key=lambda x: x[1], reverse=True))

                sign = 1
                # print(time.process_time())
                break

        if sign == 0:
            break

    new_current_cluster = sorted(current_cluster, key=lambda x: abs(x.stat), reverse=True)

    return new_current_cluster


# 删除重叠的簇
def defined_region(search_region):
    interest_region = list()

    for l in search_region:
        if l.sign == 1:
            for j in l.member_flow:
                if j[0] not in l.member_o:
                    l.member_o.append(j[0])
                if j[1] not in l.member_d:
                    l.member_d.append(j[1])
            interest_region.append(l)
    return interest_region


def delete_overlap_region(interest_region):
    region_num = len(interest_region)
    sign = [0 for i in range(region_num)]

    for i in range(region_num - 1):
        for j in range(i + 1, region_num):
            if sign[j] == 0:
                if len(list(set(interest_region[i].member_o).intersection(set(interest_region[j].member_o)))) != 0 and \
                        len(list(set(interest_region[i].member_d).intersection(set(interest_region[j].member_d)))) != 0:
                    sign[j] = 1
            else:
                continue

    final_region = list()
    for i in range(region_num):
        if sign[i] == 0:
            final_region.append(interest_region[i])
        else:
            continue

    return final_region


# 假设检验
def Monte_carlo_simulation(seed_flow, seed_flow_value, final_region, rep_num, sig_level, type_num, G_sign):
    sort_num = [i for i in range(len(seed_flow))]
    keys = [l for l in seed_flow_value.keys()]
    cluster_num = len(final_region)
    p_value = [0 for i in range(cluster_num)]
    final_true_region = list()

    for i in range(rep_num):
        # print(sort_num)
        random.shuffle(sort_num)
        for j in range(cluster_num):
            if p_value[j] <= 50:
                new_value1 = list()
                new_value2 = list()
                for e in final_region[j].member_flow:
                    location = seed_flow.index(e)
                    # print(sort_num[location])
                    new_value1.append(seed_flow_value[keys[sort_num[location]]][0])
                    new_value2.append(seed_flow_value[keys[sort_num[location]]][1])
                average_1 = np.mean(new_value1)
                average_2 = np.mean(new_value2)
                new_sign = distinguish_sign(average_1, average_2)
                if new_sign == G_sign:
                    new_stat = sqrt((type_num - 1) / (type_num / len(final_region[j].member_flow) - 1)) \
                               * average_1 * average_2
                    if abs(new_stat) >= abs(final_region[j].stat):
                        p_value[j] += 1
        # print('simulation')
        # print(time.process_time())

    for i in range(cluster_num):
        final_region[i].p_value = (p_value[i]) / (rep_num + 1)
    final_region = sorted(final_region, key=lambda x: x.p_value)

    # for i in range(cluster_num, 0, -1):
    #     if final_region[i - 1].p_value <= i * sig_level / cluster_num:
    #         final_true_region = final_region[:i]
    #     else:
    #         final_true_region = final_region[:1]
    #         break

    for i in range(0, cluster_num, 1):
        if final_region[i].p_value > sig_level:
            final_true_region = final_region[:i]
            break

    return final_true_region


# 存储显著类
def save_cluster(final_true_region, cf):
    file = open(cf, 'w')

    for l in final_true_region:
        neighbor_o = ['{}'.format(j) for j in l.member_o]
        s_o = ','.join(neighbor_o)
        file.writelines(s_o + ';')

        neighbor_d = ['{}'.format(j) for j in l.member_d]
        s_d = ','.join(neighbor_d)
        file.writelines(s_d + ';')

        file.writelines('{};'.format(l.stat))

        neighbor_flow = ['{} {}'.format(j[0], j[1]) for j in l.member_flow]
        s_flow = ','.join(neighbor_flow)
        file.writelines(s_flow + ';')

        file.writelines('{}\n'.format(l.p_value))

    file.close()


# 初始化拓展区域
def initial_start_poly(seed_flow_value, chosen_num):
    poly = list()
    weight = list()
    num = 0

    for l in seed_flow_value.keys():
        s1 = seed_flow_value[l][0]
        s2 = seed_flow_value[l][1]
        weight.append(abs(s1 * s2))
        poly.append(num)
        num += 1

    s = np.sum(weight)
    weight = np.array(weight) / s

    start_poly = np.random.choice(poly, chosen_num, p=weight)

    return start_poly


def run_simulation():
    # simulation_data:
    rep_num = 999
    sig_level = 0.01
    print("Result for data sets contain a single cluster")
    cluster_name = ['H-H Cluster:', 'H-L Cluster:', 'L-H Cluster:', 'L-L Cluster:']
    for i in range(4):
        data_dir = 'C:/Users/12449/Desktop/BiflowAMOEBA/IJGIS-code&data(2021)/Synthetic data/D{}'.format(i + 1)
        G_sign = i + 1
        nf = path.join(data_dir, 'neighbors.txt')
        cf = path.join(data_dir, 'Cluster_flow.txt')
        ff = path.join(data_dir, 'flow.txt')

        cluster_flow = read_cluster_flow(cf)
        neighbors = read_neighbor(nf)
        attribute_values = dict()
        read_simulation_flow(ff, attribute_values)

        unique_modes = normal_flow_attribute_value(attribute_values)
        seed_flows, seed_flow_values = find_seed_flows(G_sign, unique_modes, attribute_values)

        region_flag = np.zeros((len(neighbors), len(neighbors)), dtype=np.int8)
        cal_flow_flag(seed_flows, region_flag, G_sign, attribute_values)
        flow_neighbor = initial_flow_neighbor(seed_flows, neighbors, region_flag)
        cluster_regions = get_desired_cluster(seed_flows, len(neighbors), attribute_values,
                                              len(unique_modes), G_sign, flow_neighbor)

        final_regions = defined_region(cluster_regions)
        #
        final_true_regions = Monte_carlo_simulation(seed_flows, seed_flow_values,
                                                    final_regions, rep_num, sig_level, len(unique_modes), G_sign)

        final_true_regions = sorted(final_true_regions, key=lambda x: abs(x.stat), reverse=True)
        final_true_regions = delete_overlap_region(final_true_regions)
        #
        TP = 0
        final_num = 0
        for l in final_true_regions:
            final_num += len(l.member_flow)
            for n in l.member_flow:
                flow_key = '{}-{}'.format(n[0], n[1])
                if flow_key in cluster_flow[0]:
                    TP += 1
        if TP != 0:
            Precision = TP / final_num
            Recall = TP / len(cluster_flow[0])
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
    neighbors = read_neighbor(nf)
    attribute_values = dict()
    read_simulation_flow(ff, attribute_values)

    unique_modes = normal_flow_attribute_value(attribute_values)
    for i in range(2, 4):
        G_sign = i + 1
        seed_flows, seed_flow_values = find_seed_flows(G_sign, unique_modes, attribute_values)

        region_flag = np.zeros((len(neighbors), len(neighbors)), dtype=np.int8)
        cal_flow_flag(seed_flows, region_flag, G_sign, attribute_values)
        flow_neighbor = initial_flow_neighbor(seed_flows, neighbors, region_flag)
        cluster_regions = get_desired_cluster(seed_flows, len(neighbors), attribute_values,
                                              len(unique_modes), G_sign, flow_neighbor)

        final_regions = defined_region(cluster_regions)
        final_true_regions = Monte_carlo_simulation(seed_flows, seed_flow_values,
                                                    final_regions, rep_num, sig_level, len(unique_modes), G_sign)

        final_true_regions = sorted(final_true_regions, key=lambda x: abs(x.stat), reverse=True)
        final_true_regions = delete_overlap_region(final_true_regions)
        TP = 0
        final_num = 0
        for l in final_true_regions:
            final_num += len(l.member_flow)
            for n in l.member_flow:
                flow_key = '{}-{}'.format(n[0], n[1])
                if flow_key in cluster_flow[G_sign - 1]:
                    TP += 1
        if TP != 0:
            Precision = TP / final_num
            Recall = TP / len(cluster_flow[G_sign - 1])
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
    run_simulation()
