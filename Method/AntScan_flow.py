# -*- coding=utf-8 -*-

import numpy as np
from math import log
import random
from os import path


# 冷热点类
class Cluster_region(object):
    def __init__(self, c_id):
        self.id = c_id
        self.member_flow = list()
        self.member_o = list()
        self.member_d = list()
        self.alpha = 0
        self.p_value = 0
        self.p1 = 0
        self.p2 = 0
        self.q1 = 0
        self.q2 = 0
        self.Na = 0
        self.Nb = 0


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

    file.close()

    return spatial_neighbors


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


# def read_cluster_flow(cluster_region):
#     cluster_flow = [list() for i in range(int(len(cluster_region)/2))]
#     for i in range(int(len(cluster_region)/2)):
#         for l in cluster_region[2 * i]:
#             for n in cluster_region[2 * i + 1]:
#                 cluster_flow[i].append('{}-{}'.format(l, n))
#     return cluster_flow


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


def read_flows(ff, flows, sign):
    file = open(ff, 'r')
    sum_num = 0

    if sign == 1:
        while True:
            str_line = file.readline()
            if not str_line:
                break

            str_lines = str_line.split(',')
            if int(str_lines[0]) != 0 and int(str_lines[1]) != 0:
                sum_num += 1
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
            if int(str_lines[0]) != 0 and int(str_lines[1]) != 0:
                sum_num += 1
                indices = '{}-{}'.format(int(str_lines[0]), int(str_lines[1]))
                if indices not in flows.keys():
                    flows[indices] = [0, 1]
                else:
                    flows[indices][1] += 1

    file.close()

    return sum_num


def read_simulation_flow(ff, flows):
    file = open(ff, 'r')
    num_A = 0
    num_B = 0

    while True:
        str_line = file.readline()
        if not str_line:
            break
        str_lines = str_line.split(',')

        indices = '{}-{}'.format(int(str_lines[0]), int(str_lines[1]))
        flows[indices] = [float(str_lines[2]), float(str_lines[3])]
        num_A += float(str_lines[2])
        num_B += float(str_lines[3])

    file.close()

    return num_A, num_B


# 标识流
def cal_flow_flag(attribute_value, region_flag):
    for l in attribute_value.keys():
        location = l.split('-')
        region_flag[int(location[0]) - 1, int(location[1]) - 1] = 1


# 初始化拓展区域
def initial_start_poly(attribute_values, Tau_origin, Tau_destination, ant_num):
    poly = list()
    weight = list()
    polygons = list()
    num = 0

    for l in attribute_values.keys():
        location = l.split('-')
        s1 = np.sum(Tau_origin[int(location[0]) - 1, :])
        s2 = np.sum(Tau_destination[int(location[1]) - 1, :])
        weight.append(s1 + s2)
        poly.append(num)
        polygons.append([int(location[0]), int(location[1])])
        num += 1

    s = np.sum(weight)
    weight = np.array(weight) / s

    start_poly = np.random.choice(poly, ant_num, p=weight)

    return start_poly, polygons


# 初始化信息素
def initial_pheromone(Tau_origin, Tau_destination, attribute_values, neighbor):
    for l in attribute_values.keys():
        location = l.split('-')
        for j in neighbor[int(location[0]) - 1]:
            Tau_origin[j - 1, int(location[0]) - 1] += (attribute_values[l][0] + attribute_values[l][1])
            Tau_origin[int(location[0]) - 1, j - 1] += (attribute_values[l][0] + attribute_values[l][1])
        for j in neighbor[int(location[1]) - 1]:
            Tau_destination[j - 1, int(location[1]) - 1] += (attribute_values[l][0] + attribute_values[l][1])
            Tau_destination[int(location[1]) - 1, j - 1] += (attribute_values[l][0] + attribute_values[l][1])


# 寻找流邻居
def find_flow_neighbor(origin, destination, neighbor_o, neighbor_d, add_o, add_d, flow_flag,
                       Tau_origin, Tau_destination, neighbor_flow, Tau_flow):
    for l in add_o:
        for j in neighbor_d:
            if flow_flag[l - 1, j - 1] == 1:
                index1 = '{}-{}'.format(l, j)
                neighbor_flow[index1] = len(Tau_flow)
                Tau_flow.append(Tau_origin[origin, l - 1] + Tau_destination[destination, j - 1])

    for l in add_d:
        for j in neighbor_o:
            if flow_flag[j - 1, l - 1] == 1:
                index1 = '{}-{}'.format(j, l)
                neighbor_flow[index1] = len(Tau_flow)
                Tau_flow.append(Tau_origin[origin, j - 1] + Tau_destination[destination, l - 1])

    if len(neighbor_flow) != 0:
        weight = np.array([Tau_flow[neighbor_flow[l]] for l in neighbor_flow.keys()])
        s = weight.sum()
        if s == 0:
            weight = np.ones(len(weight)) / len(weight)
        else:
            weight = weight / s

        next_indices = [i for i in range(len(weight))]
        next_flows = [l for l in neighbor_flow.keys()]

        next_index = np.random.choice(next_indices, 1, p=weight)[0]

        return next_flows[next_index]

    else:
        return []


def cluster_alpha_calculate(ant_cluster, attribute_value, num_A, num_B, region_num, region_flag):
    cluster_region = list()

    for i in range(len(ant_cluster)):
        origin_flag = np.zeros(region_num)
        destination_flag = np.zeros(region_num)
        new_region = Cluster_region(i)
        local_num_A = 0
        local_num_B = 0

        for j in ant_cluster[i]:
            if origin_flag[j[0] - 1] == 0:
                origin_flag[j[0] - 1] = 1
                new_region.member_o.append(j[0])
            if destination_flag[j[1] - 1] == 0:
                destination_flag[j[1] - 1] = 1
                new_region.member_d.append(j[1])

            new_region.member_flow.append(j)

            index1 = '{}-{}'.format(j[0], j[1])
            local_num_A += attribute_value[index1][0]
            local_num_B += attribute_value[index1][1]
        new_region.p1 = local_num_A / (local_num_A + local_num_B)
        new_region.q1 = (num_A - local_num_A) / (num_A + num_B - local_num_A - local_num_B)
        new_region.p2 = local_num_B / (local_num_A + local_num_B)
        new_region.q2 = (num_B - local_num_B) / (num_A + num_B - local_num_A - local_num_B)
        new_region.Na = local_num_A
        new_region.Nb = local_num_B
        new_alpha = 0
        if new_region.p1 != 0:
            new_alpha += local_num_A * log(new_region.p1)
        if new_region.p2 != 0:
            new_alpha += local_num_B * log(new_region.p2)
        if new_region.q1 != 0:
            new_alpha += (num_A - local_num_A) * log(new_region.q1)
        if new_region.q2 != 0:
            new_alpha += (num_B - local_num_B) * log(new_region.q2)
        new_alpha -= (num_A * log(num_A / (num_A + num_B)) + num_B * log(num_B / (num_A + num_B)))
        new_region.alpha = new_alpha
        cluster_region.append(new_region)
    cluster_region = sorted(cluster_region, key=lambda x: abs(x.alpha), reverse=True)

    return cluster_region


# # 蚁群算法将流视为整体
def Ant_amoeba(neighbor, Ant_num, Iteration_time, EliteAnt_num, EvapCoeff, attribute_value, region_flag, num_A, num_B):
    poly_num = len(neighbor)

    Tau_origin = np.zeros((poly_num, poly_num))
    Tau_destination = np.zeros((poly_num, poly_num))
    initial_pheromone(Tau_origin, Tau_destination, attribute_value, neighbor)

    Best_cluster = list()
    cluster_size = np.random.rand(20) * poly_num
    cluster_num = np.mean(cluster_size)
    size_std = np.std(cluster_size) / 2

    for i in range(Iteration_time):
        ant_cluster = list()
        start_polygons, polygons = initial_start_poly(attribute_value, Tau_origin, Tau_destination, Ant_num)

        for j in range(Ant_num):
            start_poly = start_polygons[j]
            neighbor_flow = dict()
            cluster_region = list()
            neigh_flag_o = [0 for n in range(len(neighbor))]
            neigh_flag_d = [0 for n in range(len(neighbor))]
            Tau_flow = list()

            C_length = np.random.normal(loc=cluster_num, scale=size_std)

            new_add_o = list()
            neighbor_o = list()
            new_add_o.append(polygons[start_poly][0])
            neigh_flag_o[polygons[start_poly][0] - 1] = 1
            for l in neighbor[polygons[start_poly][0] - 1]:
                new_add_o.append(l)
                neighbor_o.append(l)
                neigh_flag_o[l - 1] = 1

            new_add_d = list()
            neighbor_d = list()
            new_add_d.append(polygons[start_poly][1])
            neigh_flag_d[polygons[start_poly][1] - 1] = 1
            for l in neighbor[polygons[start_poly][1] - 1]:
                new_add_d.append(l)
                neighbor_d.append(l)
                neigh_flag_d[l - 1] = 1

            l = find_flow_neighbor(polygons[start_poly][0] - 1, polygons[start_poly][1] - 1,
                                   neighbor_o, neighbor_d, new_add_o, new_add_d, region_flag,
                                   Tau_origin, Tau_destination, neighbor_flow, Tau_flow)

            cluster_region.append([polygons[start_poly][0], polygons[start_poly][1]])

            while C_length > 0:
                if len(l) == 0:
                    break
                next_flow = l.split('-')
                next_o = int(next_flow[0])
                next_d = int(next_flow[1])
                cluster_region.append([next_o, next_d])

                new_add_o = list()
                for e in neighbor[next_o - 1]:
                    if neigh_flag_o[e - 1] == 0:
                        neigh_flag_o[e - 1] = 1
                        new_add_o.append(e)
                        neighbor_o.append(e)

                new_add_d = list()
                for e in neighbor[next_d - 1]:
                    if neigh_flag_d[e - 1] == 0:
                        neigh_flag_d[e - 1] = 1
                        new_add_d.append(e)
                        neighbor_d.append(e)

                neighbor_flow.pop(l)

                l = find_flow_neighbor(next_o - 1, next_d - 1, neighbor_o, neighbor_d, new_add_o, new_add_d, region_flag
                                       , Tau_origin, Tau_destination, neighbor_flow, Tau_flow)

                C_length -= 1

            ant_cluster.append(cluster_region)
        cluster_regions = cluster_alpha_calculate(ant_cluster, attribute_value, num_A, num_B, len(neighbor), region_flag)

        for n in range(EliteAnt_num):
            Best_cluster.append(cluster_regions[n])
            for m in cluster_regions[n].member_o:
                Tau_origin[m - 1, np.array(cluster_regions[n].member_o) - 1] = (1 - EvapCoeff) * \
                                                                               Tau_origin[m - 1, np.array(
                                                                                   cluster_regions[n].member_o) - 1] \
                                                                               + (EliteAnt_num - n) ** 2
                Tau_origin[np.array(cluster_regions[n].member_o) - 1, m - 1] = (1 - EvapCoeff) * \
                                                                               Tau_origin[np.array(cluster_regions[
                                                                                                       n].member_o) - 1, m - 1] \
                                                                               + (EliteAnt_num - n) ** 2
            for m in cluster_regions[n].member_d:
                Tau_destination[m - 1, np.array(cluster_regions[n].member_d) - 1] = (1 - EvapCoeff) * \
                                                                                    Tau_destination[m - 1, np.array(
                                                                                        cluster_regions[
                                                                                            n].member_d) - 1] \
                                                                                    + (EliteAnt_num - n) ** 2
                Tau_destination[np.array(cluster_regions[n].member_d) - 1, m - 1] = (1 - EvapCoeff) * \
                                                                                    Tau_destination[np.array(
                                                                                        cluster_regions[
                                                                                            n].member_d) - 1, m - 1] \
                                                                                    + (EliteAnt_num - n) ** 2

    Best_cluster = sorted(Best_cluster, key=lambda x: abs(x.alpha), reverse=True)

    return Best_cluster


# 删除重叠的簇
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
def Monte_carlo_simulation(attribute_values, final_region, rep_num, sig_level, num_A, num_B):
    keys = [l for l in attribute_values.keys()]
    sort_num = [i for i in range(len(keys))]
    cluster_num = len(final_region)
    p_value = [0 for i in range(cluster_num)]
    final_true_region = list()

    for i in range(rep_num):
        random.shuffle(sort_num)
        for j in range(cluster_num):
            if p_value[j] <= sig_level * (rep_num + 1):
                local_num_A = 0
                local_num_B = 0
                for e in final_region[j].member_flow:
                    index1 = '{}-{}'.format(e[0], e[1])
                    location = keys.index(index1)
                    # print(sort_num[location])
                    local_num_A += attribute_values[keys[sort_num[location]]][0]
                    local_num_B += attribute_values[keys[sort_num[location]]][1]
                new_alpha = 0
                if local_num_A != 0:
                    new_alpha += local_num_A * log(local_num_A / (local_num_A + local_num_B))
                if local_num_B != 0:
                    new_alpha += local_num_B * log(local_num_B / (local_num_A + local_num_B))
                if (num_A - local_num_A) != 0:
                    new_alpha += (num_A - local_num_A) * log((num_A - local_num_A) / (num_A + num_B - local_num_A - local_num_B))
                if (num_B - local_num_B) != 0:
                    new_alpha += (num_B - local_num_B) * log((num_B - local_num_B) / (num_A + num_B - local_num_A - local_num_B))
                new_alpha -= (num_A * log(num_A / (num_A + num_B)) + num_B * log(num_B / (num_A + num_B)))
                if new_alpha >= final_region[j].alpha:
                    p_value[j] += 1

    for i in range(cluster_num):
        final_region[i].p_value = (p_value[i] + 1) / (rep_num + 1)
    final_region = sorted(final_region, key=lambda x: x.p_value)

    for i in range(cluster_num, 0, -1):
        if final_region[i - 1].p_value <= sig_level:
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

        file.writelines('{},{},{},{},{},{},{};'.format(l.alpha, l.p1, l.q1, l.p2, l.q2, l.Na, l.Nb))

        neighbor_flow = ['{} {}'.format(j[0], j[1]) for j in l.member_flow]
        s_flow = ','.join(neighbor_flow)
        file.writelines(s_flow + ';')

        file.writelines('{}\n'.format(l.p_value))

    file.close()


def run():
    rep_num = 999
    sig_level = 0.01
    print("Result for data sets contain a single cluster")
    cluster_name = ['H-H Cluster:', 'H-L Cluster:', 'L-H Cluster:', 'L-L Cluster:']
    ant_num = 200
    iteration_time = 100
    EliteAntNum = 30
    EvapCoeff = 0.1

    for i in range(4):
        data_dir = 'C:/Users/12449/Desktop/BiflowAMOEBA/IJGIS-code&data(2021)/Synthetic data/D{}'.format(i + 1)
        nf = path.join(data_dir, 'neighbors.txt')
        cf = path.join(data_dir, 'Cluster_flow.txt')
        ff = path.join(data_dir, 'flow.txt')

        cluster_flow = read_cluster_flow(cf)
        neighbors = read_neighbor(nf)
        attribute_values = dict()
        num_A, num_B = read_simulation_flow(ff, attribute_values)

        region_flag = np.zeros((len(neighbors), len(neighbors)))
        cal_flow_flag(attribute_values, region_flag)
        new_region = Ant_amoeba(neighbors, ant_num, iteration_time, EliteAntNum, EvapCoeff,
                                attribute_values, region_flag, num_A, num_B)

        final_regions = delete_overlap_region(new_region)
        final_true_regions = Monte_carlo_simulation(attribute_values, final_regions, rep_num, sig_level, num_A, num_B)

        Sign = 1
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
            Recall = TP / len(cluster_flow[Sign - 1])
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
    num_A, num_B = read_simulation_flow(ff, attribute_values)

    region_flag = np.zeros((len(neighbors), len(neighbors)))
    cal_flow_flag(attribute_values, region_flag)
    new_region = Ant_amoeba(neighbors, ant_num, iteration_time, EliteAntNum, EvapCoeff,
                            attribute_values, region_flag, num_A, num_B)

    final_regions = delete_overlap_region(new_region)

    final_true_regions = Monte_carlo_simulation(attribute_values, final_regions, rep_num, sig_level, num_A, num_B)

    for i in range(4):
        Sign = i + 1
        TP = 0
        final_num = 0
        for l in final_true_regions:
            flag = cluster_flow[i][0].split('-')
            if int(flag[0]) in l.member_o or int(flag[1]) in l.member_d:
                final_num += len(l.member_flow)
                for n in l.member_flow:
                    flow_key = '{}-{}'.format(n[0], n[1])
                    if flow_key in cluster_flow[i]:
                        TP += 1
            if TP != 0:
                Precision = TP / final_num
                Recall = TP / len(cluster_flow[Sign - 1])
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
