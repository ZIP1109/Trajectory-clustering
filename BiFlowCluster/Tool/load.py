# -*- encoding=utf-8 -*-
#
# functions to load data

def load_road_neighbors(rnf):
    """
    从指定文件中加载道路邻接信息
    :param rnf:
    :return:
    """
    road_neighbors = list()
    f = open(rnf)
    for i in f:
        i = i.rstrip('\n').split(';')
        road_neighbor = dict()
        for j in i:
            j = j.split(',')
            rid = int(j[0])
            ss = round(float(j[1]), 3)
            se = round(float(j[2]), 3)
            es = round(float(j[3]), 3)
            ee = round(float(j[4]), 3)
            road_neighbor[rid] = {'ss': ss, 'se': se, 'es': es, 'ee': ee, 'dist': [ss, se, es, ee]}
        road_neighbors.append(road_neighbor)
    return road_neighbors