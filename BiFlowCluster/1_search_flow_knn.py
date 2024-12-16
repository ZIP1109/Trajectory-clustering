# -*- encoding=utf-8 -*-
#
# searching flow knn using edge neighbors

from os import path

import pandas as pd
import os
import time

from Tool.load import load_road_neighbors
from Method.flow_knn import FlowKNN


""""
def neighbors_2_file(nf, knn_s, dists_s):
    wf = open(nf, 'w')
    for knn, dists in zip(knn_s, dists_s):
        s = ['{0},{1},{2},{3}'.format(knn[i], round(dists[i][0], 3), round(dists[i][1], 3), round(dists[i][2], 3),)
             for i in range(len(knn))]
        wf.write(';'.join(s) + '\n')
    wf.close()
"""


def neighbors_2_file(nf, knn_s, dists_s):
    """
    将邻居索引及其对应的距离写入文件。如果目标目录不存在，则自动创建。

    参数:
    - nf (str): 要写入数据的文件名（包含路径）。
    - knn_s (可迭代对象): 每个数据点的邻居索引列表。
    - dists_s (可迭代对象): 对应于每个邻居的距离列表。
    """
    # 获取文件所在的目录路径
    directory = os.path.dirname(nf)

    # 如果目录路径不为空，则创建目录（exist_ok=True 表示如果目录已存在则不报错）
    if directory:
        os.makedirs(directory, exist_ok=True)

    # 使用 with 语句确保文件正确关闭
    with open(nf, 'w') as wf:
        for knn, dists in zip(knn_s, dists_s):
            # 使用列表推导式和 f-字符串生成格式化字符串
            s = [
                f"{knn[i]},{round(dists[i][0], 3)},{round(dists[i][1], 3)},{round(dists[i][2], 3)}"
                for i in range(len(knn))
            ]
            # 将字符串列表用分号连接并写入文件
            wf.write(';'.join(s) + '\n')

def run():
    root_dir = 'F:/论文/轨迹聚类/OD流/A network-constrained clustering method/BiFlowCluster&Code&Data'

    road_dir = path.join(root_dir, 'road network')
    road_neighbor_file = path.join(road_dir, 'edges_neighbors.txt')
    road_neighbors = load_road_neighbors(road_neighbor_file)
    flow_knn = FlowKNN(road_neighbors)
    k = 15

    for i in range(1,101):
        for type in ['I','II']:
            root_dir = 'F:/论文/轨迹聚类/OD流/A network-constrained clustering method/BiFlowCluster&Code&Data'
            my_test_address = 'F:/论文/轨迹聚类/OD流/A network-constrained clustering method/mytest'

            data_dir = path.join(root_dir, 'synthetic data/'+str(i)+'/Type'+type)
            my_test_address = path.join(my_test_address, 'synthetic data/'+str(i)+'/Type'+type)
            data_file = path.join(data_dir, 'data.csv')
            data = pd.read_csv(data_file)
            knn, dists = flow_knn.run(data, k, buffer_size=250)
            kf = path.join(my_test_address, 'knn_dist.txt')
            neighbors_2_file(kf, knn, dists)

if __name__ == '__main__':
    start_time = time.time()
    run()
    print(f'运行时间：{time.time() - start_time}秒')