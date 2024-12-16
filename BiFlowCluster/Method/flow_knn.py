# -*- encoding=utf-8 -*-
#
#
class FlowKNN:
    def __init__(self, edges_neighbors):
        self._edges_neighbors = edges_neighbors # 存储路段的邻居（边近邻）

        self._o_road_contains = None
        self._d_road_contains = None # 二维列表，第i个列表存储路段编号为i的路段上的起点/终点的索引
        self._o_road_buffer_points = None
        self._d_road_buffer_points = None
        self._op_info = None
        self._dp_info = None # 列表，每个元素是字典，每个字典3个索引sd：、ed：、road_id:
        self._buffer_size = 0

    def _create_point_info(self, data):
        """
        找到每条路段（边）上的所有起点/终点，存储到_op_info.append/_od_info.append中
        其中每个由点到路段起点的距离、点到路段终点的距离、路段的id组成
        通过路段id（rid）可以查找到位于rid上所有点在_op_info.append/_od_info.append上的索引
        :param data:
        :return:
        """
        self._op_info = list()
        self._dp_info = list()
        self._o_road_contains = [list() for i in range(len(self._edges_neighbors))]
        self._d_road_contains = [list() for i in range(len(self._edges_neighbors))]
        for i in range(len(data)):
            rid, sd, ed = data['o_rid'][i], data['o_sd'][i], data['o_ed'][i]
            p_info = dict()
            p_info['sd'] = sd # 点到路段起点的距离
            p_info['ed'] = ed # 点到路段终点的距离
            p_info['road_id'] = rid
            self._op_info.append(p_info)
            self._o_road_contains[rid].append(i)

            rid, sd, ed = data['d_rid'][i], data['d_sd'][i], data['d_ed'][i]
            p_info = dict()
            p_info['sd'] = sd
            p_info['ed'] = ed
            p_info['road_id'] = rid
            self._dp_info.append(p_info)
            self._d_road_contains[rid].append(i)

    def _get_road_buffer_points(self, rid, point_type):
        """
        通过边近邻确定输入rid的近邻上的所有点
        :param rid:
        :param point_type:
        :return:
        """
        points = list()
        if point_type == 'o':
            for rid in self._edges_neighbors[rid]:
                points.extend(self._o_road_contains[rid])
        else:
            for rid in self._edges_neighbors[rid]:
                points.extend(self._d_road_contains[rid])
        return points

    def _search_buffer_flow_neighbor(self, flow_id):
        """
        确定其起点和终点的缓冲点
        返回两者的交集
        :param flow_id:
        :return:
        """
        o_rid, d_rid = self._op_info[flow_id]['road_id'], self._dp_info[flow_id]['road_id']
        o_buffer_points = set(self._get_road_buffer_points(o_rid, 'o'))
        d_buffer_points = set(self._get_road_buffer_points(d_rid, 'd'))
        buffer_flow_neighbors = o_buffer_points.intersection(d_buffer_points)
        return list(buffer_flow_neighbors) # 当一个flow的起点和终点同时在起点/终点缓冲区，则证明该flow在rid的近邻

    def _calc_point_dists(self, center, others, point_type):
        if point_type == 'o':
            points_info = self._op_info
        else:
            points_info = self._dp_info
        c_rid = points_info[center]['road_id']
        c_sd, c_ed = points_info[center]['sd'],  points_info[center]['ed']
        dists = list()
        for pid in others:
            o_rid = points_info[pid]['road_id']
            o_sd, o_ed = points_info[pid]['sd'], points_info[pid]['ed']
            if o_rid == c_rid:
                dists.append(abs(c_sd - o_sd))
                continue
            # c点到c_rid的起点+c_rid的起点到o_rid的起点+o点到o_rid的起点
            d1 = c_sd + self._edges_neighbors[c_rid][o_rid]['ss'] + o_sd
            #
            d2 = c_sd + self._edges_neighbors[c_rid][o_rid]['se'] + o_ed
            #
            d3 = c_ed + self._edges_neighbors[c_rid][o_rid]['es'] + o_sd
            #
            d4 = c_ed + self._edges_neighbors[c_rid][o_rid]['ee'] + o_ed
            dists.append(min(d1, d2, d3, d4))
        return dists

    def _calc_flow_neighbors_dist(self, flow_id, neighbors):
        o_dists = self._calc_point_dists(flow_id, neighbors, 'o')
        d_dists = self._calc_point_dists(flow_id, neighbors, 'd')
        dists = list()
        for i in range(len(neighbors)):
            if o_dists[i] >= self._buffer_size or d_dists[i] >= self._buffer_size:
                continue
            flow_dist = (o_dists[i] + d_dists[i]) / 2
            dists.append([neighbors[i], flow_dist, o_dists[i], d_dists[i]])
        return dists

    def _search_buffer_knn(self, k):
        flow_buffer_neighbors = list()
        dists = list()
        for flow_id in range(len(self._op_info)):
            buffer_flow_neighbor = self._search_buffer_flow_neighbor(flow_id)
            buffer_flow_neighbor.remove(flow_id)
            flow_dists = self._calc_flow_neighbors_dist(flow_id, buffer_flow_neighbor)
            flow_dists = sorted(flow_dists, key=lambda a: a[1])
            num = min(k, len(flow_dists))
            flow_neighbors = list()
            dist = list()
            if num > 0:
                for j in range(num):
                    flow_neighbors.append(flow_dists[j][0])
                    dist.append(flow_dists[j][1:])
            flow_buffer_neighbors.append(flow_neighbors)
            dists.append(dist)
        return flow_buffer_neighbors, dists

    def run(self, data, k, buffer_size):
        self._buffer_size = buffer_size
        self._create_point_info(data)
        flow_knn, dists = self._search_buffer_knn(k)
        return flow_knn, dists
