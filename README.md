#### Requirements
    python 3.6.5
    networkx(2.3)
    numpy(1.14.2)
    pandas(1.2.4)

#### Data instructions
    The road network are stored in the folder 'road network'
    The file 'road.shp' stores shape file of road network
    The file 'edges_neighbors.txt' stores neighboring road segments for each road segment.

### 数据说明 
    道路网络存储在名为「road network」的文件夹中。
    文件 road.shp：存储道路网络的形状文件。
    文件 edges_neighbors.txt：存储每个道路段的相邻道路段信息。


    The synthetic datasets are stored in the folder 'synthetic data'.
    floder 'TypeI' stores all files related to type I flows, floder 'TypeII' stores all files related to type II flows
    The file 'data.csv' stores synthetic dataset.
    The file 'knn_dist.txt' stores network-constrained k-nearest neighbors for each flow, you can use 1_search_flow_knn.py to create this file
    The file 'keep_o_neighbors.txt'/'keep_d_neighbors.txt' stores Mo/Md, you can use 2_prepare_for_fast_mcs.py to create this file 
    The file 'cores_range.txt' stores core flows and its directly reachable distance, you can use 3_identify_core_flows.py to create this file 
    The file 'other_neighbors' stores neighboring flows of the other type for each core flow, you can use 4_search_the_other_flow_neighbors.py to create this file 
    
    合成数据集存储在名为「synthetic data」的文件夹中。
    文件夹「TypeI」：存储所有与类型 I 流量相关的文件。
    文件夹「TypeII」：存储所有与类型 II 流量相关的文件。
    文件 data.csv：存储合成数据集。
    文件 knn_dist.txt：存储每个流的网络约束 k 近邻，可使用脚本 1_search_flow_knn.py 创建。
    文件 keep_o_neighbors.txt/keep_d_neighbors.txt：存储 Mo/Md，可使用脚本 2_prepare_for_fast_mcs.py 创建。
    文件 cores_range.txt：存储核心流及其直接可达距离，可使用脚本 3_identify_core_flows.py 创建。
    文件 other_neighbors：存储每个核心流的另一类型的邻近流，可使用脚本 4_search_the_other_flow_neighbors.py 创建。 </step3_refined_translation>

#### Run instructions
    The codes are are stored in the folder 'BiFlowCluster'.
    The consists of four steps:
        step 1. search network-constrained k-nearest neighbors for each flow by running 1_search_flow_knn.py; // You can skip this step beacuse we've already run this step
        step 2. search Mo/Md for each flow by running 2_prepare_for_fast_mcs.py; // You can skip this step beacuse we've already run this step
        step 3. find high-density flow of each type by 3_identify_core_flows.py; // You can skip this step beacuse we've already run this step
        step 4. search neighboring flows of the other type for each high-density flow by running 4_search_the_other_flow_neighbors.py; // You can skip this step beacuse we've already run this step
        step 5. construct bivariate flow clusters by running 5_create_biflow_clusters.py
    运行 1_search_flow_knn.py 搜索每个流的网络约束 k 近邻。（此步骤可以跳过，因为已完成运行。）
    运行 2_prepare_for_fast_mcs.py 搜索每个流的 Mo/Md。（此步骤可以跳过，因为已完成运行。）
    运行 3_identify_core_flows.py 找出每种类型的高密度流。（此步骤可以跳过，因为已完成运行。）
    运行 4_search_the_other_flow_neighbors.py 搜索每个高密度流的另一类型的邻近流。（此步骤可以跳过，因为已完成运行。）
    运行 5_create_biflow_clusters.py 构建双变量流量簇。

#### Output instructions
    each detected bivariate flow cluster is stroed in a single file, the file name starts with its cluster type, for example 'HH_1.txt' 
    /////////////////////HH_1.txt//////////////////////////
    data format:
    I: flowID1, flowID2,....,flowIDn. // this line only contains flow ID of type I
    Ii: flowID1, flowID2,....,flowIDn. // this line only contains flow ID of type II
    I:0,1,2,3,4,5,6,7,8,9
    II:0,1,2,3,4,5,6,7,8,9
