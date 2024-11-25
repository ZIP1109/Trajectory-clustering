The codes of these three methods were developed based on JetBrains PyCharm Community Edition 2019.2.3 x64, and the language is Python 3.8.1.

//////////////////////////////////////////////////////////////Implement of BiFlowCluster method///////////////////////////////////////////////////////////////////////////////////////////////

BiFlowCluster.py was used to detect clustering results based on  BiFlowCluster method
// running cost
               The running cost of BiFlowAMOEBA.py is 40 seconds for synthetic data.


////////////////////////////////////////////////////////////////////BiFlowCluster.py////////////////////////////////////////////////////////////////////////////////////////////////////

 data_dir =  'C:/Users/12449/Desktop/IJGIS-code&data(2021)/Synthetic data/D{}'.format(i) # the main path of all files
//the paths of input file
	nf = path.join(data_dir, 'neighbors.txt')  # the path of file that stores neighbors of each polygon
    	ff = path.join(data_dir, 'flow.txt')) # the path of file that stores bivariate flow datasets
    	cf = path.join(data_dir, 'Cluster_flow.txt') # the path of file that stores flows of clusters
   
//the input parameters
	rep_num = 999 # the number of Monte Carlo simulations
    	sig_level = 0.01 # the significant level


//////////////////////////////////////////////////////////////Implement of AntScan_flow method///////////////////////////////////////////////////////////////////////////////////////////////

AntScan_flow.py was used to detect clustering results based on  AntScan_flow method
// running cost
               The running cost of BiFlowCluster.py is 209 seconds for synthetic data.


////////////////////////////////////////////////////////////////////AntScan_flow.py////////////////////////////////////////////////////////////////////////////////////////////////////

 data_dir =  'C:/Users/12449/Desktop/IJGIS-code&data(2021)/Synthetic data/D{}'.format(i) # the main path of all files
//the paths of input file
	nf = path.join(data_dir, 'neighbors.txt')  # the path of file that stores neighbors of each polygon
    	ff = path.join(data_dir, 'flow.txt')) # the path of file that stores bivariate flow datasets
    	cf = path.join(data_dir, 'Cluster_flow.txt') # the path of file that stores flows of clusters

//the input parameters
	rep_num = 999 # the number of Monte Carlo simulations
    	sig_level = 0.01 # the significant level
	ant_num = 200 # the number of ants
  	iteration_time = 100 # the number of iterations
    	EliteAntNum = 30 # the number of elite ants
    	EvapCoeff = 0.1 # the evaporation coefficient 

//////////////////////////////////////////////////////////////Implement of BiFlowLISA method///////////////////////////////////////////////////////////////////////////////////////////////

BiFlowLISA.py was used to detect clustering results based on  BiFlowLISA method
// running cost
               The running cost of BiFlowLISA.py is 80 seconds for synthetic data.

////////////////////////////////////////////////////////////////////BiFlowLISA.py////////////////////////////////////////////////////////////////////////////////////////////////////

 data_dir =  'C:/Users/12449/Desktop/IJGIS-code&data(2021)/Synthetic data/D{}'.format(i) # the main path of all files
//the paths of input file
	nf = path.join(data_dir, 'neighbors.txt')  # the path of file that stores neighbors of each polygon
    	ff = path.join(data_dir, 'flow.txt')) # the path of file that stores bivariate flow datasets
    	cf = path.join(data_dir, 'Cluster_flow.txt') # the path of file that stores flows of clusters

//the input parameters
	rep_num = 999 # the number of Monte Carlo simulations
    	sig_level = 0.01 # the significant level


