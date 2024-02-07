import numpy as np
import heapq
from typing import Union

class Graph:

    def __init__(self, adjacency_mat: Union[np.ndarray, str]):
        """
    
        Unlike the BFS assignment, this Graph class takes an adjacency matrix as input. `adjacency_mat` 
        can either be a 2D numpy array of floats or a path to a CSV file containing a 2D numpy array of floats.

        In this project, we will assume `adjacency_mat` corresponds to the adjacency matrix of an undirected graph.
    
        """
        if type(adjacency_mat) == str:
            self.adj_mat = self._load_adjacency_matrix_from_csv(adjacency_mat)
        elif type(adjacency_mat) == np.ndarray:
            self.adj_mat = adjacency_mat
        else: 
            raise TypeError('Input must be a valid path or an adjacency matrix')
        self.mst = None

    def _load_adjacency_matrix_from_csv(self, path: str) -> np.ndarray:
        with open(path) as f:
            return np.loadtxt(f, delimiter=',')

    def construct_mst(self):
        """
    
        TODO: Given `self.adj_mat`, the adjacency matrix of a connected undirected graph, implement Prim's 
        algorithm to construct an adjacency matrix encoding the minimum spanning tree of `self.adj_mat`. 
            
        `self.adj_mat` is a 2D numpy array of floats. Note that because we assume our input graph is
        undirected, `self.adj_mat` is symmetric. Row i and column j represents the edge weight between
        vertex i and vertex j. An edge weight of zero indicates that no edge exists. 
        
        This function does not return anything. Instead, store the adjacency matrix representation
        of the minimum spanning tree of `self.adj_mat` in `self.mst`. We highly encourage the
        use of priority queues in your implementation. Refer to the heapq module, particularly the 
        `heapify`, `heappop`, and `heappush` functions.

        """

        # initialize total cost, minheap, and mst
        total_cost=[]
        minheap=[]
        mst_path=[]

        # get the set of nodes and pick the first node as the starting node, add to the minheap where the first key represents cost (0 for start), and the value represents the node
        Nodes=list(range(len(self.adj_mat)))
        start=Nodes[0]
        heapq.heappush(minheap, (0, start))

        # run while the minheap is not empty
        while len(minheap)!=0:

            # get the cheapest node and check that it is not in our mst
            cost, node=heapq.heappop(minheap)
            if node not in mst_path:

                # if node has not yet been added to mst, add node and its cost
                mst_path.append(node) 
                total_cost.append(cost)

                # for the given adjacent nodes, check if it is in the mst, if not, add to the minheap and repeat
                for adj_node, cost in zip(Nodes, self.adj_mat[node]):

                    # check that adjacent node is not in mst and that the edge exists
                    if (adj_node not in mst_path) and (cost!=0):
                        heapq.heappush(minheap, (cost, adj_node))

        # create adjacency matrix coordinates by taking adjacent pairs of nodes from our mst
        mst_node_pairs=[(i,j) for i,j in zip(mst_path, mst_path[1:]+[mst_path[0]])] # creating node pairs (https://www.geeksforgeeks.org/python-pair-iteration-in-list/)

        # initialize empty mst adjacency matrix and populate with node pairs
        mst=[[0]*len(Nodes) for _ in range(len(Nodes))] # https://stackoverflow.com/questions/13157961/2d-array-of-zeros
        for node,c in zip(mst_node_pairs, total_cost):

            # symmetric, so add both directions of node pairs
            mst[node[0]][node[1]]=c
            mst[node[1]][node[0]]=c

        self.mst=np.array(mst)

# g=Graph("data/small.csv")
# g.adj_mat

# # initialize total cost, minheap, and mst
# total_cost=[]
# minheap=[]
# mst_path=[]

# # get the set of nodes and pick the first node as the starting node, add to the minheap where the first key represents cost (0 for start), and the value represents the node
# Nodes=list(range(len(g.adj_mat)))
# start=Nodes[0]
# heapq.heappush(minheap, (0, start))

# # run while the minheap is not empty
# while len(minheap)!=0:

#     # get the cheapest node and check that it is not in our mst
#     cost, node=heapq.heappop(minheap)
#     if node not in mst_path:

#         # if node has not yet been added to mst, add node and its cost
#         mst_path.append(node) 
#         total_cost.append(cost)

#         # for the given adjacent nodes, check if it is in the mst, if not, add to the minheap and repeat
#         for adj_node, cost in zip(Nodes, g.adj_mat[node]):

#             # check that adjacent node is not in mst and that the edge exists
#             if (adj_node not in mst_path) and (cost!=0):
#                 heapq.heappush(minheap, (cost, adj_node))

# # create adjacency matrix coordinates by taking adjacent pairs of nodes from our mst
# mst_node_pairs=[(i,j) for i,j in zip(mst_path, mst_path[1:]+[mst_path[0]])] # creating node pairs (https://www.geeksforgeeks.org/python-pair-iteration-in-list/)

# # initialize empty mst adjacency matrix and populate with node pairs
# mst=[[0]*len(Nodes) for _ in range(len(Nodes))] # https://stackoverflow.com/questions/13157961/2d-array-of-zeros
# for node,c in zip(mst_node_pairs, total_cost):
#     mst[node[0]][node[1]]=c
#     mst[node[1]][node[0]]=c


# mst
# mst=np.array(mst)

# mst

# file_path = './data/small.csv'
# g = Graph(file_path)
# g.construct_mst()
# g.mst


# def check_mst(adj_mat: np.ndarray, 
#               mst: np.ndarray, 
#               expected_weight: int, 
#               allowed_error: float = 0.0001):

#     def approx_equal(a, b):
#         return abs(a - b) < allowed_error


#     total = 0
#     for i in range(mst.shape[0]):
#         for j in range(i+1):
#             total += mst[i, j]
#     return total


# check_mst(g.adj_mat, mst, 8)

# total=0
# for i in range(mst.shape[0]):
#     for j in range(i+1):
#         # total += mst[i, j]
#         print(i,j)
#         print(mst[i, j])
#         total+=mst[i, j]


# print(total)

# for i in range(mst.shape[0]):
#     for j in range(i+1):
#         print(i,j)