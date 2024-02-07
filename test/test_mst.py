import pytest
import numpy as np
from mst import Graph
from sklearn.metrics import pairwise_distances


def check_mst(adj_mat: np.ndarray, 
              mst: np.ndarray, 
              expected_weight: int, 
              allowed_error: float = 0.0001):
    """
    
    Helper function to check the correctness of the adjacency matrix encoding an MST.
    Note that because the MST of a graph is not guaranteed to be unique, we cannot 
    simply check for equality against a known MST of a graph. 

    Arguments:
        adj_mat: adjacency matrix of full graph
        mst: adjacency matrix of proposed minimum spanning tree
        expected_weight: weight of the minimum spanning tree of the full graph
        allowed_error: allowed difference between proposed MST weight and `expected_weight`

    TODO: Add additional assertions to ensure the correctness of your MST implementation. For
    example, how many edges should a minimum spanning tree have? Are minimum spanning trees
    always connected? What else can you think of?

    """

    def approx_equal(a, b):
        return abs(a - b) < allowed_error

    total = 0
    for i in range(mst.shape[0]):
        for j in range(i+1):
            total += mst[i, j]
    assert approx_equal(total, expected_weight), 'Proposed MST has incorrect expected weight'

    # check that the number of edges is equal to the number of nodes minus one; the number of edges should be the number of nonzero elements in the matrix divided by two (since it is symmetric)
    assert (np.count_nonzero(mst)/2)==(mst.shape[0]-1)

    # check that the number of nodes in the adjacency matrix is equal to the number of nodes in the mst (we can add one to the number of nonzero edges)
    # this is to check if the adjacency matrix is fully connected, if it isn't, it is a forest instead of a mst
    assert(((np.count_nonzero(mst)/2)+1)==adj_mat.shape[0])


def test_mst_small():
    """
    
    Unit test for the construction of a minimum spanning tree on a small graph.
    
    """
    file_path = './data/small.csv'
    g = Graph(file_path)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 8)


def test_mst_single_cell_data():
    """
    
    Unit test for the construction of a minimum spanning tree using single cell
    data, taken from the Slingshot R package.

    https://bioconductor.org/packages/release/bioc/html/slingshot.html

    """
    file_path = './data/slingshot_example.txt'
    coords = np.loadtxt(file_path) # load coordinates of single cells in low-dimensional subspace
    dist_mat = pairwise_distances(coords) # compute pairwise distances to form graph
    g = Graph(dist_mat)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 57.263561605571695)


def test_mst_student():
    """
    
    TODO: Write at least one unit test for MST construction.
    
    """
    # check that an error is raised when we provide an unsymmetric/directed network
    with pytest.raises(TypeError) as g:
        file_path = './data/small_directed.csv' # created from small.csv
        Graph(file_path)
    assert g.type==TypeError

    # check that an error is raised when we provide a uniform adjacency matrix
    with pytest.raises(TypeError) as g:
        file_path = './data/small_uniform_weights.csv' # created from small.csv
        Graph(file_path)
    assert g.type==TypeError

    # check that our output mst adjacency matrix is symmetric (undirected output)
    file_path = './data/small.csv'
    g = Graph(file_path)
    g.construct_mst()
    assert np.allclose(g.mst, g.mst.T)

    # pass
