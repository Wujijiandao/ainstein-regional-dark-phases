import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from analysis_descendant_transport_history import transported_partition


def test_transported_partition_is_full_and_uses_known_labels():
    n=4;box=4.0;nlab=2
    # two particle groups at opposite x positions
    pos=np.array([[0.5,0.5,0.5],[0.6,0.5,0.5],[2.5,2.5,2.5],[2.6,2.5,2.5]])
    labs=np.array([1,1,2,2],dtype=np.int32)
    grid,occ=transported_partition(pos,labs,n,box,nlab)
    assert grid.shape==(n,n,n)
    assert np.all(grid>0)
    assert set(np.unique(grid)) <= {1,2}
    assert 0<occ<=1


def test_periodic_nearest_fill_wraps():
    n=4;box=4.0;nlab=2
    pos=np.array([[0.1,0.1,0.1],[2.1,2.1,2.1]])
    labs=np.array([1,2],dtype=np.int32)
    grid,_=transported_partition(pos,labs,n,box,nlab)
    # cell at x=3 is periodically closer to x=0 group than x=2 for matching y,z=0
    assert grid[3,0,0]==1
