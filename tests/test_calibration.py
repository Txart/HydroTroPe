import pytest
import numpy as np

from ..calibration import mse_two_vectors_with_nans

 
def test_mse_two_vectors_with_nans():
    v = np.array([1, np.nan, 1])
    w1 = np.array([0, 0, np.nan])
    w2 = np.array([0, np.nan, np.nan])
    
    np.testing.assert_allclose(mse_two_vectors_with_nans(v, w1), 1)
    np.testing.assert_allclose(mse_two_vectors_with_nans(v, w2), 1)