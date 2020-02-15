import os
import numpy as np
from tridiag_mat import inversion as inv
import nose.tools as nt

def build_matrix(n):
    mat = np.zeros((n, n))
    i = np.arange(n - 1)
    print (i)
    mat [i, i ] = 1
    mat[i, i+1] = 0.5
    mat[i+1, i] = 0.4
    mat[-1, -1] = 1
    return mat

def test_theta():
    mat = build_matrix(10)
    theta = inv.theta(mat, 0)
    nt.assert_equal(theta, 1.)
    theta = inv.theta(mat, 1)
    nt.assert_equal(theta, mat[0, 0])
    theta = inv.theta(mat, 2)
    nt.assert_equal(theta, 0.8)

def test_phi():
    mat = build_matrix(10)
    phi = inv.phi(mat, 10)
    nt.assert_equal(phi, 1.)
    phi = inv.phi(mat, 11)
    nt.assert_equal(phi, 1.)
    phi = inv.phi(mat, 9)
    nt.assert_almost_equal(phi, 0.8)

def test_get_inverse_elm():
    mat = build_matrix(10)
    mat_inv = np.linalg.inv(mat)
    mat_inv_00 = inv.get_inverse_elm(mat, 0, 0)
    nt.assert_almost_equal(mat_inv_00, mat_inv[0,0])
    mat_inv_10 = inv.get_inverse_elm(mat, 1, 0)
    nt.assert_almost_equal(mat_inv_10, mat_inv[1,0])
    mat_inv_13 = inv.get_inverse_elm(mat, 1, 3)
    nt.assert_almost_equal(mat_inv_13, mat_inv[1,3])

def test_get_inverse_mat():
    mat = build_matrix(10)
    mat_inv = np.linalg.inv(mat)
    Minv = inv.get_inverse_mat(mat)
    np.testing.assert_almost_equal(Minv, mat_inv)
