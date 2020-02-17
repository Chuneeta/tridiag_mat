import numpy as np

def theta(mat, i):
    r, c = mat.shape
    assert r == c, 'Square matrix is required'
    if i == 0: return 1.
    elif i == 1: return mat[i-1, i-1]
    else: return mat[i-1, i-1] * theta(mat, i-1) - mat[i - 2 , i - 1] * mat[i - 1, i - 2] * theta(mat, i-2)

def phi(mat, i):
    r, c = mat.shape
    assert r == c, 'Square matrix is required'
    n = r
    if i == n: 
        return mat[-1, -1]
    elif i == n + 1: 
        return 1.
    else:
        return mat[i, i] * phi(mat, i + 1) - mat[i - 1 , i] * mat[i, i - 1] * phi(mat, i + 2)

def get_inverse_elm(mat, i, j):
    r, c = mat.shape
    assert r == c, 'Square matrix is required'
    ii = i + 1
    jj = j + 1
    if i < j:
        k = np.arange(i, j) 
        prod = np.prod(mat[k, k+1])
        T_ij = ((-1)**(i + j) * prod * theta(mat, ii - 1) * phi(mat, jj + 1) ) / theta(mat, r)
    if i == j:
        T_ij = (theta(mat, ii - 1) * phi(mat, jj + 1)) / theta(mat, r)
    if i > j:
        k = np.arange(j, i)
        prod = np.prod(mat[k + 1, k])
        T_ij = ((-1)**(i + j) * prod * theta(mat, jj - 1) * phi(mat, ii + 1) ) / theta(mat, r)
    
    return T_ij

def _find_convergence(mat, inv_mat, max_iter=25, tol=1e-11):
    I_true = np.identity(len(mat))
    I_est = np.dot(mat, inv_mat)
    chisq = np.sum(I_est - I_true)**2
    iter_n = 0
    while chisq > tol or iter_n < max_iter:
        I_est_inv = _eval_approx_inverse(I_est)
        I_est = np.dot(I_est, I_est_inv)
        iter_n += 1
    return np.dot(inv_mat, I_est)

def _eval_approx_inverse(mat):
    r, c = mat.shape
    M_inv = np.zeros((r, c))
    i = np.arange(c - 1)
    M_inv[i, i] = 1 / mat[i, i]
    M_inv[i, i + 1] = -1 * mat[i, i + 1] / (mat[i, i] * mat[i + 1, i + 1])
    M_inv[i + 1, i] = -1 * mat[i + 1, i] / (mat[i, i] * mat[i + 1, i + 1])
    # filling in last element
    M_inv[c-1, c-1] = get_inverse_elm(mat, c-1, c-1)
    
    return M_inv

def get_approx_inverse(mat):
    r, c = mat.shape
    assert r == c, 'Square matrix is required'
    M_inv = _eval_approx_inverse(mat)
    M_inv = _find_convergence(mat, M_inv, tol=1e-3)
    
    return M_inv

def get_inverse_mat(mat):
    r, c = mat.shape
    assert r == c, 'Square matrix is required'
    # decomposing into upper and lower triangular matrices
    tril_ind = np.tril_indices(mat.shape[0], k=-1)
    triu_ind = np.triu_indices(mat.shape[0], k=1)
    diag_ind = np.diag_indices(mat.shape[0])
    # evaluating the inverse of each element
    Minv = np.zeros((r, c))
    Minv[tril_ind] = [get_inverse_elm(mat, i, j) for (i,j) in zip(tril_ind[0], tril_ind[1])]
    Minv[triu_ind] = [get_inverse_elm(mat, i, j) for (i,j) in zip(triu_ind[0], triu_ind[1])]
    Minv[diag_ind] = [get_inverse_elm(mat, i, i) for i in np.arange(mat.shape[0])]
   
    return Minv

