import numpy as np
import scipy.sparse as sp

dim = 50000
mat = sp.dia_array((np.ones(1), 0), shape=(dim,dim),dtype=bool).tocsr()
# res = sp.block_array([[mat for i in range(50000)]])
# print(res.shape)
wide = sp.block_diag(mat).tocoo()
mult = mat.tocoo() @ wide
print(mult.shape)
