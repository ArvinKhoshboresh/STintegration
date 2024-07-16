# import numpy as np
# from scipy.sparse import random as sparse_random
# import sslap
#
# size = 1000
# density = 0.01
# loop_count = 100
#
# for i in range(loop_count):
#     sparse_matrix = sparse_random(size, size, density=density, format='coo', data_rvs=np.random.rand)
#     print(sparse_matrix)
#     result = sslap.auction_solve(coo_mat=sparse_matrix, cardinality_check=False)
#     print(result["sol"])
#     print(f"Loop {i + 1} completed")

import numpy as np
from scipy.sparse import random as sparse_random
import sslap

size = 1000
density = 0.01
loop_count = 100  # Adjust loop count as needed

for i in range(loop_count):
    dense_matrix = sparse_random(size, size, density=density, format='coo', data_rvs=np.random.rand).todense()
    result = sslap.auction_solve(dense_matrix, cardinality_check=False, fast=True)
    print(f"Loop {i + 1} completed")

