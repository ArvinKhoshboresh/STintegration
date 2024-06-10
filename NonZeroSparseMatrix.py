import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix


class NonZeroSparseMatrix:
    # Wrapper for sparse matrices to store sparse values while returning some default value instead of zero for
    # non-stored indices
    # Takes some default value to return instead of zero, an existing matrix, or shape to create new matrix

    def __init__(self, default_value, sparse_matrix=None, shape=None):
        self.default_value = default_value
        if sparse_matrix is not None:
            self.sparse_matrix = sparse_matrix
        elif shape is not None:
            self.sparse_matrix = coo_matrix(shape)
        else:
            raise ValueError("Must provide either a sparse_matrix or a shape")

    def to_csr(self):
        return self.sparse_matrix.tocsr()

    def to_csc(self):
        return self.sparse_matrix.tocsc()

    def to_coo(self):
        return self.sparse_matrix.tocoo()

    def __getitem__(self, indices):
        row, col = indices
        try:
            return self.sparse_matrix[row, col] if (row, col) in zip(self.sparse_matrix.row,
                                                                     self.sparse_matrix.col) else self.default_value
        except AttributeError:  # Handle formats other than COO
            if row < self.sparse_matrix.shape[0] and col < self.sparse_matrix.shape[1]:
                return self.sparse_matrix[row, col] if self.sparse_matrix[row, col] != 0 else self.default_value
            else:
                raise IndexError("Index out of bounds")

    # No idea if this works
    def __setitem__(self, indices, value):
        row, col = indices
        if value != self.default_value:
            if isinstance(self.sparse_matrix, coo_matrix):
                self.sparse_matrix = coo_matrix(
                    (np.append(self.sparse_matrix.data, value),
                     (np.append(self.sparse_matrix.row, row),
                      np.append(self.sparse_matrix.col, col))),
                    shape=self.sparse_matrix.shape)
            else:
                self.sparse_matrix[row, col] = value
        else:
            if isinstance(self.sparse_matrix, coo_matrix):
                mask = ~((self.sparse_matrix.row == row) & (self.sparse_matrix.col == col))
                self.sparse_matrix = coo_matrix(
                    (self.sparse_matrix.data[mask],
                     (self.sparse_matrix.row[mask],
                      self.sparse_matrix.col[mask])),
                    shape=self.sparse_matrix.shape)
            else:
                self.sparse_matrix[row, col] = 0

    # def __str__(self):
    #     dense_repr = np.full(self.sparse_matrix.shape, self.default_value)
    #     if isinstance(self.sparse_matrix, coo_matrix):
    #         for r, c, v in zip(self.sparse_matrix.row, self.sparse_matrix.col, self.sparse_matrix.data):
    #             dense_repr[r, c] = v
    #     else:
    #         dense_repr += self.sparse_matrix.toarray()
    #     return str(dense_repr)
