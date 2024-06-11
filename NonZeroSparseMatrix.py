from scipy.sparse import csr_matrix


class NonZeroSparseMatrix:
    def __init__(self, matrix, default_value=0):
        if not isinstance(matrix, csr_matrix):
            raise ValueError("matrix must be an instance of scipy.sparse.csr_matrix")
        self.matrix = matrix
        self.default_value = default_value

    def __getitem__(self, index):
        value = self.matrix[index]
        # if isinstance(value, csr_matrix.__class__):  # Handle slicing cases
        #     return value
        if value == 0:
            return self.default_value
        return value

    def __setitem__(self, index, value):
        self.matrix[index] = value

    def toarray(self):
        # Doesn't work if done with too large array as it would need too much RAM
        arr = self.matrix.toarray()
        arr[arr == 0] = self.default_value
        return arr

    def to_csr(self):
        return self.matrix
