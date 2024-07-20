import numpy as np
from ortools.linear_solver import pywraplp
import sys
import scanpy as sc
import time
from scipy.spatial import KDTree
import math

def extract_coords(adata):
    x = adata.obs['CCF_AP_axis'].values
    y = adata.obs['CCF_ML_axis'].values
    z = adata.obs['CCF_DV_axis'].values
    return np.vstack((x, y, z)).T


def weighted_distance(distance, scale, max_distance):
    # Apply exponential growth function
    return max_distance * (1 - math.exp(-distance / (scale * max_distance)))


def distance_matrix(adata1, adata2, num_neighbours):
    print(f"Started distance calculations for region {adata1[0].obs['CCFname']}...")

    coords1 = extract_coords(adata1)
    coords2 = extract_coords(adata2)
    kdtree2 = KDTree(coords2)
    num_valid_neighbours = num_neighbours
    if num_neighbours > len(coords2): num_valid_neighbours = len(coords2)
    distance_to_neighbours, neighbour_indices = kdtree2.query(coords1, k=num_valid_neighbours)
    distance_thresholds = distance_to_neighbours[:, -1]

    expression_matrix1 = adata1.X.toarray()
    expression_matrix2 = adata2.X.toarray()

    distances_matrix = np.full((len(coords1), len(coords2)), 2147483647, dtype=np.uint32)

    print("Calculating Distances...")
    time2 = time.time()
    for cell_idx in range(0, len(coords1)):
        neighbours_physical_matrix = coords2[neighbour_indices[cell_idx]]
        neighbours_expression_matrix = expression_matrix2[neighbour_indices[cell_idx]]

        physical_vector_norm = np.sum(coords1[cell_idx] ** 2)
        physical_norm2 = np.sum(neighbours_physical_matrix ** 2, axis=1)
        physical_dot_product = np.dot(neighbours_physical_matrix, coords1[cell_idx])
        physical_distances = np.sqrt(physical_vector_norm + physical_norm2 - 2 * physical_dot_product)

        expression_vector_norm = np.sum(expression_matrix1[cell_idx] ** 2)
        expression_norm2 = np.sum(neighbours_expression_matrix ** 2, axis=1)
        expression_dot_product = np.dot(neighbours_expression_matrix, expression_matrix1[cell_idx])
        expression_distances = np.sqrt(expression_vector_norm + expression_norm2 - 2 * expression_dot_product)

        for idx in range(0, len(neighbour_indices[cell_idx])):
            spatial_distance = physical_distances[idx]
            expression_distance = expression_distances[idx]
            if spatial_distance > 0 and expression_distance > 0:
                adjusted_spt_distance = weighted_distance(spatial_distance, 1, (distance_thresholds[cell_idx])) * 2
                adjusted_expression_distance = expression_distance * 100
                distances_matrix[cell_idx, neighbour_indices[cell_idx][idx]] = np.add(adjusted_spt_distance,
                                                                                      adjusted_expression_distance)

                # # Print results for debugging
                # print(f"\nCell {cell_idx} from adata1 to Cell {neighbour_indices[cell_idx][idx]} from adata2:\n"
                #     f"Calced Distance: {distances_matrix[cell_idx, neighbour_indices[cell_idx][idx]]}\n"
                #     # f"Exp Coords1: {','.join(map(str, expression_matrix1[cell_idx]))}\n"
                #     # f"Exp Coords2: {','.join(map(str, expression_matrix2[neighbour_indices[cell_idx, idx]]))}\n"
                #     f"Calced Exp Distance: {expression_distances[idx]}\n"
                #     f"Adjusted Exp Distance: {adjusted_expression_distance}\n"
                #     f"Spt Coords1: {','.join(map(str, coords1[cell_idx]))}\n"
                #     f"Spt Coords2: {','.join(map(str, coords2[neighbour_indices[cell_idx][idx]]))}\n"
                #     f"Calced Spt Distance: {physical_distances[idx]}\n"
                #     f"Adjusted Spt Distance: {adjusted_spt_distance}\n")

    print(f"Time to Calculate Euclidian Distances: {time.time() - time2}s")

    return distances_matrix


def create_data_model(adata1, adata2, adata3, adata4):
    # Example cost matrices for 4 datasets with 3 samples each MUST BE SQUARE FOR NOW
    cost_matrix_1_2 = np.random.rand(3, 3)
    cost_matrix_1_3 = np.random.rand(3, 3)
    cost_matrix_1_4 = np.random.rand(3, 3)
    cost_matrix_2_3 = np.random.rand(3, 3)
    cost_matrix_2_4 = np.random.rand(3, 3)
    cost_matrix_3_4 = np.random.rand(3, 3)

    size = cost_matrix_1_2.shape[0]
    combined_distances_matrix = np.zeros((size, size, size, size))

    for i in range(size):
        for j in range(size):
            for k in range(size):
                for l in range(size):
                    combined_distances_matrix[i, j, k, l] = (
                            cost_matrix_1_2[i, j] +
                            cost_matrix_1_3[i, k] +
                            cost_matrix_1_4[i, l] +
                            cost_matrix_2_3[j, k] +
                            cost_matrix_2_4[j, l] +
                            cost_matrix_3_4[k, l]
                    )

    return combined_distances_matrix, size


def solve_assignment(combined_distances_matrix, size):
    # Create the mip solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SCIP')

    # Create decision variables
    x = {}
    for i in range(size):
        for j in range(size):
            for k in range(size):
                for l in range(size):
                    x[i, j, k, l] = solver.BoolVar(f'x[{i},{j},{k},{l}]')

    # Define the objective function
    solver.Minimize(solver.Sum(combined_distances_matrix[i, j, k, l] * x[i, j, k, l]
                               for i in range(size) for j in range(size) for k in range(size) for l in range(size)))

    # Constraints: each sample should be assigned exactly once in each dimension
    for i in range(size):
        solver.Add(solver.Sum(x[i, j, k, l] for j in range(size) for k in range(size) for l in range(size)) == 1)
    for j in range(size):
        solver.Add(solver.Sum(x[i, j, k, l] for i in range(size) for k in range(size) for l in range(size)) == 1)
    for k in range(size):
        solver.Add(solver.Sum(x[i, j, k, l] for i in range(size) for j in range(size) for l in range(size)) == 1)
    for l in range(size):
        solver.Add(solver.Sum(x[i, j, k, l] for i in range(size) for j in range(size) for k in range(size)) == 1)

    # Solve the problem
    status = solver.Solve()

    if status == pywraplp.Solver.OPTIMAL:
        print('Solution:')
        assignments = []
        for i in range(size):
            for j in range(size):
                for k in range(size):
                    for l in range(size):
                        if x[i, j, k, l].solution_value() == 1:
                            assignments.append((i, j, k, l))
                            print(f'Sample {i} from Dataset A assigned to Sample {j} from Dataset B, '
                                  f'Sample {k} from Dataset C, Sample {l} from Dataset D')
        return assignments
    else:
        print('The problem does not have an optimal solution.')


# Load AnnData objects
adata1_path = sys.argv[1]
adata2_path = sys.argv[2]
adata3_path = sys.argv[3]
adata4_path = sys.argv[4]

full_adata1 = sc.read_h5ad(adata1_path)
full_adata2 = sc.read_h5ad(adata2_path)
full_adata3 = sc.read_h5ad(adata3_path)
full_adata4 = sc.read_h5ad(adata4_path)

num_neighbours = 150

combined_distances_matrix, size = create_data_model(full_adata1, full_adata2, full_adata3, full_adata4)
solve_assignment(combined_distances_matrix, size)


