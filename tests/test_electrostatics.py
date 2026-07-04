import numpy as np
from rot_eo_model.electrostatics import build_parallel_plate_mesh, ElectrostaticModel2D, EPSILON_0
from rot_eo_model.tensors import CapacitanceMatrix

def test_parallel_plate_capacitance():
    width = 10e-3
    gap = 1e-3
    eps_r = 44.0
    
    mesh = build_parallel_plate_mesh(width, gap, res=50)
    def eps_func(p):
        return np.full(p.shape[1], eps_r)
        
    model = ElectrostaticModel2D(mesh, eps_func)
    
    electrode_names = ['top', 'bottom']
    C_node_prime = model.extract_capacitance_matrix_per_length(electrode_names)
    
    # Analytical: C' = eps0 * eps_r * width / gap
    C_analytical = EPSILON_0 * eps_r * width / gap
    
    # The diagonal elements should be C_analytical and off-diagonal should be -C_analytical
    assert np.isclose(C_node_prime[0, 0], C_analytical, rtol=0.05)
    assert np.isclose(C_node_prime[0, 1], -C_analytical, rtol=0.05)

def test_differential_transformation():
    # 4x4 dummy C_node matrix
    C_node = np.array([
        [ 2.0, -1.0, -0.1, -0.1],
        [-1.0,  2.0, -0.1, -0.1],
        [-0.1, -0.1,  2.0, -1.0],
        [-0.1, -0.1, -1.0,  2.0],
    ])
    
    T = 0.5 * np.array([
        [ 1,  0],
        [-1,  0],
        [ 0,  1],
        [ 0, -1],
    ])
    
    cap_mat = CapacitanceMatrix(C_node=C_node, T=T)
    assert cap_mat.verify_symmetry()
    
    C_diff = cap_mat.C_diff
    # C_diff[0, 0] = T[:, 0].T @ C_node @ T[:, 0]
    # T[:, 0] = [0.5, -0.5, 0, 0]
    # C_diff[0, 0] = 0.25*(2.0) + 0.25*(2.0) - 2*0.25*(-1.0) = 1.0 + 0.5 = 1.5
    assert np.isclose(C_diff[0, 0], 1.5)
    
    # Energy test
    V_node = np.array([1.0, -1.0, 0.0, 0.0])
    U = cap_mat.electrostatic_energy(V_node)
    assert U > 0
