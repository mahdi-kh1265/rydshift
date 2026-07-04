import numpy as np
from typing import Tuple, List, Callable
from skfem import MeshTri, Basis, ElementTriP1, BilinearForm, LinearForm, solve, condense
from skfem.helpers import dot, grad
from rot_eo_model.tensors import CapacitanceMatrix, FieldPerVoltMatrix, DielectricTensor
import math

# Use vacuum permittivity in F/m
EPSILON_0 = 8.8541878128e-12

@BilinearForm
def laplace(u, v, w):
    """Bilinear form for Laplace equation with variable permittivity."""
    # If dielectric is anisotropic, we could use a tensor, but for z-cut LN 
    # the transverse permittivity is isotropic eps11 = eps22 = eps_perp.
    eps = w.eps_perp
    return eps * dot(grad(u), grad(v))

class ElectrostaticModel2D:
    def __init__(self, mesh: MeshTri, eps_perp_func: Callable[[np.ndarray], np.ndarray]):
        """Initialize with a 2D mesh and a function returning eps_perp at given coordinates.
        
        Parameters
        ----------
        mesh : skfem.MeshTri
            The 2D triangular mesh.
        eps_perp_func : callable
            A function taking mesh.p (coordinates) and returning relative permittivity.
        """
        self.mesh = mesh
        self.element = ElementTriP1()
        self.basis = Basis(self.mesh, self.element)
        
        # Evaluate permittivity at mesh nodes
        p = self.mesh.p
        eps_vals = eps_perp_func(p)
        
        # Assemble the stiffness matrix
        self.A = laplace.assemble(self.basis, eps_perp=eps_vals)
        
    def solve_potential(self, boundary_conditions: dict) -> np.ndarray:
        """Solve Laplace equation given Dirichlet boundary conditions.
        
        Parameters
        ----------
        boundary_conditions : dict
            Dictionary mapping boundary names (defined in mesh.boundaries) to voltage values.
            
        Returns
        -------
        u : np.ndarray
            Nodal potential values.
        """
        D = {}
        for bnd_name, val in boundary_conditions.items():
            if bnd_name in self.mesh.boundaries:
                dofs = self.basis.get_dofs(bnd_name)
                D.update({dof: val for dof in dofs.all()})
        
        if not D:
            raise ValueError("No matching boundaries found in mesh for the given boundary conditions.")
            
        dofs_idx = np.array(list(D.keys()), dtype=int)
        dofs_vals = np.array(list(D.values()), dtype=float)
        
        # Free DOFs are all nodes minus Dirichlet nodes
        I = self.basis.complement_dofs(dofs_idx)
        
        # Right hand side is zero (Laplace equation)
        b = np.zeros(self.basis.N)
        x = np.zeros(self.basis.N)
        x[dofs_idx] = dofs_vals
        
        x_solved = solve(*condense(self.A, b, x=x, D=dofs_idx))
        return x_solved
    
    def compute_charge_per_length(self, u: np.ndarray, boundary_name: str) -> float:
        """Compute the total charge per unit length on a boundary.
        
        Q' = - epsilon_0 * int_S (epsilon_r * grad(V) . n) dS
        Using the variational formulation, Q' corresponds to the reaction force at the boundary.
        """
        # Reaction forces vector
        F = self.A @ u
        
        dofs = self.basis.get_dofs(boundary_name).all()
        # Sum of reaction forces on the boundary gives total flux * epsilon_0
        q_prime = np.sum(F[dofs]) * EPSILON_0
        return q_prime
    
    def extract_capacitance_matrix_per_length(self, electrode_names: List[str]) -> np.ndarray:
        """Extract the node capacitance matrix C'_node (F/m).
        
        For N electrodes, C_ij is extracted by setting V_j = 1 and V_k = 0 (k != j),
        and computing charge Q_i on electrode i.
        
        Handles shared corner DOFs: if a node belongs to two electrode boundaries,
        its reaction force is split equally between them.
        """
        N = len(electrode_names)
        
        # Pre-compute DOF sets and a weight vector for shared nodes
        dof_sets = []
        for name in electrode_names:
            dof_sets.append(set(self.basis.get_dofs(name).all()))
        
        # Build per-DOF weight: 1.0 for nodes in one boundary, 0.5 for nodes in two
        dof_weight = {}  # dof_index -> {electrode_index: weight}
        for i, dofs_i in enumerate(dof_sets):
            for d in dofs_i:
                if d not in dof_weight:
                    dof_weight[d] = {}
                dof_weight[d][i] = 1.0  # will be adjusted below
        
        for d, owners in dof_weight.items():
            n_owners = len(owners)
            if n_owners > 1:
                for i in owners:
                    owners[i] = 1.0 / n_owners
        
        C_node_prime = np.zeros((N, N))
        
        for j, driver in enumerate(electrode_names):
            bc = {name: (1.0 if name == driver else 0.0) for name in electrode_names}
            u = self.solve_potential(bc)
            F = self.A @ u  # reaction forces
            
            for i, sensor in enumerate(electrode_names):
                q = 0.0
                for d in dof_sets[i]:
                    w = dof_weight[d].get(i, 1.0)
                    q += F[d] * w
                C_node_prime[i, j] = q * EPSILON_0
                
        return C_node_prime
    
    def extract_electric_field(self, u: np.ndarray, points: np.ndarray) -> np.ndarray:
        """Evaluate electric field E = -grad(V) at given points.
        
        Parameters
        ----------
        u : np.ndarray
            Nodal potentials.
        points : np.ndarray (2, M)
            Coordinates to evaluate the field.
            
        Returns
        -------
        E : np.ndarray (2, M)
            Electric field vectors (Ex, Ey) at the given points.
        """
        from matplotlib.tri import Triangulation, LinearTriInterpolator
        tri = Triangulation(self.mesh.p[0], self.mesh.p[1], self.mesh.t.T)
        interp = LinearTriInterpolator(tri, u)
        grad_x, grad_y = interp.gradient(points[0], points[1])
        
        # Where points are outside the mesh, gradient returns masked arrays or NaNs
        if np.ma.is_masked(grad_x):
            grad_x = grad_x.filled(0.0)
            grad_y = grad_y.filled(0.0)
        else:
            grad_x = np.nan_to_num(grad_x, nan=0.0)
            grad_y = np.nan_to_num(grad_y, nan=0.0)
            
        return np.vstack((-grad_x, -grad_y))


def build_parallel_plate_mesh(width: float, gap: float, res: int = 20) -> MeshTri:
    """Build a simple parallel plate mesh for benchmark."""
    # Domain is a bit larger than the plates to reduce fringing
    # Actually, for a strict parallel plate benchmark, we can make the domain exactly width x gap
    # and set side boundaries to periodic or natural to get exact C = eps A / d.
    # To keep it simple, we'll just make the domain exactly width x gap.
    # Left and right boundaries will have zero Neumann (grad V . n = 0), so field is uniform.
    x = np.linspace(-width/2, width/2, res)
    y = np.linspace(-gap/2, gap/2, res)
    
    mesh = MeshTri.init_tensor(x, y)
    
    # Tag boundaries
    mesh = mesh.with_boundaries({
        'top': lambda p: p[1] == gap/2,
        'bottom': lambda p: p[1] == -gap/2,
    })
    return mesh


def build_four_electrode_mesh(crystal_size: float, gap: float, electrode_width: float, res: int = 30) -> MeshTri:
    """Build a four-electrode square cross-section mesh with partial electrodes.
    
    crystal_size: Size of the square crystal.
    gap: Distance between opposite electrodes.
    electrode_width: Width of the electrodes on the faces.
    
    NOTE: For full-face electrodes, use build_full_face_electrode_mesh instead.
    """
    # Simple tensor mesh for the whole crystal
    x = np.linspace(-crystal_size/2, crystal_size/2, res)
    y = np.linspace(-crystal_size/2, crystal_size/2, res)
    mesh = MeshTri.init_tensor(x, y)
    
    # Tag boundaries for the 4 electrodes
    # Assuming electrodes are centered on the faces
    hw = electrode_width / 2
    mesh = mesh.with_boundaries({
        '+x': lambda p: (p[0] == crystal_size/2) & (np.abs(p[1]) <= hw),
        '-x': lambda p: (p[0] == -crystal_size/2) & (np.abs(p[1]) <= hw),
        '+y': lambda p: (p[1] == crystal_size/2) & (np.abs(p[0]) <= hw),
        '-y': lambda p: (p[1] == -crystal_size/2) & (np.abs(p[0]) <= hw),
    })
    return mesh


def build_full_face_electrode_mesh(crystal_size: float, res: int = 50) -> MeshTri:
    """Build a four-electrode square cross-section mesh with full-face electrodes.
    
    Each electrode covers the entire corresponding side of the square.
    Corner nodes appear in two adjacent boundary groups simultaneously.
    This matches a z-cut crystal bar where full 3 mm x 30 mm face electrodes
    are applied to the left, right, top, and bottom faces.
    
    Parameters
    ----------
    crystal_size : float
        Side length of the square cross-section (m).  For a 3 mm crystal, pass 3e-3.
    res : int
        Number of grid points along each axis.
    
    Returns
    -------
    mesh : MeshTri
        Tagged with boundaries '+x', '-x', '+y', '-y'.
        Corner nodes belong to two boundaries each.
    """
    half = crystal_size / 2
    x = np.linspace(-half, half, res)
    y = np.linspace(-half, half, res)
    mesh = MeshTri.init_tensor(x, y)
    
    tol = crystal_size / (res - 1) * 0.1
    # Include corner nodes in every matching side
    mesh = mesh.with_boundaries({
        '+x': lambda p: np.abs(p[0] - half) < tol,
        '-x': lambda p: np.abs(p[0] + half) < tol,
        '+y': lambda p: np.abs(p[1] - half) < tol,
        '-y': lambda p: np.abs(p[1] + half) < tol,
    })
    return mesh


def gaussian_beam_weights(points: np.ndarray, w0: float) -> np.ndarray:
    """Compute Gaussian beam intensity weights for a set of points.
    
    I(x,y) = exp(-2 * (x^2 + y^2) / w0^2)
    """
    r2 = points[0]**2 + points[1]**2
    w = np.exp(-2 * r2 / w0**2)
    return w / np.sum(w)


def extract_field_per_volt_matrix(
    model: ElectrostaticModel2D, 
    electrode_names: List[str], 
    points: np.ndarray, 
    w0: float
) -> FieldPerVoltMatrix:
    """Extract the 2x2 field-per-volt matrix A and RMS nonuniformity."""
    weights = gaussian_beam_weights(points, w0)
    
    A = np.zeros((2, 2))
    rms_nonuniformity = np.zeros((2, 2))
    
    # We want A such that [Ex, Ey] = A @ [Vx, Vy] (differential drives)
    # Ex from Vx drive
    bc_x = {'+x': 0.5, '-x': -0.5, '+y': 0.0, '-y': 0.0}
    ux = model.solve_potential(bc_x)
    Ex_vec = model.extract_electric_field(ux, points)
    
    # Ey from Vy drive
    bc_y = {'+x': 0.0, '-x': 0.0, '+y': 0.5, '-y': -0.5}
    uy = model.solve_potential(bc_y)
    Ey_vec = model.extract_electric_field(uy, points)
    
    # A[0,0] = <Ex> from Vx
    A[0, 0] = np.sum(Ex_vec[0] * weights)
    # A[1,0] = <Ey> from Vx
    A[1, 0] = np.sum(Ex_vec[1] * weights)
    
    # A[0,1] = <Ex> from Vy
    A[0, 1] = np.sum(Ey_vec[0] * weights)
    # A[1,1] = <Ey> from Vy
    A[1, 1] = np.sum(Ey_vec[1] * weights)
    
    # RMS nonuniformity
    rms_nonuniformity[0, 0] = np.sqrt(np.sum((Ex_vec[0] - A[0, 0])**2 * weights))
    rms_nonuniformity[1, 0] = np.sqrt(np.sum((Ex_vec[1] - A[1, 0])**2 * weights))
    rms_nonuniformity[0, 1] = np.sqrt(np.sum((Ey_vec[0] - A[0, 1])**2 * weights))
    rms_nonuniformity[1, 1] = np.sqrt(np.sum((Ey_vec[1] - A[1, 1])**2 * weights))
    
    return FieldPerVoltMatrix(A=A, rms_nonuniformity=rms_nonuniformity)
