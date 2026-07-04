import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple

@dataclass(frozen=True)
class PockelsTensor:
    """Electro-optic tensor (Pockels tensor) r_ijk.
    
    Can be stored as a 6x3 contracted matrix or 3x3x3 full tensor.
    Default uses the 6x3 contracted Voigt notation commonly used.
    """
    contracted: np.ndarray  # 6x3 array, r_ij in m/V
    
    @property
    def full(self) -> np.ndarray:
        """Returns the full 3x3x3 tensor r_ijk.
        
        Mapping from Voigt (contracted) indices to full:
        1 -> 11, 2 -> 22, 3 -> 33, 4 -> 23, 5 -> 13, 6 -> 12
        Note: Python uses 0-based indexing (0..2).
        """
        r = np.zeros((3, 3, 3))
        # 1: 00
        r[0, 0, :] = self.contracted[0, :]
        # 2: 11
        r[1, 1, :] = self.contracted[1, :]
        # 3: 22
        r[2, 2, :] = self.contracted[2, :]
        # 4: 12, 21
        r[1, 2, :] = self.contracted[3, :]
        r[2, 1, :] = self.contracted[3, :]
        # 5: 02, 20
        r[0, 2, :] = self.contracted[4, :]
        r[2, 0, :] = self.contracted[4, :]
        # 6: 01, 10
        r[0, 1, :] = self.contracted[5, :]
        r[1, 0, :] = self.contracted[5, :]
        return r

@dataclass(frozen=True)
class DielectricTensor:
    """Relative permittivity tensor (3x3)."""
    tensor: np.ndarray  # 3x3 array

    @classmethod
    def from_diagonal(cls, eps11: float, eps22: float, eps33: float):
        return cls(np.diag([eps11, eps22, eps33]))

@dataclass(frozen=True)
class CapacitanceMatrix:
    """Maxwell node capacitance matrix (F) and differential capacitance matrix (F)."""
    C_node: np.ndarray  # NxN matrix
    T: np.ndarray       # NxM transformation matrix, e.g. 4x2
    
    @property
    def C_diff(self) -> np.ndarray:
        """Differential capacitance matrix C_diff = T.T @ C_node @ T"""
        return self.T.T @ self.C_node @ self.T
    
    def verify_symmetry(self, atol=1e-12) -> bool:
        """Verify the node capacitance matrix is symmetric."""
        return np.allclose(self.C_node, self.C_node.T, atol=atol)
    
    def electrostatic_energy(self, V_node: np.ndarray) -> float:
        """Compute positive electrostatic energy U = 0.5 * V.T @ C_node @ V"""
        V = np.asarray(V_node)
        return float(0.5 * V.T @ self.C_node @ V)

@dataclass(frozen=True)
class FieldPerVoltMatrix:
    """Field-per-volt matrix A (e.g. 2x2 for [Ex, Ey] from [Vx, Vy])."""
    A: np.ndarray
    rms_nonuniformity: Optional[np.ndarray] = None
