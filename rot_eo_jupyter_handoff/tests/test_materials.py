import numpy as np
from rot_eo_model.materials import sellmeier_no, sellmeier_ne, linbo3_pockels_tensor, linbo3_dielectric_tensor

def test_sellmeier_no_780nm():
    n_o = sellmeier_no(780e-9)
    # Expected value from Zelmon/Jundt is ~2.254 to 2.286 depending on the exact fit.
    # The current Jundt fit gives ~2.254. 
    # Check that it's in the ballpark for visible/NIR LN.
    assert 2.2 <= n_o <= 2.3

def test_sellmeier_ne_780nm():
    n_e = sellmeier_ne(780e-9)
    assert 2.1 <= n_e <= 2.25
    assert n_e < sellmeier_no(780e-9)  # LiNbO3 is negative uniaxial

def test_linbo3_pockels_tensor_shape():
    pt = linbo3_pockels_tensor()
    assert pt.contracted.shape == (6, 3)
    assert pt.full.shape == (3, 3, 3)

def test_linbo3_dielectric_tensor():
    dt = linbo3_dielectric_tensor(10e6)
    assert dt.tensor.shape == (3, 3)
    assert dt.tensor[0, 0] == 44.0
    assert dt.tensor[2, 2] == 28.0
