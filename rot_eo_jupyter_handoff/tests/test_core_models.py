import numpy as np

from rot_eo_model.eo import half_wave_field_V_per_m, vpi_from_gap_V, transverse_eta_matrix, eigen_indices
from rot_eo_model.rf_lumped import cap_current_peak_A, reflected_capacitance_F


def test_vpi_order_of_magnitude():
    vpi = vpi_from_gap_V(3e-3)
    assert 300 < vpi < 800


def test_half_wave_field_positive():
    assert half_wave_field_V_per_m() > 0


def test_rf_current_13pf_100mhz_480v():
    I = cap_current_peak_A(100e6, 13e-12, 480)
    assert np.isclose(I, 3.92, rtol=0.05)


def test_reflected_capacitance():
    assert np.isclose(reflected_capacitance_F(13e-12, 2), 52e-12)


def test_eo_eigen_split_nonzero():
    from rot_eo_model.materials import linbo3_pockels_tensor
    eta = transverse_eta_matrix(1e5, 0.0, pockels=linbo3_pockels_tensor())
    n, _ = eigen_indices(eta)
    assert abs(n[1] - n[0]) > 0
