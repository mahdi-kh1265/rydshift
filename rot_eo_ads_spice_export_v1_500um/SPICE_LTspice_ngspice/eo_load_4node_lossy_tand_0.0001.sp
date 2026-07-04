* EO Load Subcircuit: eo_load_4node_lossy_tand_0.0001.sp
* Loss tangent: 0.0001
* Original Maxwell C_node matrix (F):
* 3.0816e-11  -2.9435e-12  -1.3936e-11  -1.3936e-11
* -2.9435e-12  3.0816e-11  -1.3936e-11  -1.3936e-11
* -1.3936e-11  -1.3936e-11  3.0116e-11  -2.2430e-12
* -1.3936e-11  -1.3936e-11  -2.2430e-12  3.0116e-11

.SUBCKT EO_LOAD_4NODE p_x m_x p_y m_y
C1 p_x m_x 2.943527e-12
R1 p_x m_x 5.406946e+06 * evaluated at 100 MHz for tand=0.0001
C2 p_x p_y 1.393641e-11
R2 p_x p_y 1.142008e+06 * evaluated at 100 MHz for tand=0.0001
C3 p_x m_y 1.393641e-11
R3 p_x m_y 1.142008e+06 * evaluated at 100 MHz for tand=0.0001
C4 m_x p_y 1.393641e-11
R4 m_x p_y 1.142008e+06 * evaluated at 100 MHz for tand=0.0001
C5 m_x m_y 1.393641e-11
R5 m_x m_y 1.142008e+06 * evaluated at 100 MHz for tand=0.0001
C6 p_y m_y 2.243040e-12
R6 p_y m_y 7.095502e+06 * evaluated at 100 MHz for tand=0.0001
.ENDS
