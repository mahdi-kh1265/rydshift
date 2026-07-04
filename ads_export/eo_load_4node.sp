* EO Load Subcircuit: eo_load_4node.sp
* Loss tangent: 0.0
* Original Maxwell C_node matrix (F):
* 3.0816e-11  -2.9435e-12  -1.3936e-11  -1.3936e-11
* -2.9435e-12  3.0816e-11  -1.3936e-11  -1.3936e-11
* -1.3936e-11  -1.3936e-11  3.0116e-11  -2.2430e-12
* -1.3936e-11  -1.3936e-11  -2.2430e-12  3.0116e-11

.SUBCKT EO_LOAD_4NODE p_x m_x p_y m_y
C1 p_x m_x 2.943527e-12
C2 p_x p_y 1.393641e-11
C3 p_x m_y 1.393641e-11
C4 m_x p_y 1.393641e-11
C5 m_x m_y 1.393641e-11
C6 p_y m_y 2.243040e-12
.ENDS
