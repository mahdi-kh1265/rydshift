# Electro-Optic Model

## Unperturbed index ellipsoid

For LiNbO3 in crystal axes:

```math
\eta_o x^2 + \eta_o y^2 + \eta_e z^2 = 1
```

or matrix form:

```math
\boldsymbol{\eta}^{(0)} = \begin{pmatrix}
\eta_o & 0 & 0 \\
0 & \eta_o & 0 \\
0 & 0 & \eta_e
\end{pmatrix}
```

where:

```math
\eta_o = \frac{1}{n_o^2}, \quad \eta_e = \frac{1}{n_e^2}
```

For z-propagation, the transverse block is:

```math
\boldsymbol{\eta}_\perp^{(0)} = \begin{pmatrix}\eta_o & 0 \\ 0 & \eta_o\end{pmatrix}
```

so the unperturbed `x-y` slice is a circle. No static retardance.

## EO perturbation for z-propagation

Using the common LiNbO3 point-group `3m` convention used in earlier derivation, the transverse EO block for propagation along `z` is:

```math
\boldsymbol{\eta}_{\perp} =
\begin{pmatrix}
\eta_o + r_{13}E_z - r_{22}E_y & -r_{22}E_x \\
-r_{22}E_x & \eta_o + r_{13}E_z + r_{22}E_y
\end{pmatrix}
```

The common term `r13 Ez` shifts both transverse eigenvalues equally and does not create transverse retardance by itself.

## Eigenvalue splitting

For a symmetric transverse matrix:

```math
\boldsymbol{\eta}_\perp = \begin{pmatrix}\eta_{xx} & \eta_{xy} \\ \eta_{xy} & \eta_{yy}\end{pmatrix}
```

Eigenvalues are:

```math
\lambda_\pm = \frac{\eta_{xx}+\eta_{yy}}{2}
\pm \sqrt{\left(\frac{\eta_{xx}-\eta_{yy}}{2}\right)^2 + \eta_{xy}^2}
```

and:

```math
n_\pm = \frac{1}{\sqrt{\lambda_\pm}}
```

For the EO block above:

```math
\lambda_\pm = \eta_o + r_{13}E_z \pm |r_{22}|\sqrt{E_x^2 + E_y^2}
```

For small perturbations:

```math
\Delta n \approx n_o^3 |r_{22}|\sqrt{E_x^2 + E_y^2}
```

## Retardance and half-wave field

```math
\Gamma = \frac{2\pi}{\lambda_0} L \Delta n
```

Thus:

```math
\Gamma \approx \frac{2\pi}{\lambda_0} L n_o^3 |r_{22}| E_\perp
```

Half-wave condition:

```math
\Gamma = \pi
```

so:

```math
E_\pi = \frac{\lambda_0}{2 L n_o^3 |r_{22}|}
```

If a simulated electrode geometry gives field-per-differential-volt `alpha = E_beam / Vdiff`, then:

```math
V_\pi = \frac{E_\pi}{\alpha}
```

For ideal plates separated by `d`, `alpha ≈ 1/d`.

## Eigenaxis angle

The instantaneous transverse principal-axis angle satisfies:

```math
\tan(2\theta) = \frac{2\eta_{xy}}{\eta_{xx}-\eta_{yy}}
```

Using the EO terms:

```math
\tan(2\theta) \sim \frac{E_x}{E_y}
```

up to sign convention. Quadrature `E_x` and `E_y` rotate the ellipse/eigenaxes.

## Rotating drive target

```math
E_x(t)=E_0\cos\phi(t), \quad E_y(t)=E_0\sin\phi(t)
```

Then `sqrt(Ex^2+Ey^2)=E0` gives constant retardance, while the eigenaxis angle rotates. For near half-wave retardance, circular input can be converted with geometric phase. Frequency shift comes from time derivative of the geometric phase.

