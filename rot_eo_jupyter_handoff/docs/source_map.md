# Source Map

This file lists source material Claude/the developer should consult while implementing the model. It does not embed full copyrighted documents; it records source URLs and the specific reason each source matters.

## Electro-optic / high-voltage RF driver precedent

### Holger Müller, "Fast high-voltage amplifiers for driving electro-optic modulators"
- URL: https://pubs.aip.org/aip/rsi/article/76/8/084701/823516/Fast-high-voltage-amplifiers-for-driving-electro
- Also: https://arxiv.org/abs/physics/0506050
- Why it matters: direct precedent for high-voltage, high-speed linear amplifiers for EO modulators. Reported examples include 60-550 Vpp and 1.3-300 MHz bandwidth into capacitive or resistive loads. Use as conceptual justification that a custom broadband HV RF amplifier is a real path.

### Newport / New Focus, Practical Uses and Applications of Electro-Optic Modulators
- URL: https://www.newport.com/n/practical-uses-and-applications-of-electro-optic-modulators/
- Why it matters: discusses broadband bulk EOMs dominated by crystal capacitance in the DC-to-~100 MHz range; gives capacitance examples in the 10-30 pF range. This supports our 13 pF load assumption as realistic.

### MKS/New Focus 400X and 406X User Manual
- URL: https://api.p1.mks.com/medias/sys_master/images/images/h65/hcc/8797007839262/400X-and-406X-User-Manual-Rev-J.pdf
- Why it matters: manual for DC-250 MHz phase modulators. Discusses resonant and broadband EO drive architectures; useful comparison between resonant voltage enhancement and broadband drive.

### Conoptics EO modulation systems
- URL: https://www.conoptics.com/eo-modulation-systems/
- Why it matters: example of broadband EO system using RF amplifier and terminated RF/electrode structure; illustrates the broadband-vs-power tradeoff.

## LDMOS / RF PA parts and references

### NXP MRF300AN/BN datasheet
- URL: https://www.nxp.com/docs/en/data-sheet/MRF300AN.pdf
- Why it matters: 300 W CW, 1.8-250 MHz, 50 V wideband RF power LDMOS. Candidate serious one-axis PA device. Mirror AN/BN pinouts simplify push-pull/two-up layouts.

### NXP MRF300 product page
- URL: https://www.nxp.com/products/radio-frequency-rf/legacy-rf/legacy-rf-power/300-w-cw-over-1-8-250-mhz-50-v-wideband-rf-power-ldmos-transistor%3AMRF300AN
- Why it matters: current vendor product page and design-resource entry point.

### NXP MRF101AN/BN datasheet
- URL: https://www.nxp.com/docs/en/data-sheet/MRF101AN.pdf
- Why it matters: 100 W CW, 1.8-250 MHz, 50 V wideband RF LDMOS. Good lower-cost learning/prototype part.

### NXP MRF101 product page
- URL: https://www.nxp.com/products/radio-frequency-rf/legacy-rf/legacy-rf-power/100-w-cw-over-1-8-250-mhz-50-v-wideband-rf-power-ldmos-transistor%3AMRF101AN
- Why it matters: current vendor product page and reference-circuit entry point.

### NXP MRF101 reference circuits
- URL: https://www.nxp.com/products/radio-frequency-rf/legacy-rf/legacy-rf-power/mrf101an-reference-circuits%3AMRF101AN-REFCirc
- Why it matters: concrete layouts/reference designs across HF/VHF frequencies; useful for biasing, layout, and matching style.

### Ampleon BLF645 broadband amplifier reference
- URL: https://www.ampleon.com/documents/application-note/AN10953.pdf
- Why it matters: example of broadband LDMOS PA design over a wide frequency range. Even if we do not use this exact part, it is useful for topology/layout thinking.

## RF matching / network modeling

### Kerr, "Some Fundamental and Practical Limits on Broadband Matching to Capacitive Devices"
- URL: https://www.gb.nrao.edu/electronics/edir/edir295.pdf
- Why it matters: Bode-Fano-type matching limits for capacitive loads. Useful for explaining why broadband high-voltage passive matching cannot be magic.

### LibreTexts Fano-Bode limits
- URL: https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Electronics/Microwave_and_RF_Design_III_-_Networks_%28Steer%29/07%3A_Chapter_7/7.2%3A_Fano-Bode_Limits
- Why it matters: accessible explanation: more reactive energy stored in a load means narrower achievable match bandwidth.

### Mini-Circuits RF transformer notes
- URL: https://blog.minicircuits.com/application-note-on-transformers-an-20-002/
- Why it matters: practical transformer/balun voltage, current, and impedance-ratio relationships.

### scikit-rf docs
- URL: https://scikit-rf.readthedocs.io/en/latest/tutorials/Networks.html
- Why it matters: Python network modeling using S, Z, Y, and ABCD/network representations; needed for measured S-parameter integration and DPD model.

## Python / FEM / optimization tooling

### scikit-fem docs
- URL: https://scikit-fem.readthedocs.io/
- Why it matters: pure-Python finite-element assembly library, useful for 2D electrostatic simulation and capacitance extraction.

### SciPy least_squares docs
- URL: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
- Why it matters: optimization routine for fitting network parameters, component values, and model-to-measurement calibration.

### SimPEG / discretize
- URL: https://simpeg.xyz/
- URL: https://github.com/simpeg/discretize
- Why it matters: alternative finite-volume / inverse-problem ecosystem if scikit-fem becomes inconvenient for capacitance extraction or inverse geometry optimization.

## Topics to search more deeply later

- LiNbO3 Sellmeier equations and wavelength-dependent EO coefficients near 780 nm.
- Dielectric tensor / RF permittivity of LN at 6-100 MHz.
- Loss tangent of LN at RF, especially for expected field strengths.
- Breakdown field / safe RF electrode voltage for bulk LN and air/polyimide/epoxy fixtures.
- RF high-voltage capacitive dividers/pickups for hundreds of volts at 100 MHz.
- Push-pull LDMOS class-AB PA output-network design with non-50-ohm capacitive load.

