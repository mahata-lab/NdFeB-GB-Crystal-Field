# Figure 5 contrast and evidence

The six archived orbital cubes, integration direction (c), a-b coordinate plane, one-electron normalization and all orbital radii are unchanged. The higher-contrast figure retains the existing panel layout and ion choices. No new simulation or ground-state density reconstruction is involved.

## Display mapping

The plotted quantity is the projected density P(a,b)=integral |psi(a,b,c)|^2 dc/N, in Å^-2. One fixed density range 0–7 Å^-2 is used for every ion and environment. Only the density-to-color coordinate changes from P/7 to u=sqrt(P/7). The same sequential white-to-blue palette and mapping apply to all six panels. There is no independent ion/panel rescaling, clipping of a density peak, or logarithmic floor. The colorbar labels show P rather than u.

The square-root transfer is monotonic and invertible: P=7u^2. For example, P=0.25 Å^-2 occupies 3.6% of a linear color range but 18.9% of the square-root range, making weak lobes visible. Zero remains zero. Color distances no longer represent equal density increments; this is why the nonlinear mapping is explicit in the figure, caption and supplement. These colors do not represent signed orbital amplitude or a newly calculated observable.

The quantitative evidence remains the full 3D cube integrals in Table S10: Ce 0.8374→0.6095 Å, Pr 0.6612→0.5849 Å, Dy 0.4715→0.4646 Å. These are individual eigenorbital radii, not total correlated-shell radii or a unique mechanism of the crystal-field crossover. The source numerical tables and other six figure topics were compared with the pre-change package and are unchanged.

## Literature and technical basis

- Crameri, Shephard and Heron, *The misuse of colour in science communication*, Nature Communications 11, 5444 (2020), [DOI](https://doi.org/10.1038/s41467-020-19160-7). This supports transparent, ordered color encoding and inclusion of a colorbar; it is not evidence for our physical orbital interpretation or a specific square-root rule.
- Rahat and Kraisler, *Plateaus in the Potentials of Density-Functional Theory: Analytical Derivation and Useful Approximations*, J. Chem. Theory Comput. 21, 3476–3492 (2025), [DOI](https://doi.org/10.1021/acs.jctc.4c01771). Figure 3 presents calculated FCI densities logarithmically to expose spatial decay. This is an example of explicitly transformed density display; it does not calculate the lanthanide systems here or prescribe our transfer function.
- [Matplotlib's primary documentation on power-law color normalization](https://matplotlib.org/stable/tutorials/colors/colormapnorms.html#power-law) gives gamma=0.5 as an example and explicitly cautions that power-law mapping needs justification and disclosure. We choose the square-root coordinate to reveal weak nonnegative density while retaining zero without a floor. The plot is generated directly from this mathematical mapping, independently of Matplotlib's implementation.

The main discussion cites Crameri et al.; the response explains the display choice. The supplementary Methods define u, P and the common upper limit. Stronger contrast is a presentation change, not additional scientific evidence for a model-dependent mechanism.

## Orbital labels, ticks and contours

Every panel title identifies the ion, Bulk/GB environment and archived file index MO27. The spatial axes have labeled ticks every 0.5 Å and minor ticks every 0.25 Å. Thin black contours use the same projected-density levels (0.25, 0.75 and 2.0 Å^-2) in all six panels. Bicubic interpolation refines only the contour display; all numerical cube integrals and radii are unchanged. The lines represent isodensity, not coordination boundaries.
