---
title: 'AtmosphericModels.jl: Wind profiles, air density and turbulent wind fields for wind energy simulations'
tags:
  - Julia
  - atmospheric boundary layer
  - wind shear
  - wind profile
  - turbulence
  - Mann model
  - airborne wind energy
authors:
  - name: Uwe Fechner
    orcid: 0009-0008-2532-9458
    affiliation: 1
affiliations:
  - name: Delft University of Technology, The Netherlands
    index: 1
date: 11 September 2026
bibliography: paper.bib
---

# Summary

Any simulation of a wind energy system needs a description of the wind it operates in. For
a conventional wind turbine that description covers a rotor disc of a hundred metres or
so, centred on a fixed hub height. An airborne wind energy (AWE) system is different: a
tethered kite or aircraft sweeps through the atmosphere on a trajectory that spans several
hundred metres in height and moves with the wind, so the simulation has to know the mean
wind speed at every height the kite and the tether reach, the air density at those heights,
and, for realistic control and load studies, the turbulent wind vector at any point of the
trajectory at any time.

`AtmosphericModels.jl` provides exactly these three ingredients as a small Julia
[@Bezanson2017] package. It returns the air density as a function of height, the mean wind
speed as a function of height under a selectable wind profile law, and a three-dimensional
turbulent wind vector at any position and time, generated once with the spectral tensor
model of @Mann1994 [-@Mann1998] and read back in a few tens of nanoseconds per call. The
package was written as the atmosphere component of a family of open-source tools for the
simulation and control of kite power systems [@Fechner2015; @Fechner2016; @KiteModels],
but it depends on none of them and can be used on its own.

![Three of the built-in wind profile laws for a near-shore site: power law ($v_\mathrm{w,exp}$),
logarithmic law ($v_\mathrm{w,log}$) and their combination ($v_\mathrm{w,fit}$), fitted to
wind speeds measured at three heights [@Fechner2015].\label{fig:profile}](wind_profile.png){ width=65% }

# Statement of need

The two standard sources of synthetic turbulence for wind energy, TurbSim [@Jonkman2009]
and the Mann turbulence generator, are stand-alone programs that write a wind field file
which a separate aeroelastic code then reads. The grids they are designed for are the size
of a rotor disc, and the file formats and reading routines are tied to the wind turbine
codes they were written for. AWE researchers who want a turbulent inflow for a kite
simulation written in Julia therefore end up either coupling their model to one of those
tool chains, or writing their own field generator and lookup, which is what happened in
the author's earlier work [@Fechner2016]. Likewise, the height-dependent mean wind and the
air density are simple enough that every simulation code re-implements them, usually with
one hard-coded profile law and no way to fit a profile to measurements.

`AtmosphericModels.jl` packages these pieces so that they can be re-used and tested. It is
an in-process library rather than a file-based tool chain: a model is constructed from a
settings file, and `get_wind(am, x, y, z, t)` returns the wind vector at that point. The
wind field is a box of 8.1 km along the wind by 200 m across by 1 km in height at 2 m
resolution by default, which is long enough that the turbulence seen by the kite repeats
only every 13.5 minutes at 10 m/s. The field is advected past the kite with the mean wind
according to Taylor's frozen-turbulence hypothesis [@Taylor1938], so the position of the
kite is shifted along the wind direction by the product of time and mean wind speed before
the lookup, and the field itself is never recomputed during a simulation. Because the
inverse FFT makes the field periodic, the lookup wraps around in the two horizontal
directions and a simulation can run for hours without leaving the box. The wind direction
is a parameter of the lookup: the query position is rotated into the frame of the stored
field, so its long axis is always aligned with the mean wind. An asymmetric, rectangular
field, long along the wind and narrow across it, is therefore sufficient to obtain the same
long repetition period for any wind direction, without storing a square field that would be
forty times larger. The turbulent component is scaled
at read time by the turbulence intensity requested in the settings, so one stored field per
ground wind speed serves every intensity.

For the mean wind, the package implements the power law and the logarithmic law
[@Stull2000; @Burton2001], and the combination of both that was used in @Fechner2015 to fit
a profile to wind speeds measured at three heights (\autoref{fig:profile}). Three further
laws fit a profile directly to a list of measured height/speed pairs from the settings, e.g.
from a met mast or lidar: a logarithmic and a power-law fit by linear least squares, and a
low-level-jet profile, a power-law background plus a Gaussian bump, fitted by a
Levenberg-Marquardt solve. The jet profile
matters for AWE because the height range of 200 to 500 m that kites exploit is exactly
where nocturnal low-level jets appear, and the wind resource at those heights differs
markedly from what a power law extrapolated from 10 m would predict [@Bechtle2019]. The
air density follows an exponential decay with height with a scale height of 8550 m,
corrected for the reference temperature at ground level.

The default settings shipped with the package describe two real sites in the Netherlands.
For the inland site Cabauw, the shear exponent and the turbulence intensity at three
reference wind speeds were calibrated in @Fechner2016 against one year of measurements from
the KNMI meteorological tower [@KNMI2011]; the turbulence intensity of the IEC normal
turbulence model [@IEC61400-1] is scaled by these calibration factors. For the near-shore
site Maasvlakte, the profile in \autoref{fig:profile} is used.

# Implementation

The turbulence generator is a Julia port of a MATLAB implementation of the Mann model by
René Bos. The spectral tensor is evaluated on the wave-number grid, multiplied by a white
noise vector with random phase, and transformed to physical space by a three-dimensional
inverse FFT [@Frigo2005]. Generation of the default field takes about 30 s and is
deterministic, because a fixed-seed stable random number generator is used; the field is
written once to a scratch space as a NumPy `.npz` archive and loaded on subsequent runs. Every setting that influences the field is hashed into the file name, so changing a
parameter causes a new field to be generated rather than a stale one to be loaded.

The read path is designed to sit inside the inner loop of a simulation. The mean wind
profile is evaluated by a method that dispatches on the profile law at compile time, the
turbulence lookup reads either the nearest grid point or a trilinear interpolation of the
eight surrounding points, and a vectorised method evaluates the wind at many positions --
for example, all particles of a discretised tether -- while computing the position
independent quantities once. On a laptop CPU the air density costs about 3 ns, the
nearest-point turbulence lookup about 30 ns and the interpolated one about 50 ns per
position. The `.npz` field files are shared between all packages that use
`AtmosphericModels.jl` on the same machine via a Julia scratch space, so that a downstream
project does not accumulate its own copies of files that are of the order of a gigabyte.

# Functionality

The package exports the `AtmosphericModel` type, constructed from a `Settings` object of
the companion package `KiteUtils.jl`; `calc_rho` for the air density; `calc_wind_factor`
for the ratio of the wind speed at a height to the ground wind speed under any of the seven
profile laws (`CONSTANT`, `EXP`, `LOG`, `EXPLOG`, `CUSTOM_LOG`, `CUSTOM_EXP`, `CUSTOM_JET`);
`get_wind` and `get_wind!` for the turbulent wind vector at a point or at a vector of
points; `calc_turbulent_wind`, which returns the wind at the kite and at half its height
for the tether; and `new_windfield`/`new_windfields` to (re)generate the stored fields.
Examples reachable from an interactive menu plot the profile laws, the custom fits, the
turbulent wind as a function of time and slices of the wind field. The test suite checks
the profile laws and the air density against closed-form values, the custom fits against
synthetic profiles, and the wind field generation, storage and lookup end to end.
Documentation with annotated settings for both sites is published online, and the package
is archived on Zenodo [@AtmosphericModels].

# AI usage disclosure

Generative AI tools were used in the development of this software and in the preparation
of this paper, as follows.

*Software.* The original Python implementation of the atmospheric model and its port to
Julia, including the Mann turbulence generator, the profile laws and the air density model,
were written by the author without AI assistance. From 2026 onwards, Claude Code
(Anthropic) was used as a coding assistant for later additions and refactorings, among them
the custom profile fits, the trilinear interpolation and vectorised lookup, the move of the
wind field files to a scratch space with a settings digest in the file name, and parts of
the accompanying tests and documentation. Every AI-assisted change was reviewed by the
author and verified by the test suite.

*Paper.* This paper was drafted with Claude Code from the package's source code,
documentation and change log, and then reviewed and edited by the author, who takes full
responsibility for its content.

# Acknowledgements

The Mann model implementation is based on a MATLAB module written by René Bos. The
calibration of the Cabauw scenario used data from the CESAR observatory of the Royal
Netherlands Meteorological Institute. The author thanks the contributors to the package
and the developers of the Julia packages it builds on, in particular FFTW.jl.

# References
