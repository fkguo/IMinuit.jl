# IMinuit

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fkguo.github.io/IMinuit.jl/stable) -->
<!-- [![Build Status](https://github.com/fkguo/IMinuit.jl/workflows/CI/badge.svg)](https://github.com/fkguo/IMinuit.jl/actions) -->

Julia wrapper of the Python package [`iminuit`](https://github.com/scikit-hep/iminuit), which is the interface to the C++ MINUIT2, widely used in fitting in the high-energy physics community. 
Supported `iminuit` versions: 2.0.0 or later

The `minuit` object in `iminuit` is defined as an `AbstractFit`:
if using array parameters, then `ArrayFit`;
if using individual parameters, then `Fit`.

Install by `]add https://github.com/fkguo/IMinuit.jl`
> Remember to point your `pycall` to the python enviroment which has `iminuit` installed, and then build your `pycall`, for more information, go check [PyCall](https://github.com/JuliaPy/PyCall.jl)

For functions defined, click [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fkguo.github.io/IMinuit.jl/dev)

For interactive examples, click
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fkguo/IMinuit.jl/master?urlpath=lab%2Ftree%2Fdocs%2Fexample.ipynb)
