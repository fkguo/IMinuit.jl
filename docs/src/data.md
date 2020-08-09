# Useful functions


```
Data(x::T, y::T, err::T) where {T<:Vector{Real}}
Data(df::DataFrame)
```
* Fields: `x, y, err, ndata`. This defines a type for data with three columns:` x, y, err`;
* `ndata` is the number of data rows, automatically counted.
* Different `Data` sets can be concatenated as `vat(dat1, dat2, dat3)`.
* `getindex` defined for `Data`, e.g., `data[1:10]` gives the first 10 points in `data`.



```
@plt_data(data, kws...)
@plt_data!(data, kws...)
@plt_best(dist, fit, data, kws...)
@plt_best!(dist, fit, data, kws...)
```
convenient macros for plotting the data, and a comparison of the best-fit curve
with the data (needs `using Plots`). `kwds...` in `@plt_data` adjusts settings of the `scatter` plot, and that in `@plt_best` adjusts settings of the best-fit curve `plot`.
For `@plt_best` and `@plt_best!`, the ordering of `dist`, `fit`, and `data` does not matter, i.e., `@plt_best dist fit data`,  `@plt_best fit dist data`, `@plt_best fit data dist` and so on all work.