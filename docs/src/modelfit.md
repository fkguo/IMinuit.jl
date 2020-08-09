# chisq and model_fit


```
chisq(dist::Function, data, par; fitrange = ())
```
defines the χ² function: `fun` the function to be fitted to the data given by `data`.
The parameters are collected into `par`, given as an array or a tuple.
* `dist` should be defined as `dist(x, par)`.
* `data` can be either given as the `Data` type, or of the form `(xdata, ydata [, err])`.
If no `err` is given explicitly, the errors are assumed to be 1 for all data points.
* `fitrange`: default to the whole data set; may given as, e.g., `2:10`,
which means only fitting to the 2nd to the 10th data points.


If there is only one distribution to be fitted, then one can simply use:

```
model_fit(model::Function, data::Data, start_values; kws...)
@model_fit model data start_values kws...
```
the returning type is `ArrayFit`, which can be passed
to `migrad`, `minos` etc. to perform the fit and error analysis.
 `model` is the function to be fitted to `data`; it should be of the form
`model(x, params)` with `params` given either as an array or a tuple.

