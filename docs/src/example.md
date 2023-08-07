# Example

```@example 1
using IMinuit # hide
fcn(x) = x[1]^2 + (x[2]-2)^2 + (x[3]-3.2)^4
pnames = ["a", "b", "c"]
m = Minuit(fcn, [1,0,1]; name=pnames, error=0.1*ones(3), fix_a=true, limit_b=(0, 50) )
migrad(m)
hesse(m)
```

Here the parameters are collected in an array `x` which is the argument of `fcn(x)`,
and `typeof(m) = ArrayFit`. In this case, one can use external code to compute the gradient as
`gradfun(par) = gradient(fcn, par)` (the exported `gradient` function is from `ForwardDiff`),
and include `grad = gradfun` as a keyword argument.

If `fcn` is defined as `fcn(a, b)`, then the starting values need to be set as
`Minuit(fcn, a = 1, b = 0)`.

The parameters from the previous fit can be used as input to a new fit:
```@example 1
m1 = Minuit(fcn, m, name=pnames, fix_a=false, limit_b=(-10,10))
migrad(m1)
```
In the above example, the renaming keyword `name=pnames` is needed since the parameters 
have been named in the previous fit `m`.

The asymmetric errors can be obtained by using `minos`:
```@example 1
minos(m1)
```

```@example 1
migrad(m1)
```
