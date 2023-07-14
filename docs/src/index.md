```@meta
CurrentModule = IMinuit
```

# IMinuit.jl


Julia wrapper of the Python package [`iminuit`](https://github.com/scikit-hep/iminuit), which is the interface to the C++ MINUIT2, widely used in fitting in the high-energy physics community.

> **IMPORTANT**
> In the v2.0 or later, huge changes are introduced to the interface, in this package you can still use the old usage like passing `limit_`, `fix_` and the like to the `Minuit` function, but the code to implement these usage is over complex, if you have better implementation you can pull a request, or you can learn the new usage if you want to keep your code robust and reusable.
> Also keep in mind that the introduction in this doc is different from that pulled from the python package, this doc would not remind you of that over and over again.
    
    
```@index
```

```@autodocs
Modules = [IMinuit]
``` 
