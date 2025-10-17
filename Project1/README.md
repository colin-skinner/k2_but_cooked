Useful commands

```
  julia --project=. src/run.jl data/small.csv 100 1000
  julia --project=. src/run.jl data/medium.csv 10 1000
  julia --project=. src/run.jl data/large.csv 1 300
```


```
julia --project=. src/run.jl data/..
```

```
julia --project=. -e "using Pkg; Pkg.test()"
```


```
julia --project=. -e "using Pkg; Pkg.test(test_args=[\"Score\"])"
```

4037.451155424665	-107581.2553374516	-479214.8037299974