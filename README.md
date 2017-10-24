# GraphChallenge

[![Build Status](https://travis-ci.org/rohitvarkey/GraphChallenge.jl.svg?branch=master)](https://travis-ci.org/rohitvarkey/GraphChallenge.jl)

[![Coverage Status](https://coveralls.io/repos/rohitvarkey/GraphChallenge.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/rohitvarkey/GraphChallenge.jl?branch=master)

[![codecov.io](http://codecov.io/github/rohitvarkey/GraphChallenge.jl/coverage.svg?branch=master)](http://codecov.io/github/rohitvarkey/GraphChallenge.jl?branch=master)


### Usage

An experiment can be launched on the static graph with 50 nodes using an array
as the interblock edge matrix using

```julia
using GraphChallenge
static_partition_experiment(Array{Int64, 2}, 50)
```


