# MaxwellSALT

[![Build Status](https://travis-ci.org/wsshin/MaxwellSALT.jl.svg?branch=master)](https://travis-ci.org/wsshin/MaxwellSALT.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/k5u98x6ud03f7ygh?svg=true)](https://ci.appveyor.com/project/wsshin/maxwellsalt-jl)
[![codecov](https://codecov.io/gh/wsshin/MaxwellSALT.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/wsshin/MaxwellSALT.jl)

## Installation
1. Install Julia 1.0.  The binary is available at https://julialang.org/downloads/.
1. Start the Julia REPL (similar to the Python Interpreter or MATLAB Command Line, where you can input Julia commands).  On LINUX this is done by executing `julia` in a terminal.
1. Install registered packages.  To do this, type `]` in REPL to enter the package mode (the prompt changes from `julia>` to `pkg>`), and execute the followings line-by-line:
    - `add Arpack`
    - `add PyPlot`
1. Install unregistered packages by executing the followings line-by-line
    - `add https://github.com/stevengj/GeometryPrimitives.jl`
    - `add https://github.com/wsshin/MaxwellFDM.jl`
    - `add https://github.com/wsshin/SALTBase.jl`
    - `add https://github.com/wsshin/MaxwellSALT.jl`
1. Press backspace in REPL to exit the package mode (the prompt changes from `pkg>` to `julia>`).
