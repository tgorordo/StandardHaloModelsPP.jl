mod docs

default:
  just --list

format:
  julia --color=yes --project=test -e 'using JuliaFormatter; format(["src", "test"])'

test:
  julia --color=yes --project -e "using Pkg; Pkg.test()"

aqua:
  julia --color=yes --project=test test/runaqua.jl

bench:
  julia --color=yes --project=test -e 'include("test/runbench.jl"); runbench()'
