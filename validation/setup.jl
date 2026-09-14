using Pkg
Pkg.activate(@__DIR__)
Pkg.Registry.add("General")
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.instantiate()
