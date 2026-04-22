# SPDX-FileCopyrightText: 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

# build and display the html documentation locally
# you must have installed the package LiveServer in your global environment

using Pkg

if !("Documenter" ∈ keys(Pkg.project().dependencies))
    Pkg.activate(joinpath(@__DIR__, "..", "docs"))
end
using LiveServer; servedocs(launch_browser=true)
