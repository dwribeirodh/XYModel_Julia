# This script will install all necessary packages for use of the repository
# XYModel_Julia.
# Any questions, email ribei040@umn.edu

using Pkg


pkg_names = [
                "Distributions";
                "Colors";
                "ProgressBars";
                "Plots";
                "FiniteDifferences";
                "LinearAlgebra";
                "SpecialFunctions";
                "DelimitedFiles";
                "QuadGK";
                "Dates"
        ]

for name in pkg_names
        Pkg.add(name)
end

check = keys(Pkg.installed())

for name in pkg_names
        !in(name, check) ? println("Package < "*name *" > failed to install") : print("")
end
