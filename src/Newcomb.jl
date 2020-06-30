# ==================================================================================================================================
# Package: https://github.com/mkretlow/Newcomb.jl, Mike Kretlow [astrodynamics.de], Start 2018, MIT License
# ==================================================================================================================================

module Newcomb

using Libdl

export sunXYZ


# Load shared library

const lib = find_library(["libnewcomb"],[dirname(@__FILE__), joinpath(dirname(@__FILE__), "..", "deps"), "deps",])

if isempty(lib)
	error("Could not find shared library")
end

@info("Using shared library $lib")



"""*********************************************************************************************************************************
    pos = sunXYZ(yd::Float64, equin::Float64 = 2000.0)

Computes equatorial cartesian Solar coordinates using Newcomb's full theory of the Sun.
*********************************************************************************************************************************"""
function sunXYZ(yd::Float64, equin::Float64 = 2000.0)

   pos = zeros(3)

   ccall((:sonne_,lib),Nothing,(Ref{Float64},Ref{Float64},Ref{Float64}), yd,equin,pos)

   return pos

end

end  # module


#==================================================================================================================================#
