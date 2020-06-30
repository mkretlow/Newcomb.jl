#
# Build shared lib from Fortran source
#

function getcompiler()

    haskey(ENV, "FC") && return ENV["FC"]

    if success(`which ifort`)
        return "ifort"

    elseif success(`which gfortran`)
        return "gfortran"

    else
        error("No compatible Fortran compiler found.")
    end

end


function unixbuild(compiler, path, ext)

    if compiler == "ifort"
        run(`ifort -$(Sys.isapple() ? "dynamiclib" : "shared") -O3 -xHost -ipo -fpic -o $path/libnewcomb.$ext newcomb_jl.for`)

    elseif compiler == "gfortran"
        run(`gfortran -O2 -cpp -fPIC -shared -o $path/libnewcomb.$ext newcomb_jl.for`)
    end

    println("Unix build finished: $path/libnewcomb.$ext")

end


function windowsbuild(path)
    error("Uups. Windowsbuild not implemented. Sorry.")
end


# Main starts here:

compiler = Sys.isunix() ? getcompiler() : ""

path = dirname(@__FILE__)

ext = Sys.isapple() ? "dylib" : "so"

build() = Sys.isunix() ? unixbuild(compiler, path, ext) : windowsbuild(path)

build()
