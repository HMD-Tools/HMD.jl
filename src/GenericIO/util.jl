
asserting() = false # when set to true, this will enable all `@myassert`s
macro myassert(test)
    esc(:(if $(@__MODULE__).asserting()
      @assert($test)
    end))
end

function assume_element(mass::AbstractString)
    mass = parse(Float64, mass)
    assume_element(mass)
end

function assume_element(mass::Real)
    if mass > 239
        error("TRans-Uranium is not supported. Detected: $mass Da")
    elseif mass <= 0
        error("Atomic mass must be positive. Detected: $mass Da")
    end
    for e in elements
        if ≈(ustrip(e.atomic_mass), mass, atol=1e-2)
            return Symbol(e.symbol)
        end
    end

    error("Unknown mass of atoms $mass detected. Atomic mass must be in Dalton and have 1e-2 presicion or more.")
end

function getlines(filename)
    fp = open(filename, "fp")
    lines = readlines(fp)
    close(fp)

    lines
end

macro scan(format::AbstractString, str)
    fmt = split(format)
    if any(s ->  match(r"^%[fd]$", s) == nothing , fmt)
        error("""fmt must only contain "%f", "%d" and space  """)
    end

    ex = Expr(:tuple)
    ex.args = map( 1:length(fmt) ) do i
        if fmt[i] == "%f"
            :(parse(Float64, arr[$i]))
        else
            :(parse(Int64, arr[$i]))
        end
    end
    quote
        arr = split($(esc(str)))
        $(ex)
    end
end

# @scan "%i %i %f %f" line