#--------------------------------------------------------------------------------# Datum
"""
    Datum(name, a, f_inv)

A geodetic datum defined by an ellipsoid of revolution.

# Fields
- `name::String`: Human-readable name (e.g. `"WGS 84"`).
- `a::Float64`: Semi-major axis in meters (equatorial radius).
- `f_inv::Float64`: Inverse flattening (1/f). Use `Inf` for a perfect sphere.

### Examples
```julia
julia> WGS84
Datum("WGS 84", a=6.378137e6, 1/f=298.257223563)

julia> semi_major_axis(WGS84)
6.378137e6

julia> semi_minor_axis(WGS84)
6.356752314245179e6
```
"""
struct Datum
    name::String
    a::Float64
    f_inv::Float64
end

function Base.show(io::IO, d::Datum)
    print(io, "Datum(\"", d.name, "\", a=", d.a, ", 1/f=", d.f_inv, ")")
end

#--------------------------------------------------------------------------------# Derived quantities
"""
    semi_major_axis(d::Datum) -> Float64

Return the semi-major axis (equatorial radius) in meters.
"""
semi_major_axis(d::Datum) = d.a

"""
    flattening(d::Datum) -> Float64

Return the flattening `f = 1 / f_inv`.
"""
flattening(d::Datum) = 1.0 / d.f_inv

"""
    semi_minor_axis(d::Datum) -> Float64

Return the semi-minor axis (polar radius) in meters: `b = a * (1 - f)`.
"""
semi_minor_axis(d::Datum) = d.a * (1.0 - flattening(d))

"""
    eccentricity_sq(d::Datum) -> Float64

Return the first eccentricity squared: `e² = 2f - f²`.
"""
function eccentricity_sq(d::Datum)
    f = flattening(d)
    return 2f - f^2
end

"""
    eccentricity(d::Datum) -> Float64

Return the first eccentricity: `e = √(2f - f²)`.
"""
eccentricity(d::Datum) = sqrt(eccentricity_sq(d))

#--------------------------------------------------------------------------------# Common datums
"""
    WGS84

World Geodetic System 1984. The standard datum used by GPS.
"""
const WGS84 = Datum("WGS 84", 6_378_137.0, 298.257223563)

"""
    GRS80

Geodetic Reference System 1980. Used by NAD83 and ITRF.
"""
const GRS80 = Datum("GRS 80", 6_378_137.0, 298.257222101)

"""
    CLARKE1866

Clarke 1866 ellipsoid. Used by the original NAD27.
"""
const CLARKE1866 = Datum("Clarke 1866", 6_378_206.4, 294.978698214)

"""
    WGS72

World Geodetic System 1972.
"""
const WGS72 = Datum("WGS 72", 6_378_135.0, 298.26)
