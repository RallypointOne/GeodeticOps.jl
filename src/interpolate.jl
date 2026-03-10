#--------------------------------------------------------------------------------# Geodesic interpolation

"""
    intermediate_point(a::AbstractFixedCoords, b::AbstractFixedCoords, fraction; datum=DEFAULT_DATUM[]) -> LLA

Return the point at the given `fraction` (0 to 1) along the geodesic from `a` to `b`.

### Examples
```julia
julia> intermediate_point(LLA(0.0, 0.0), LLA(10.0, 0.0), 0.5)
LLA(lon=5.0°, lat=0.0°, alt=0.0m)
```
"""
function intermediate_point(a::AbstractFixedCoords, b::AbstractFixedCoords, fraction; datum::Datum=DEFAULT_DATUM[])
    a_lla = _to_lla(a, datum)
    inv = vincenty_inverse(a_lla, b; datum)
    forward_geodesic(a_lla, inv.azimuth1, inv.distance * fraction; datum)
end

"""
    midpoint(a::AbstractFixedCoords, b::AbstractFixedCoords; datum=DEFAULT_DATUM[]) -> LLA

Return the midpoint along the geodesic from `a` to `b`.

### Examples
```julia
julia> midpoint(LLA(0.0, 0.0), LLA(10.0, 0.0))
LLA(lon=5.0°, lat=0.0°, alt=0.0m)
```
"""
midpoint(a::AbstractFixedCoords, b::AbstractFixedCoords; datum::Datum=DEFAULT_DATUM[]) =
    intermediate_point(a, b, 0.5; datum)
