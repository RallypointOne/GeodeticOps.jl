#--------------------------------------------------------------------------------# Haversine distance

"""
    haversine(a::AbstractFixedCoords, b::AbstractFixedCoords; datum=DEFAULT_DATUM[]) -> Float64

Compute the great-circle distance in meters between two coordinates using the
haversine formula. Uses the mean radius `R = (2a + b) / 3` of the ellipsoid defined by
`datum`. Altitude is ignored.

Accepts any `AbstractFixedCoords` (LLA, LLARad, ECEF, UTM).

### Examples
```julia
julia> haversine(LLA(-74.006, 40.7128), LLA(-0.1278, 51.5074))  # ≈ 5570 km
5.570222e6
```
"""
haversine(a::AbstractFixedCoords, b::AbstractFixedCoords; datum::Datum=DEFAULT_DATUM[]) =
    haversine(_to_lla(a, datum), _to_lla(b, datum); datum)

function haversine(a::LLA, b::LLA; datum::Datum=DEFAULT_DATUM[])
    R = (2 * semi_major_axis(datum) + semi_minor_axis(datum)) / 3

    lat1 = deg2rad(a.lat)
    lat2 = deg2rad(b.lat)
    dlat = lat2 - lat1
    dlon = deg2rad(b.lon - a.lon)

    h = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
    return 2 * R * asin(sqrt(h))
end
