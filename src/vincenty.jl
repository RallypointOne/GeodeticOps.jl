#--------------------------------------------------------------------------------# Vincenty inverse (geodesic distance + azimuths)

"""
    vincenty_inverse(a::AbstractFixedCoords, b::AbstractFixedCoords; datum=DEFAULT_DATUM[])

Solve the inverse geodesic problem using Vincenty's iterative method.
Returns a `NamedTuple` with:
- `distance`: geodesic distance in meters.
- `azimuth1`: forward azimuth at `a` (degrees, clockwise from north).
- `azimuth2`: reverse azimuth at `b` (degrees, clockwise from north).

Accepts any `AbstractFixedCoords` (LLA, LLARad, ECEF, UTM).

### Examples
```julia
julia> vincenty_inverse(LLA(-74.006, 40.7128), LLA(-0.1278, 51.5074))
(distance = 5.5702e6, azimuth1 = 51.2, azimuth2 = 107.9)
```
"""
function vincenty_inverse(a::AbstractFixedCoords, b::AbstractFixedCoords; datum::Datum=DEFAULT_DATUM[])
    p1 = _to_lla(a, datum)
    p2 = _to_lla(b, datum)

    aa = semi_major_axis(datum)
    bb = semi_minor_axis(datum)
    f = flattening(datum)

    phi1 = deg2rad(p1.lat)
    phi2 = deg2rad(p2.lat)
    L = deg2rad(p2.lon - p1.lon)

    U1 = atan((1 - f) * tan(phi1))
    U2 = atan((1 - f) * tan(phi2))
    sinU1 = sin(U1); cosU1 = cos(U1)
    sinU2 = sin(U2); cosU2 = cos(U2)

    lambda = L
    sinlam = coslam = 0.0
    sin_sigma = cos_sigma = sigma = 0.0
    sin_alpha = cos2_alpha = cos_2sigma_m = 0.0
    converged = false

    for _ in 1:200
        sinlam = sin(lambda)
        coslam = cos(lambda)

        sin_sigma = sqrt(
            (cosU2 * sinlam)^2 +
            (cosU1 * sinU2 - sinU1 * cosU2 * coslam)^2
        )
        sin_sigma == 0.0 && return (distance=0.0, azimuth1=0.0, azimuth2=0.0)

        cos_sigma = sinU1 * sinU2 + cosU1 * cosU2 * coslam
        sigma = atan(sin_sigma, cos_sigma)

        sin_alpha = cosU1 * cosU2 * sinlam / sin_sigma
        cos2_alpha = 1 - sin_alpha^2
        cos_2sigma_m = cos2_alpha == 0.0 ? 0.0 : cos_sigma - 2 * sinU1 * sinU2 / cos2_alpha

        C = f / 16 * cos2_alpha * (4 + f * (4 - 3 * cos2_alpha))
        lambda_prev = lambda
        lambda = L + (1 - C) * f * sin_alpha * (
            sigma + C * sin_sigma * (
                cos_2sigma_m + C * cos_sigma * (-1 + 2 * cos_2sigma_m^2)
            )
        )

        if abs(lambda - lambda_prev) < 1e-12
            converged = true
            break
        end
    end

    if !converged
        d = haversine(p1, p2; datum)
        return (distance=d, azimuth1=NaN, azimuth2=NaN)
    end

    u_sq = cos2_alpha * (aa^2 - bb^2) / bb^2
    A, B = _vincenty_AB(u_sq)
    delta_sigma = _vincenty_delta_sigma(B, sin_sigma, cos_sigma, cos_2sigma_m)

    dist = bb * A * (sigma - delta_sigma)
    az1 = mod(rad2deg(atan(cosU2 * sinlam, cosU1 * sinU2 - sinU1 * cosU2 * coslam)), 360.0)
    az2 = mod(rad2deg(atan(cosU1 * sinlam, -sinU1 * cosU2 + cosU1 * sinU2 * coslam)) + 180.0, 360.0)

    return (distance=dist, azimuth1=az1, azimuth2=az2)
end

"""
    distance(a::AbstractFixedCoords, b::AbstractFixedCoords; datum=DEFAULT_DATUM[]) -> Float64

Compute the geodesic distance in meters between two points using Vincenty's formula.

### Examples
```julia
julia> distance(LLA(-74.006, 40.7128), LLA(-0.1278, 51.5074))
5.5702e6
```
"""
distance(a::AbstractFixedCoords, b::AbstractFixedCoords; datum::Datum=DEFAULT_DATUM[]) =
    vincenty_inverse(a, b; datum).distance

"""
    forward_azimuth(a::AbstractFixedCoords, b::AbstractFixedCoords; datum=DEFAULT_DATUM[]) -> Float64

Compute the forward azimuth (degrees, clockwise from north) from `a` to `b`.

### Examples
```julia
julia> forward_azimuth(LLA(0.0, 0.0), LLA(1.0, 1.0))
44.99563948...
```
"""
forward_azimuth(a::AbstractFixedCoords, b::AbstractFixedCoords; datum::Datum=DEFAULT_DATUM[]) =
    vincenty_inverse(a, b; datum).azimuth1

"""
    bearing(a::AbstractFixedCoords, b::AbstractFixedCoords; datum=DEFAULT_DATUM[]) -> Float64

Alias for [`forward_azimuth`](@ref). Compute the bearing (degrees, clockwise from north) from `a` to `b`.

### Examples
```julia
julia> bearing(LLA(0.0, 0.0), LLA(1.0, 1.0))
44.99563948...
```
"""
const bearing = forward_azimuth
