#--------------------------------------------------------------------------------# Forward geodesic (Vincenty direct)

"""
    forward_geodesic(origin::AbstractFixedCoords, azimuth, dist; datum=DEFAULT_DATUM[]) -> LLA

Solve the forward (direct) geodesic problem: given a starting point, an azimuth
(degrees, clockwise from north), and a distance (meters), compute the destination point.

### Examples
```julia
julia> forward_geodesic(LLA(0.0, 0.0), 90.0, 111_320.0)
LLA(lon=1.0°, lat=0.0°, alt=0.0m)
```
"""
function forward_geodesic(origin::AbstractFixedCoords, azimuth, dist; datum::Datum=DEFAULT_DATUM[])
    p = _to_lla(origin, datum)
    aa = semi_major_axis(datum)
    bb = semi_minor_axis(datum)
    f = flattening(datum)

    phi1 = deg2rad(p.lat)
    alpha1 = deg2rad(azimuth)

    U1 = atan((1 - f) * tan(phi1))
    sinU1 = sin(U1); cosU1 = cos(U1)
    sin_alpha1 = sin(alpha1); cos_alpha1 = cos(alpha1)

    sigma1 = atan(tan(U1), cos_alpha1)
    sin_alpha = cosU1 * sin_alpha1
    cos2_alpha = 1 - sin_alpha^2

    u_sq = cos2_alpha * (aa^2 - bb^2) / bb^2
    A, B = _vincenty_AB(u_sq)

    sigma = dist / (bb * A)
    for _ in 1:200
        cos_2sigma_m = cos(2 * sigma1 + sigma)
        sin_sigma = sin(sigma)
        cos_sigma = cos(sigma)

        delta_sigma = _vincenty_delta_sigma(B, sin_sigma, cos_sigma, cos_2sigma_m)

        sigma_new = dist / (bb * A) + delta_sigma
        abs(sigma_new - sigma) < 1e-12 && break
        sigma = sigma_new
    end

    sin_sigma = sin(sigma)
    cos_sigma = cos(sigma)
    cos_2sigma_m = cos(2 * sigma1 + sigma)

    lat2 = atan(
        sinU1 * cos_sigma + cosU1 * sin_sigma * cos_alpha1,
        (1 - f) * sqrt(sin_alpha^2 + (sinU1 * sin_sigma - cosU1 * cos_sigma * cos_alpha1)^2),
    )

    lam = atan(
        sin_sigma * sin_alpha1,
        cosU1 * cos_sigma - sinU1 * sin_sigma * cos_alpha1,
    )

    C = f / 16 * cos2_alpha * (4 + f * (4 - 3 * cos2_alpha))
    L = lam - (1 - C) * f * sin_alpha * (
        sigma + C * sin_sigma * (
            cos_2sigma_m + C * cos_sigma * (-1 + 2 * cos_2sigma_m^2)
        )
    )

    lon2 = deg2rad(p.lon) + L
    LLA(rad2deg(lon2), rad2deg(lat2), 0.0)
end
