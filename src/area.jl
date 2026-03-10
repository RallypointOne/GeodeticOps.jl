#--------------------------------------------------------------------------------# Polygon area on the ellipsoid

"""
    polygon_area(points; datum=DEFAULT_DATUM[]) -> Float64

Compute the area of a polygon on the ellipsoid in square meters. The polygon is defined by
an ordered sequence of `AbstractFixedCoords` vertices. The polygon is implicitly closed (the last
vertex connects back to the first).

Uses the spherical excess method with the authalic (equal-area) sphere radius.

### Examples
```julia
julia> # Approximate area of a 1°×1° square at the equator (~12,300 km²)
       pts = [LLA(0, 0), LLA(1, 0), LLA(1, 1), LLA(0, 1)];

julia> polygon_area(pts) / 1e6  # km²
≈ 12308
```
"""
function polygon_area(points; datum::Datum=DEFAULT_DATUM[])
    n = length(points)
    n < 3 && return 0.0

    # Authalic sphere radius (equal-area)
    a = semi_major_axis(datum)
    esq = eccentricity_sq(datum)
    e = sqrt(esq)
    R = a * sqrt((1 + (1 - esq) / (2e) * log((1 + e) / (1 - e))) / 2)

    # Convert to radians in a single pass (no intermediate llas array)
    lats = Vector{Float64}(undef, n)
    lons = Vector{Float64}(undef, n)
    for (i, p) in enumerate(points)
        lla = _to_lla(p, datum)
        lats[i] = deg2rad(lla.lat)
        lons[i] = deg2rad(lla.lon)
    end

    # Spherical excess via the shoelace-style formula on the sphere
    # Uses the formula: 2A = Σ (λ_{i+1} - λ_{i-1}) * sin(φ_i)
    total = 0.0
    for i in 1:n
        j = mod1(i + 1, n)
        k = mod1(i - 1, n)
        dlon = lons[j] - lons[k]
        # Normalize to [-π, π]
        dlon = atan(sin(dlon), cos(dlon))
        total += dlon * sin(lats[i])
    end

    return abs(total) * R^2 / 2
end
