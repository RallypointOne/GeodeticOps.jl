#--------------------------------------------------------------------------------# Bounding box

"""
    BoundingBox(min_lon, min_lat, max_lon, max_lat)

An axis-aligned bounding box in geodetic coordinates (degrees).

When `min_lon > max_lon`, the box wraps across the antimeridian.

### Examples
```julia
julia> bbox = BoundingBox(-75.0, 40.0, -73.0, 41.0)
BoundingBox(lon=[-75.0, -73.0], lat=[40.0, 41.0])
```
"""
struct BoundingBox
    min_lon::Float64
    min_lat::Float64
    max_lon::Float64
    max_lat::Float64
end

function Base.show(io::IO, b::BoundingBox)
    print(io, "BoundingBox(lon=[", b.min_lon, ", ", b.max_lon, "], lat=[", b.min_lat, ", ", b.max_lat, "])")
end

"""
    Base.in(point::AbstractFixedCoords, bbox::BoundingBox; datum=DEFAULT_DATUM[]) -> Bool

Test whether `point` lies within the bounding box.

### Examples
```julia
julia> bbox = BoundingBox(-75.0, 40.0, -73.0, 41.0);

julia> in(LLA(-74.0, 40.5), bbox)
true
```
"""
function Base.in(point::AbstractFixedCoords, bbox::BoundingBox; datum::Datum=DEFAULT_DATUM[])
    p = _to_lla(point, datum)
    lat_ok = bbox.min_lat <= p.lat <= bbox.max_lat
    if bbox.min_lon <= bbox.max_lon
        lon_ok = bbox.min_lon <= p.lon <= bbox.max_lon
    else
        # Wraps across antimeridian
        lon_ok = p.lon >= bbox.min_lon || p.lon <= bbox.max_lon
    end
    return lat_ok && lon_ok
end

"""
    boundingbox(points; datum=DEFAULT_DATUM[]) -> BoundingBox

Compute the axis-aligned bounding box for a collection of `AbstractFixedCoords`.

### Examples
```julia
julia> pts = [LLA(-74.0, 40.0), LLA(-73.0, 41.0), LLA(-75.0, 40.5)];

julia> boundingbox(pts)
BoundingBox(lon=[-75.0, -73.0], lat=[40.0, 41.0])
```
"""
function boundingbox(points; datum::Datum=DEFAULT_DATUM[])
    min_lon = min_lat = Inf
    max_lon = max_lat = -Inf
    for p in points
        lla = _to_lla(p, datum)
        lon, lat = lla.lon, lla.lat
        lon < min_lon && (min_lon = lon)
        lon > max_lon && (max_lon = lon)
        lat < min_lat && (min_lat = lat)
        lat > max_lat && (max_lat = lat)
    end
    BoundingBox(min_lon, min_lat, max_lon, max_lat)
end

"""
    destination_boundingbox(origin::AbstractFixedCoords, dist; datum=DEFAULT_DATUM[]) -> BoundingBox

Compute the bounding box of all points within `dist` meters of `origin`.

### Examples
```julia
julia> destination_boundingbox(LLA(0.0, 0.0), 100_000.0)
BoundingBox(lon=[-0.9, 0.9], lat=[-0.9, 0.9])
```
"""
function destination_boundingbox(origin::AbstractFixedCoords, dist; datum::Datum=DEFAULT_DATUM[])
    north = forward_geodesic(origin, 0.0, dist; datum)
    east = forward_geodesic(origin, 90.0, dist; datum)
    south = forward_geodesic(origin, 180.0, dist; datum)
    west = forward_geodesic(origin, 270.0, dist; datum)
    BoundingBox(west.lon, south.lat, east.lon, north.lat)
end
