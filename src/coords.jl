#--------------------------------------------------------------------------------# AbstractCoords

"""
    AbstractCoords{T}

Abstract supertype for all coordinate representations.
"""
abstract type AbstractCoords{T} end

"""
    AbstractFixedCoords{T} <: AbstractCoords{T}

Abstract supertype for coordinate types that represent a fixed position on the ellipsoid
and can be converted to `LLA` using only a `Datum` (no reference point required).

Subtypes: `LLA`, `LLARad`, `ECEF`, `UTM`.
"""
abstract type AbstractFixedCoords{T} <: AbstractCoords{T} end

"""
    AbstractRelativeCoords{T} <: AbstractCoords{T}

Abstract supertype for coordinate types that are relative to a reference point.
These require both a `Datum` and a reference `LLA` for conversion.

Subtypes: `ENU`, `NED`, `AER`.
"""
abstract type AbstractRelativeCoords{T} <: AbstractCoords{T} end

#-----------------------------------------------------------------------------# LLA
"""
    LLA{T}(lon, lat, alt=0)

Geodetic coordinates: longitude, latitude (degrees), and altitude (meters above ellipsoid).

# Fields
- `lon::T`: Longitude in degrees (−180, 180].
- `lat::T`: Latitude in degrees [−90, 90].
- `alt::T`: Altitude above the ellipsoid in meters.

### Examples
```julia
julia> LLA(0.0, 0.0)
LLA(lon=0.0°, lat=0.0°, alt=0.0m)

julia> LLA{Float32}(1, 2, 3)
LLA(lon=1.0°, lat=2.0°, alt=3.0m)
```
"""
struct LLA{T} <: AbstractFixedCoords{T}
    lon::T
    lat::T
    alt::T
    LLA{T}(lon, lat, alt) where {T} = new{T}(convert(T, lon), convert(T, lat), convert(T, alt))
end

LLA(lon::T, lat::T, alt::T) where {T} = LLA{T}(lon, lat, alt)
LLA(lon, lat, alt) = LLA(promote(lon, lat, alt)...)
LLA{T}(lon, lat) where {T} = LLA{T}(lon, lat, zero(T))
LLA(lon::T, lat::T) where {T} = LLA{T}(lon, lat)
LLA(lon, lat) = LLA(promote(lon, lat)...)

function Base.show(io::IO, c::LLA)
    print(io, "LLA(lon=", c.lon, "°, lat=", c.lat, "°, alt=", c.alt, "m)")
end

#-----------------------------------------------------------------------------# LLARad
"""
    LLARad{T}(lon, lat, alt=0)

Geodetic coordinates: longitude, latitude (radians), and altitude (meters above ellipsoid).

# Fields
- `lon::T`: Longitude in radians (−π, π].
- `lat::T`: Latitude in radians [−π/2, π/2].
- `alt::T`: Altitude above the ellipsoid in meters.

### Examples
```julia
julia> LLARad(0.0, 0.0)
LLARad(lon=0.0 rad, lat=0.0 rad, alt=0.0m)

julia> LLARad(LLA(180.0, 45.0))
LLARad(lon=3.141592653589793 rad, lat=0.7853981633974483 rad, alt=0.0m)
```
"""
struct LLARad{T} <: AbstractFixedCoords{T}
    lon::T
    lat::T
    alt::T
    LLARad{T}(lon, lat, alt) where {T} = new{T}(convert(T, lon), convert(T, lat), convert(T, alt))
end

LLARad(lon::T, lat::T, alt::T) where {T} = LLARad{T}(lon, lat, alt)
LLARad(lon, lat, alt) = LLARad(promote(lon, lat, alt)...)
LLARad{T}(lon, lat) where {T} = LLARad{T}(lon, lat, zero(T))
LLARad(lon::T, lat::T) where {T} = LLARad{T}(lon, lat)
LLARad(lon, lat) = LLARad(promote(lon, lat)...)

LLARad(c::LLA) = LLARad(deg2rad(c.lon), deg2rad(c.lat), c.alt)
LLA(c::LLARad) = LLA(rad2deg(c.lon), rad2deg(c.lat), c.alt)

function Base.show(io::IO, c::LLARad)
    print(io, "LLARad(lon=", c.lon, " rad, lat=", c.lat, " rad, alt=", c.alt, "m)")
end

#-----------------------------------------------------------------------------# ECEF
"""
    ECEF{T}(x, y, z)

Earth-Centered, Earth-Fixed Cartesian coordinates in meters.

# Fields
- `x::T`: X-axis (intersection of prime meridian and equator).
- `y::T`: Y-axis (90° east on the equator).
- `z::T`: Z-axis (north pole).

### Examples
```julia
julia> ECEF(6_378_137.0, 0.0, 0.0)
ECEF(x=6.378137e6, y=0.0, z=0.0)
```
"""
struct ECEF{T} <: AbstractFixedCoords{T}
    x::T
    y::T
    z::T
    ECEF{T}(x, y, z) where {T} = new{T}(convert(T, x), convert(T, y), convert(T, z))
end

ECEF(x::T, y::T, z::T) where {T} = ECEF{T}(x, y, z)
ECEF(x, y, z) = ECEF(promote(x, y, z)...)

function Base.show(io::IO, c::ECEF)
    print(io, "ECEF(x=", c.x, ", y=", c.y, ", z=", c.z, ")")
end

#-----------------------------------------------------------------------------# ENU
"""
    ENU{T}(e, n, u)

East-North-Up local tangent plane coordinates in meters, relative to a reference point.

# Fields
- `e::T`: East displacement in meters.
- `n::T`: North displacement in meters.
- `u::T`: Up displacement in meters.
"""
struct ENU{T} <: AbstractRelativeCoords{T}
    e::T
    n::T
    u::T
    ENU{T}(e, n, u) where {T} = new{T}(convert(T, e), convert(T, n), convert(T, u))
end

ENU(e::T, n::T, u::T) where {T} = ENU{T}(e, n, u)
ENU(e, n, u) = ENU(promote(e, n, u)...)

function Base.show(io::IO, c::ENU)
    print(io, "ENU(e=", c.e, ", n=", c.n, ", u=", c.u, ")")
end

#-----------------------------------------------------------------------------# NED
"""
    NED{T}(n, e, d)

North-East-Down local tangent plane coordinates in meters, relative to a reference point.
Standard aviation/aerospace convention.

# Fields
- `n::T`: North displacement in meters.
- `e::T`: East displacement in meters.
- `d::T`: Down displacement in meters.
"""
struct NED{T} <: AbstractRelativeCoords{T}
    n::T
    e::T
    d::T
    NED{T}(n, e, d) where {T} = new{T}(convert(T, n), convert(T, e), convert(T, d))
end

NED(n::T, e::T, d::T) where {T} = NED{T}(n, e, d)
NED(n, e, d) = NED(promote(n, e, d)...)

NED(c::ENU) = NED(c.n, c.e, -c.u)
ENU(c::NED) = ENU(c.e, c.n, -c.d)

function Base.show(io::IO, c::NED)
    print(io, "NED(n=", c.n, ", e=", c.e, ", d=", c.d, ")")
end

#-----------------------------------------------------------------------------# AER
"""
    AER{T}(az, el, range)

Azimuth-Elevation-Range coordinates relative to a reference point.

# Fields
- `az::T`: Azimuth in degrees, clockwise from north [0, 360).
- `el::T`: Elevation in degrees, up from horizon [−90, 90].
- `range::T`: Slant range in meters.
"""
struct AER{T} <: AbstractRelativeCoords{T}
    az::T
    el::T
    range::T
    AER{T}(az, el, range) where {T} = new{T}(convert(T, az), convert(T, el), convert(T, range))
end

AER(az::T, el::T, range::T) where {T} = AER{T}(az, el, range)
AER(az, el, range) = AER(promote(az, el, range)...)

function Base.show(io::IO, c::AER)
    print(io, "AER(az=", c.az, "°, el=", c.el, "°, range=", c.range, "m)")
end

#-----------------------------------------------------------------------------# UTM
"""
    UTM{T}(zone, isnorth, easting, northing)

Universal Transverse Mercator projected coordinates.

# Fields
- `zone::Int`: UTM zone number (1–60).
- `isnorth::Bool`: `true` for northern hemisphere, `false` for southern.
- `easting::T`: Easting in meters.
- `northing::T`: Northing in meters.

!!! note
    Standard zone determination only. Norway/Svalbard exceptions (zones 32V, 31X–37X) are
    not handled.
"""
struct UTM{T} <: AbstractFixedCoords{T}
    zone::Int
    isnorth::Bool
    easting::T
    northing::T
    UTM{T}(zone, isnorth, easting, northing) where {T} =
        new{T}(Int(zone), Bool(isnorth), convert(T, easting), convert(T, northing))
end

UTM(zone::Integer, isnorth::Bool, easting::T, northing::T) where {T} = UTM{T}(zone, isnorth, easting, northing)
UTM(zone::Integer, isnorth::Bool, easting, northing) = UTM(zone, isnorth, promote(easting, northing)...)

function Base.show(io::IO, c::UTM)
    h = c.isnorth ? "N" : "S"
    print(io, "UTM(zone=", c.zone, h, ", easting=", c.easting, ", northing=", c.northing, ")")
end

#-----------------------------------------------------------------------------# Coordinate conversions

"""
    ECEF(c::LLA, datum::Datum)

Convert geodetic coordinates (degrees) to ECEF Cartesian coordinates.
"""
function ECEF(c::LLA, datum::Datum)
    lon = deg2rad(c.lon)
    lat = deg2rad(c.lat)
    a = semi_major_axis(datum)
    e2 = eccentricity_sq(datum)
    sinlat = sin(lat)
    coslat = cos(lat)
    N = a / sqrt(1 - e2 * sinlat^2)
    ECEF(
        (N + c.alt) * coslat * cos(lon),
        (N + c.alt) * coslat * sin(lon),
        (N * (1 - e2) + c.alt) * sinlat,
    )
end

"""
    ECEF(c::LLARad, datum::Datum)

Convert geodetic coordinates (radians) to ECEF Cartesian coordinates.
"""
ECEF(c::LLARad, datum::Datum) = ECEF(LLA(c), datum)

"""
    LLA(c::ECEF, datum::Datum)

Convert ECEF Cartesian coordinates to geodetic coordinates (degrees).
Uses Bowring's iterative method.
"""
function LLA(c::ECEF, datum::Datum)
    a = semi_major_axis(datum)
    b = semi_minor_axis(datum)
    e2 = eccentricity_sq(datum)
    e2b = (a^2 - b^2) / b^2  # second eccentricity squared

    p = sqrt(c.x^2 + c.y^2)
    lon = atan(c.y, c.x)

    # Bowring's initial estimate
    theta = atan(c.z * a, p * b)
    lat = atan(
        c.z + e2b * b * sin(theta)^3,
        p - e2 * a * cos(theta)^3,
    )

    # Iterate to convergence
    for _ in 1:10
        sinlat = sin(lat)
        N = a / sqrt(1 - e2 * sinlat^2)
        lat_new = atan(c.z + e2 * N * sinlat, p)
        abs(lat_new - lat) < 1e-14 && (lat = lat_new; break)
        lat = lat_new
    end

    sinlat = sin(lat)
    N = a / sqrt(1 - e2 * sinlat^2)
    alt = p / cos(lat) - N

    LLA(rad2deg(lon), rad2deg(lat), alt)
end

"""
    LLARad(c::ECEF, datum::Datum)

Convert ECEF Cartesian coordinates to geodetic coordinates (radians).
"""
LLARad(c::ECEF, datum::Datum) = LLARad(LLA(c, datum))

# Rotation matrix from ECEF to ENU at a given reference latitude/longitude (radians)
function _ecef_to_enu_matrix(lon, lat)
    sinlon = sin(lon)
    coslon = cos(lon)
    sinlat = sin(lat)
    coslat = cos(lat)
    return (
        -sinlon,         coslon,          0.0,
        -sinlat * coslon, -sinlat * sinlon, coslat,
         coslat * coslon,  coslat * sinlon, sinlat,
    )
end

"""
    ENU(c::ECEF, ref::LLA, datum::Datum)

Convert ECEF coordinates to ENU relative to the reference point `ref`.
"""
function ENU(c::ECEF, ref::LLA, datum::Datum)
    ref_ecef = ECEF(ref, datum)
    dx = c.x - ref_ecef.x
    dy = c.y - ref_ecef.y
    dz = c.z - ref_ecef.z
    lon = deg2rad(ref.lon)
    lat = deg2rad(ref.lat)
    R = _ecef_to_enu_matrix(lon, lat)
    ENU(
        R[1] * dx + R[2] * dy + R[3] * dz,
        R[4] * dx + R[5] * dy + R[6] * dz,
        R[7] * dx + R[8] * dy + R[9] * dz,
    )
end

"""
    ENU(c::LLA, ref::LLA, datum::Datum)

Convert LLA coordinates to ENU relative to the reference point `ref`.
"""
ENU(c::LLA, ref::LLA, datum::Datum) = ENU(ECEF(c, datum), ref, datum)

"""
    ECEF(c::ENU, ref::LLA, datum::Datum)

Convert ENU coordinates back to ECEF given the reference point `ref`.
"""
function ECEF(c::ENU, ref::LLA, datum::Datum)
    ref_ecef = ECEF(ref, datum)
    lon = deg2rad(ref.lon)
    lat = deg2rad(ref.lat)
    # Transpose of the ENU rotation matrix = inverse rotation (ECEF to ENU → ENU to ECEF)
    R = _ecef_to_enu_matrix(lon, lat)
    ECEF(
        ref_ecef.x + R[1] * c.e + R[4] * c.n + R[7] * c.u,
        ref_ecef.y + R[2] * c.e + R[5] * c.n + R[8] * c.u,
        ref_ecef.z + R[3] * c.e + R[6] * c.n + R[9] * c.u,
    )
end

"""
    LLA(c::ENU, ref::LLA, datum::Datum)

Convert ENU coordinates to LLA given the reference point `ref`.
"""
LLA(c::ENU, ref::LLA, datum::Datum) = LLA(ECEF(c, ref, datum), datum)

"""
    NED(c::ECEF, ref::LLA, datum::Datum)

Convert ECEF coordinates to NED relative to the reference point `ref`.
"""
NED(c::ECEF, ref::LLA, datum::Datum) = NED(ENU(c, ref, datum))

"""
    NED(c::LLA, ref::LLA, datum::Datum)

Convert LLA coordinates to NED relative to the reference point `ref`.
"""
NED(c::LLA, ref::LLA, datum::Datum) = NED(ENU(c, ref, datum))

"""
    ECEF(c::NED, ref::LLA, datum::Datum)

Convert NED coordinates back to ECEF given the reference point `ref`.
"""
ECEF(c::NED, ref::LLA, datum::Datum) = ECEF(ENU(c), ref, datum)

"""
    LLA(c::NED, ref::LLA, datum::Datum)

Convert NED coordinates to LLA given the reference point `ref`.
"""
LLA(c::NED, ref::LLA, datum::Datum) = LLA(ECEF(c, ref, datum), datum)

#-----------------------------------------------------------------------------# AER ↔ ENU conversions

"""
    ENU(c::AER)

Convert AER to ENU coordinates (no datum needed — pure geometry).
"""
function ENU(c::AER)
    az = deg2rad(c.az)
    el = deg2rad(c.el)
    r_cos_el = c.range * cos(el)
    ENU(
        r_cos_el * sin(az),
        r_cos_el * cos(az),
        c.range * sin(el),
    )
end

"""
    AER(c::ENU)

Convert ENU to AER coordinates (no datum needed — pure geometry).
"""
function AER(c::ENU)
    range = sqrt(c.e^2 + c.n^2 + c.u^2)
    az = mod(rad2deg(atan(c.e, c.n)), 360.0)
    el = rad2deg(atan(c.u, sqrt(c.e^2 + c.n^2)))
    AER(az, el, range)
end

AER(c::NED) = AER(ENU(c))
NED(c::AER) = NED(ENU(c))

"""
    AER(c::LLA, ref::LLA, datum::Datum)

Convert LLA to AER relative to the reference point `ref`.
"""
AER(c::LLA, ref::LLA, datum::Datum) = AER(ENU(c, ref, datum))

"""
    AER(c::ECEF, ref::LLA, datum::Datum)

Convert ECEF to AER relative to the reference point `ref`.
"""
AER(c::ECEF, ref::LLA, datum::Datum) = AER(ENU(c, ref, datum))

"""
    LLA(c::AER, ref::LLA, datum::Datum)

Convert AER to LLA given the reference point `ref`.
"""
LLA(c::AER, ref::LLA, datum::Datum) = LLA(ENU(c), ref, datum)

"""
    ECEF(c::AER, ref::LLA, datum::Datum)

Convert AER to ECEF given the reference point `ref`.
"""
ECEF(c::AER, ref::LLA, datum::Datum) = ECEF(ENU(c), ref, datum)

#-----------------------------------------------------------------------------# UTM conversions

const _UTM_K0 = 0.9996
const _UTM_FALSE_EASTING = 500_000.0
const _UTM_FALSE_NORTHING_SOUTH = 10_000_000.0

"""
    _utm_zone(lon) -> Int

Compute the standard UTM zone number from longitude in degrees.
"""
_utm_zone(lon) = floor(Int, (lon + 180) / 6) + 1

"""
    _utm_central_meridian(zone) -> Float64

Return the central meridian (degrees) for a UTM zone.
"""
_utm_central_meridian(zone) = (zone - 1) * 6.0 - 180.0 + 3.0

# Meridional arc length from equator to latitude φ
function _meridional_arc(a, esq, lat)
    esq2 = esq^2
    esq3 = esq^3
    return a * (
        (1 - esq / 4 - 3 * esq2 / 64 - 5 * esq3 / 256) * lat -
        (3 * esq / 8 + 3 * esq2 / 32 + 45 * esq3 / 1024) * sin(2lat) +
        (15 * esq2 / 256 + 45 * esq3 / 1024) * sin(4lat) -
        (35 * esq3 / 3072) * sin(6lat)
    )
end

"""
    UTM(c::LLA, datum::Datum; zone=nothing)

Convert LLA to UTM. Zone is auto-detected from longitude unless specified.
"""
function UTM(c::LLA, datum::Datum; zone::Union{Int,Nothing}=nothing)
    lat = deg2rad(c.lat)
    z = something(zone, _utm_zone(c.lon))
    lon0 = deg2rad(_utm_central_meridian(z))
    dlon = deg2rad(c.lon) - lon0

    a = semi_major_axis(datum)
    esq = eccentricity_sq(datum)
    epsq = esq / (1 - esq)  # second eccentricity squared

    sinlat = sin(lat)
    coslat = cos(lat)
    tanlat = tan(lat)
    N = a / sqrt(1 - esq * sinlat^2)
    T = tanlat^2
    C = epsq * coslat^2
    A = dlon * coslat
    M = _meridional_arc(a, esq, lat)

    A2 = A^2
    A4 = A2^2
    A6 = A4 * A2

    easting = _UTM_K0 * N * (
        A + (1 - T + C) * A2 * A / 6 +
        (5 - 18T + T^2 + 72C - 58epsq) * A4 * A / 120
    ) + _UTM_FALSE_EASTING

    northing = _UTM_K0 * (
        M + N * tanlat * (
            A2 / 2 +
            (5 - T + 9C + 4C^2) * A4 / 24 +
            (61 - 58T + T^2 + 600C - 330epsq) * A6 / 720
        )
    )

    isnorth = c.lat >= 0
    if !isnorth
        northing += _UTM_FALSE_NORTHING_SOUTH
    end

    UTM(z, isnorth, easting, northing)
end

"""
    LLA(c::UTM, datum::Datum)

Convert UTM to LLA.
"""
function LLA(c::UTM, datum::Datum)
    a = semi_major_axis(datum)
    esq = eccentricity_sq(datum)
    epsq = esq / (1 - esq)
    ee = (1 - sqrt(1 - esq)) / (1 + sqrt(1 - esq))

    x = c.easting - _UTM_FALSE_EASTING
    y = c.isnorth ? c.northing : c.northing - _UTM_FALSE_NORTHING_SOUTH

    M = y / _UTM_K0
    esq2 = esq^2
    esq3 = esq^3
    mu = M / (a * (1 - esq / 4 - 3 * esq2 / 64 - 5 * esq3 / 256))

    ee2 = ee^2
    ee3 = ee^3
    ee4 = ee^4
    phi1 = mu +
        (3 * ee / 2 - 27 * ee3 / 32) * sin(2mu) +
        (21 * ee2 / 16 - 55 * ee4 / 32) * sin(4mu) +
        (151 * ee3 / 96) * sin(6mu) +
        (1097 * ee4 / 512) * sin(8mu)

    sinphi1 = sin(phi1)
    cosphi1 = cos(phi1)
    tanphi1 = tan(phi1)
    N1 = a / sqrt(1 - esq * sinphi1^2)
    T1 = tanphi1^2
    C1 = epsq * cosphi1^2
    R1 = a * (1 - esq) / (1 - esq * sinphi1^2)^1.5
    D = x / (N1 * _UTM_K0)

    D2 = D^2
    D4 = D2^2
    D6 = D4 * D2

    lat = phi1 - (N1 * tanphi1 / R1) * (
        D2 / 2 -
        (5 + 3T1 + 10C1 - 4C1^2 - 9epsq) * D4 / 24 +
        (61 + 90T1 + 298C1 + 45T1^2 - 252epsq - 3C1^2) * D6 / 720
    )

    lon0 = deg2rad(_utm_central_meridian(c.zone))
    lon = lon0 + (
        D -
        (1 + 2T1 + C1) * D2 * D / 6 +
        (5 - 2C1 + 28T1 - 3C1^2 + 8epsq + 24T1^2) * D4 * D / 120
    ) / cosphi1

    LLA(rad2deg(lon), rad2deg(lat), 0.0)
end

ECEF(c::UTM, datum::Datum) = ECEF(LLA(c, datum), datum)
UTM(c::ECEF, datum::Datum; zone::Union{Int,Nothing}=nothing) = UTM(LLA(c, datum), datum; zone=zone)
