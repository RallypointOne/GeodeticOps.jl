[![CI](https://github.com/RallypointOne/GeodeticOps.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/RallypointOne/GeodeticOps.jl/actions/workflows/CI.yml)
[![Docs Build](https://github.com/RallypointOne/GeodeticOps.jl/actions/workflows/Docs.yml/badge.svg)](https://github.com/RallypointOne/GeodeticOps.jl/actions/workflows/Docs.yml)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue)](https://RallypointOne.github.io/GeodeticOps.jl/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue)](https://RallypointOne.github.io/GeodeticOps.jl/dev/)

# GeodeticOps.jl

A pure-Julia package for geodetic computations on Earth's ellipsoid. No external dependencies beyond Julia's standard library.

## Features

- **7 coordinate types** with automatic conversions: `LLA`, `LLARad`, `ECEF`, `ENU`, `NED`, `AER`, `UTM`
- **Geodesic distances**: Haversine (great-circle) and Vincenty (ellipsoidal, sub-millimeter accuracy)
- **Forward and inverse geodesic problems** via Vincenty's formulae
- **Geodesic interpolation**: intermediate points and midpoints along geodesics
- **Polygon area** on the ellipsoid using spherical excess
- **Bounding boxes**: axis-aligned and distance-based
- **Multiple datums**: WGS84 (default), GRS80, Clarke 1866, WGS72

## Quick Start

```julia
using GeodeticOps

nyc = LLA(-74.006, 40.7128)
london = LLA(-0.1278, 51.5074)

# Geodesic distance
distance(nyc, london)  # 5.585e6 meters

# Bearing
bearing(nyc, london)  # 51.2 degrees

# Midpoint
midpoint(nyc, london)  # LLA(lon=-38.8°, lat=52.0°, alt=0.0m)

# Forward geodesic: 1000 km due north from the equator
forward_geodesic(LLA(0.0, 0.0), 0.0, 1_000_000.0)

# Polygon area
pts = [LLA(0, 0), LLA(1, 0), LLA(1, 1), LLA(0, 1)]
polygon_area(pts) / 1e6  # ~12,308 km²

# Coordinate conversions
ecef = ECEF(nyc, WGS84)
utm = UTM(nyc, WGS84)
```

All operations default to WGS84. Override per-call with `datum=GRS80` or change the global default:

```julia
DEFAULT_DATUM[] = GRS80
```

## Coordinate Types

```
AbstractCoords{T}
├── AbstractFixedCoords{T}       # convertible to LLA with only a Datum
│   ├── LLA{T}                   # longitude, latitude (degrees), altitude
│   ├── LLARad{T}                # longitude, latitude (radians), altitude
│   ├── ECEF{T}                  # Earth-Centered, Earth-Fixed (x, y, z)
│   └── UTM{T}                   # Universal Transverse Mercator
└── AbstractRelativeCoords{T}    # require a reference point
    ├── ENU{T}                   # East, North, Up
    ├── NED{T}                   # North, East, Down
    └── AER{T}                   # Azimuth, Elevation, Range
```

All geodetic operations accept any `AbstractFixedCoords`. Relative coordinate types require a reference `LLA` point for conversion.
