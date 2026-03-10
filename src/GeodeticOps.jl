module GeodeticOps

export Datum, semi_major_axis, semi_minor_axis, flattening, eccentricity_sq, eccentricity,
    WGS84, GRS80, CLARKE1866, WGS72

export DEFAULT_DATUM

export AbstractCoords, AbstractFixedCoords, AbstractRelativeCoords,
    LLA, LLARad, ECEF, ENU, NED, AER, UTM

export haversine
export vincenty_inverse, distance, forward_azimuth, bearing
export forward_geodesic
export intermediate_point, midpoint
export BoundingBox, boundingbox, destination_boundingbox
export polygon_area

include("datum.jl")
include("coords.jl")
include("utils.jl")
include("haversine.jl")
include("vincenty.jl")
include("forward_geodesic.jl")
include("interpolate.jl")
include("boundingbox.jl")
include("area.jl")

end # module
