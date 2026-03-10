#--------------------------------------------------------------------------------# Default datum

"""
    DEFAULT_DATUM

A `Ref{Datum}` holding the default datum used when `datum` is not specified.
Initialized to `WGS84`.

### Examples
```julia
julia> DEFAULT_DATUM[]
Datum("WGS 84", 6.378137e6, 298.257223563)

julia> DEFAULT_DATUM[] = GRS80  # change the default
```
"""
const DEFAULT_DATUM = Ref{Datum}(WGS84)

#--------------------------------------------------------------------------------# Utilities

_to_lla(c::LLA, ::Datum) = c
_to_lla(c::LLARad, ::Datum) = LLA(c)
_to_lla(c::ECEF, datum::Datum) = LLA(c, datum)
_to_lla(c::UTM, datum::Datum) = LLA(c, datum)

#--------------------------------------------------------------------------------# Shared Vincenty helpers

function _vincenty_AB(u_sq)
    A = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    B = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))
    return A, B
end

function _vincenty_delta_sigma(B, sin_sigma, cos_sigma, cos_2sigma_m)
    return B * sin_sigma * (
        cos_2sigma_m + B / 4 * (
            cos_sigma * (-1 + 2 * cos_2sigma_m^2) -
            B / 6 * cos_2sigma_m * (-3 + 4 * sin_sigma^2) * (-3 + 4 * cos_2sigma_m^2)
        )
    )
end
