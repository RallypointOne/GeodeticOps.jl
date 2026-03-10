using GeodeticOps
using Test

@testset "GeodeticOps.jl" begin
    @testset "Datum construction" begin
        d = Datum("Test", 6_378_137.0, 298.257223563)
        @test d.name == "Test"
        @test d.a == 6_378_137.0
        @test d.f_inv == 298.257223563
    end

    @testset "WGS84 derived quantities" begin
        @test semi_major_axis(WGS84) == 6_378_137.0
        @test flattening(WGS84) ≈ 1.0 / 298.257223563
        @test semi_minor_axis(WGS84) ≈ 6_356_752.314245179 atol=1e-4
        @test eccentricity_sq(WGS84) ≈ 0.00669437999014 atol=1e-11
        @test eccentricity(WGS84) ≈ 0.0818191908426 atol=1e-10
    end

    @testset "GRS80 vs WGS84" begin
        # GRS80 and WGS84 share the same semi-major axis but differ in flattening
        @test semi_major_axis(GRS80) == semi_major_axis(WGS84)
        @test GRS80.f_inv != WGS84.f_inv
    end

    @testset "CLARKE1866 derived quantities" begin
        @test semi_major_axis(CLARKE1866) == 6_378_206.4
        @test flattening(CLARKE1866) ≈ 1.0 / 294.978698214
        @test semi_minor_axis(CLARKE1866) < semi_major_axis(CLARKE1866)
        @test eccentricity(CLARKE1866) > 0.0
    end

    @testset "WGS72 derived quantities" begin
        @test semi_major_axis(WGS72) == 6_378_135.0
        @test WGS72.f_inv == 298.26
        @test semi_minor_axis(WGS72) < semi_major_axis(WGS72)
    end

    @testset "show method" begin
        s = sprint(show, WGS84)
        @test contains(s, "WGS 84")
        @test contains(s, "6.378137e6")
    end

    @testset "DEFAULT_DATUM" begin
        @test DEFAULT_DATUM[] === WGS84
        DEFAULT_DATUM[] = GRS80
        @test DEFAULT_DATUM[] === GRS80
        DEFAULT_DATUM[] = WGS84  # restore
        @test DEFAULT_DATUM[] === WGS84
    end

    @testset "LLA" begin
        @testset "construction" begin
            c = LLA(1.0, 2.0, 3.0)
            @test c.lon == 1.0
            @test c.lat == 2.0
            @test c.alt == 3.0
            @test c isa LLA{Float64}
        end
        @testset "default alt=0" begin
            c = LLA(1.0, 2.0)
            @test c.alt == 0.0
        end
        @testset "type promotion" begin
            c = LLA(1, 2, 3)
            @test c isa LLA{Int}
            c_mixed = LLA(1, 2.0, 3)
            @test c_mixed isa LLA{Float64}
            c32 = LLA(1.0f0, 2.0f0)
            @test c32 isa LLA{Float32}
        end
        @testset "explicit type parameter" begin
            c = LLA{Float32}(1, 2, 3)
            @test c isa LLA{Float32}
            @test c.lon === 1.0f0
        end
        @testset "show" begin
            s = sprint(show, LLA(1.5, -2.5, 100.0))
            @test contains(s, "lon=1.5")
            @test contains(s, "lat=-2.5")
            @test contains(s, "alt=100.0")
        end
    end

    @testset "LLARad" begin
        @testset "construction" begin
            c = LLARad(1.0, 2.0, 3.0)
            @test c.lon == 1.0
            @test c.lat == 2.0
            @test c.alt == 3.0
            @test c isa LLARad{Float64}
        end
        @testset "default alt=0" begin
            c = LLARad(1.0, 2.0)
            @test c.alt == 0.0
        end
        @testset "type promotion" begin
            c = LLARad(1, 2, 3)
            @test c isa LLARad{Int}
            c32 = LLARad(1.0f0, 2.0f0)
            @test c32 isa LLARad{Float32}
        end
        @testset "explicit type parameter" begin
            c = LLARad{Float32}(1, 2, 3)
            @test c isa LLARad{Float32}
        end
        @testset "convert from LLA" begin
            lla = LLA(180.0, 45.0, 100.0)
            r = LLARad(lla)
            @test r.lon ≈ π
            @test r.lat ≈ π / 4
            @test r.alt == 100.0
        end
        @testset "convert to LLA" begin
            r = LLARad(π, π / 4, 100.0)
            lla = LLA(r)
            @test lla.lon ≈ 180.0
            @test lla.lat ≈ 45.0
            @test lla.alt == 100.0
        end
        @testset "round-trip" begin
            orig = LLA(-74.006, 40.7128, 50.0)
            rt = LLA(LLARad(orig))
            @test rt.lon ≈ orig.lon
            @test rt.lat ≈ orig.lat
            @test rt.alt == orig.alt
        end
        @testset "show" begin
            s = sprint(show, LLARad(1.5, -0.5, 100.0))
            @test contains(s, "lon=1.5 rad")
            @test contains(s, "lat=-0.5 rad")
            @test contains(s, "alt=100.0")
        end
    end

    @testset "ENU" begin
        @testset "construction" begin
            c = ENU(1.0, 2.0, 3.0)
            @test c.e == 1.0
            @test c.n == 2.0
            @test c.u == 3.0
            @test c isa ENU{Float64}
        end
        @testset "type promotion" begin
            c = ENU(1, 2, 3)
            @test c isa ENU{Int}
            c32 = ENU(1.0f0, 2.0f0, 3.0f0)
            @test c32 isa ENU{Float32}
        end
        @testset "explicit type parameter" begin
            c = ENU{Float32}(1, 2, 3)
            @test c isa ENU{Float32}
            @test c.e === 1.0f0
        end
        @testset "show" begin
            s = sprint(show, ENU(1.0, 2.0, 3.0))
            @test contains(s, "e=1.0")
            @test contains(s, "n=2.0")
            @test contains(s, "u=3.0")
        end
    end

    @testset "NED" begin
        @testset "construction" begin
            c = NED(1.0, 2.0, 3.0)
            @test c.n == 1.0
            @test c.e == 2.0
            @test c.d == 3.0
            @test c isa NED{Float64}
        end
        @testset "type promotion" begin
            c = NED(1, 2, 3)
            @test c isa NED{Int}
            c_mixed = NED(1, 2.0, 3)
            @test c_mixed isa NED{Float64}
        end
        @testset "explicit type parameter" begin
            c = NED{Float32}(1, 2, 3)
            @test c isa NED{Float32}
            @test c.n === 1.0f0
        end
        @testset "ENU ↔ NED" begin
            enu = ENU(10.0, 20.0, 30.0)
            ned = NED(enu)
            @test ned.n == 20.0
            @test ned.e == 10.0
            @test ned.d == -30.0
            enu2 = ENU(ned)
            @test enu2.e == enu.e
            @test enu2.n == enu.n
            @test enu2.u == enu.u
        end
        @testset "show" begin
            s = sprint(show, NED(1.0, 2.0, 3.0))
            @test contains(s, "n=1.0")
            @test contains(s, "e=2.0")
            @test contains(s, "d=3.0")
        end
    end

    @testset "AER" begin
        @testset "construction" begin
            c = AER(45.0, 30.0, 1000.0)
            @test c.az == 45.0
            @test c.el == 30.0
            @test c.range == 1000.0
            @test c isa AER{Float64}
        end
        @testset "type promotion" begin
            c = AER(45, 30, 1000)
            @test c isa AER{Int}
            c_mixed = AER(45, 30.0, 1000)
            @test c_mixed isa AER{Float64}
        end
        @testset "explicit type parameter" begin
            c = AER{Float32}(45, 30, 1000)
            @test c isa AER{Float32}
            @test c.az === 45.0f0
        end
        @testset "show" begin
            s = sprint(show, AER(45.0, 30.0, 1000.0))
            @test contains(s, "az=45.0")
            @test contains(s, "el=30.0")
            @test contains(s, "range=1000.0")
        end
    end

    @testset "UTM" begin
        @testset "construction" begin
            c = UTM(18, true, 583960.0, 4507523.0)
            @test c.zone == 18
            @test c.isnorth == true
            @test c.easting == 583960.0
            @test c.northing == 4507523.0
            @test c isa UTM{Float64}
        end
        @testset "type promotion" begin
            c = UTM(18, true, 583960, 4507523)
            @test c isa UTM{Int}
            c_mixed = UTM(18, true, 583960, 4507523.0)
            @test c_mixed isa UTM{Float64}
        end
        @testset "explicit type parameter" begin
            c = UTM{Float32}(18, true, 583960, 4507523)
            @test c isa UTM{Float32}
            @test c.easting === 583960.0f0
        end
        @testset "show" begin
            s = sprint(show, UTM(18, true, 583960.0, 4507523.0))
            @test contains(s, "zone=18N")
            @test contains(s, "easting=583960.0")
            s2 = sprint(show, UTM(18, false, 583960.0, 4507523.0))
            @test contains(s2, "zone=18S")
        end
    end

    @testset "AbstractCoords" begin
        @test LLA{Float64} <: AbstractCoords{Float64}
        @test LLARad{Float64} <: AbstractCoords{Float64}
        @test ECEF{Float32} <: AbstractCoords{Float32}
        @test ENU{Float64} <: AbstractCoords{Float64}
        @test NED{Float64} <: AbstractCoords{Float64}
        @test AER{Float64} <: AbstractCoords{Float64}
        @test UTM{Float64} <: AbstractCoords{Float64}
        @test LLA{Float64} <: AbstractCoords
        @test ECEF{Float64} <: AbstractCoords
    end

    @testset "AbstractFixedCoords / AbstractRelativeCoords" begin
        @test LLA{Float64} <: AbstractFixedCoords{Float64}
        @test LLARad{Float64} <: AbstractFixedCoords{Float64}
        @test ECEF{Float64} <: AbstractFixedCoords{Float64}
        @test UTM{Float64} <: AbstractFixedCoords{Float64}
        @test ENU{Float64} <: AbstractRelativeCoords{Float64}
        @test NED{Float64} <: AbstractRelativeCoords{Float64}
        @test AER{Float64} <: AbstractRelativeCoords{Float64}
        @test !(ENU{Float64} <: AbstractFixedCoords)
        @test !(LLA{Float64} <: AbstractRelativeCoords)
    end

    @testset "ECEF" begin
        @testset "construction" begin
            c = ECEF(1.0, 2.0, 3.0)
            @test c.x == 1.0
            @test c.y == 2.0
            @test c.z == 3.0
            @test c isa ECEF{Float64}
        end
        @testset "type promotion" begin
            c = ECEF(1, 2, 3)
            @test c isa ECEF{Int}
            c_mixed = ECEF(1, 2.0, 3)
            @test c_mixed isa ECEF{Float64}
            c32 = ECEF(1.0f0, 2.0f0, 3.0f0)
            @test c32 isa ECEF{Float32}
        end
        @testset "explicit type parameter" begin
            c = ECEF{Float32}(1, 2, 3)
            @test c isa ECEF{Float32}
            @test c.x === 1.0f0
        end
        @testset "show" begin
            s = sprint(show, ECEF(1.0, 2.0, 3.0))
            @test contains(s, "x=1.0")
            @test contains(s, "y=2.0")
            @test contains(s, "z=3.0")
        end
    end

    @testset "Coordinate conversions" begin
        @testset "LLA → ECEF known point (equator/prime meridian)" begin
            # On equator at prime meridian, alt=0: x=a, y=0, z=0
            lla = LLA(0.0, 0.0, 0.0)
            ecef = ECEF(lla, WGS84)
            @test ecef.x ≈ semi_major_axis(WGS84)
            @test abs(ecef.y) < 1e-6
            @test abs(ecef.z) < 1e-6
        end
        @testset "LLA → ECEF north pole" begin
            lla = LLA(0.0, 90.0, 0.0)
            ecef = ECEF(lla, WGS84)
            @test abs(ecef.x) < 1e-6
            @test abs(ecef.y) < 1e-6
            @test ecef.z ≈ semi_minor_axis(WGS84)
        end
        @testset "LLA → ECEF south pole" begin
            lla = LLA(0.0, -90.0, 0.0)
            ecef = ECEF(lla, WGS84)
            @test abs(ecef.x) < 1e-6
            @test abs(ecef.y) < 1e-6
            @test ecef.z ≈ -semi_minor_axis(WGS84)
        end
        @testset "LLA → ECEF 90° longitude" begin
            lla = LLA(90.0, 0.0, 0.0)
            ecef = ECEF(lla, WGS84)
            @test abs(ecef.x) < 1e-6
            @test ecef.y ≈ semi_major_axis(WGS84)
            @test abs(ecef.z) < 1e-6
        end
        @testset "LLA → ECEF with non-default datum" begin
            lla = LLA(0.0, 0.0, 0.0)
            ecef = ECEF(lla, CLARKE1866)
            @test ecef.x ≈ semi_major_axis(CLARKE1866)
            @test ecef.x != semi_major_axis(WGS84)
        end
        @testset "LLA → ECEF → LLA round-trip" begin
            orig = LLA(-74.006, 40.7128, 10.0)
            rt = LLA(ECEF(orig, WGS84), WGS84)
            @test rt.lon ≈ orig.lon atol=1e-10
            @test rt.lat ≈ orig.lat atol=1e-10
            @test rt.alt ≈ orig.alt atol=1e-6
        end
        @testset "LLA → ECEF → LLA round-trip (high altitude)" begin
            orig = LLA(139.6917, 35.6895, 35_786_000.0)  # geostationary
            rt = LLA(ECEF(orig, WGS84), WGS84)
            @test rt.lon ≈ orig.lon atol=1e-10
            @test rt.lat ≈ orig.lat atol=1e-10
            @test rt.alt ≈ orig.alt atol=1e-4
        end
        @testset "LLARad → ECEF" begin
            lla = LLA(0.0, 0.0, 0.0)
            ecef_from_lla = ECEF(lla, WGS84)
            ecef_from_rad = ECEF(LLARad(lla), WGS84)
            @test ecef_from_rad.x ≈ ecef_from_lla.x
            @test ecef_from_rad.y ≈ ecef_from_lla.y
            @test ecef_from_rad.z ≈ ecef_from_lla.z
        end
        @testset "ECEF → LLARad" begin
            orig = LLA(-74.006, 40.7128, 10.0)
            ecef = ECEF(orig, WGS84)
            r = LLARad(ecef, WGS84)
            @test r.lon ≈ deg2rad(orig.lon) atol=1e-10
            @test r.lat ≈ deg2rad(orig.lat) atol=1e-10
            @test r.alt ≈ orig.alt atol=1e-6
        end
        @testset "ENU: reference point is origin" begin
            ref = LLA(-74.006, 40.7128, 10.0)
            enu = ENU(ref, ref, WGS84)
            @test abs(enu.e) < 1e-6
            @test abs(enu.n) < 1e-6
            @test abs(enu.u) < 1e-6
        end
        @testset "ENU: point due east" begin
            # Small offset in longitude at equator ≈ east displacement
            ref = LLA(0.0, 0.0, 0.0)
            east_point = LLA(0.001, 0.0, 0.0)  # ~111 m east
            enu = ENU(east_point, ref, WGS84)
            @test enu.e > 100.0   # should be ~111 m east
            @test abs(enu.n) < 1.0
            @test abs(enu.u) < 1.0
        end
        @testset "ENU: point due north" begin
            ref = LLA(0.0, 0.0, 0.0)
            north_point = LLA(0.0, 0.001, 0.0)  # ~111 m north
            enu = ENU(north_point, ref, WGS84)
            @test abs(enu.e) < 1.0
            @test enu.n > 100.0   # should be ~111 m north
            @test abs(enu.u) < 1.0
        end
        @testset "LLA → ENU → LLA round-trip" begin
            ref = LLA(-74.006, 40.7128, 0.0)
            orig = LLA(-73.996, 40.7228, 50.0)
            rt = LLA(ENU(orig, ref, WGS84), ref, WGS84)
            @test rt.lon ≈ orig.lon atol=1e-10
            @test rt.lat ≈ orig.lat atol=1e-10
            @test rt.alt ≈ orig.alt atol=1e-4
        end
        @testset "NED conversions" begin
            ref = LLA(0.0, 0.0, 0.0)
            point = LLA(0.001, 0.001, 100.0)
            enu = ENU(point, ref, WGS84)
            ned = NED(point, ref, WGS84)
            @test ned.n ≈ enu.n
            @test ned.e ≈ enu.e
            @test ned.d ≈ -enu.u
        end
        @testset "LLA → NED → LLA round-trip" begin
            ref = LLA(-74.006, 40.7128, 0.0)
            orig = LLA(-73.996, 40.7228, 50.0)
            rt = LLA(NED(orig, ref, WGS84), ref, WGS84)
            @test rt.lon ≈ orig.lon atol=1e-10
            @test rt.lat ≈ orig.lat atol=1e-10
            @test rt.alt ≈ orig.alt atol=1e-4
        end
        @testset "AER ↔ ENU" begin
            # Due north, 45° elevation, 1000m range
            aer = AER(0.0, 45.0, 1000.0)
            enu = ENU(aer)
            @test abs(enu.e) < 1e-10
            @test enu.n ≈ 1000.0 * cos(deg2rad(45.0))
            @test enu.u ≈ 1000.0 * sin(deg2rad(45.0))
            # Round-trip
            aer2 = AER(enu)
            @test aer2.az ≈ aer.az atol=1e-10
            @test aer2.el ≈ aer.el atol=1e-10
            @test aer2.range ≈ aer.range atol=1e-10
        end
        @testset "AER ↔ ENU (due east)" begin
            aer = AER(90.0, 0.0, 500.0)
            enu = ENU(aer)
            @test enu.e ≈ 500.0
            @test abs(enu.n) < 1e-10
            @test abs(enu.u) < 1e-10
        end
        @testset "AER ↔ NED" begin
            aer = AER(90.0, 0.0, 500.0)
            ned = NED(aer)
            @test abs(ned.n) < 1e-10
            @test ned.e ≈ 500.0
            @test abs(ned.d) < 1e-10
        end
        @testset "LLA → AER → LLA round-trip" begin
            ref = LLA(-74.006, 40.7128, 0.0)
            orig = LLA(-73.996, 40.7228, 50.0)
            rt = LLA(AER(orig, ref, WGS84), ref, WGS84)
            @test rt.lon ≈ orig.lon atol=1e-10
            @test rt.lat ≈ orig.lat atol=1e-10
            @test rt.alt ≈ orig.alt atol=1e-4
        end
        @testset "UTM: NYC (zone 18N)" begin
            nyc = LLA(-74.006, 40.7128)
            utm = UTM(nyc, WGS84)
            @test utm.zone == 18
            @test utm.isnorth == true
            @test utm.easting ≈ 583959.4 atol=1.0
            @test utm.northing ≈ 4_507_351.0 atol=1.0
        end
        @testset "UTM: southern hemisphere" begin
            sydney = LLA(151.2093, -33.8688)
            utm = UTM(sydney, WGS84)
            @test utm.zone == 56
            @test utm.isnorth == false
        end
        @testset "LLA → UTM → LLA round-trip" begin
            orig = LLA(-74.006, 40.7128)
            utm = UTM(orig, WGS84)
            rt = LLA(utm, WGS84)
            @test rt.lon ≈ orig.lon atol=1e-6
            @test rt.lat ≈ orig.lat atol=1e-6
        end
        @testset "LLA → UTM → LLA round-trip (southern)" begin
            orig = LLA(151.2093, -33.8688)
            utm = UTM(orig, WGS84)
            rt = LLA(utm, WGS84)
            @test rt.lon ≈ orig.lon atol=1e-6
            @test rt.lat ≈ orig.lat atol=1e-6
        end
        @testset "UTM: explicit zone override" begin
            lla = LLA(-74.006, 40.7128)
            utm17 = UTM(lla, WGS84; zone=17)
            @test utm17.zone == 17
        end
        @testset "ECEF → UTM" begin
            orig = LLA(-74.006, 40.7128)
            ecef = ECEF(orig, WGS84)
            utm = UTM(ecef, WGS84)
            @test utm.zone == 18
            @test utm.easting ≈ UTM(orig, WGS84).easting atol=1e-6
        end
        @testset "UTM → ECEF" begin
            orig = LLA(-74.006, 40.7128)
            utm = UTM(orig, WGS84)
            ecef_from_utm = ECEF(utm, WGS84)
            ecef_from_lla = ECEF(orig, WGS84)
            @test ecef_from_utm.x ≈ ecef_from_lla.x atol=1.0
            @test ecef_from_utm.y ≈ ecef_from_lla.y atol=1.0
            @test ecef_from_utm.z ≈ ecef_from_lla.z atol=1.0
        end
        @testset "UTM zone at 180° longitude" begin
            lla = LLA(179.5, 10.0)
            utm = UTM(lla, WGS84)
            @test utm.zone == 60
        end
        @testset "ENU: point due up" begin
            ref = LLA(0.0, 0.0, 0.0)
            up_point = LLA(0.0, 0.0, 1000.0)
            enu = ENU(up_point, ref, WGS84)
            @test abs(enu.e) < 1.0
            @test abs(enu.n) < 1.0
            @test enu.u ≈ 1000.0 atol=1.0
        end
        @testset "ECEF → ENU direct" begin
            ref = LLA(-74.006, 40.7128, 0.0)
            point = LLA(-73.996, 40.7228, 50.0)
            ecef = ECEF(point, WGS84)
            enu_from_lla = ENU(point, ref, WGS84)
            enu_from_ecef = ENU(ecef, ref, WGS84)
            @test enu_from_ecef.e ≈ enu_from_lla.e atol=1e-6
            @test enu_from_ecef.n ≈ enu_from_lla.n atol=1e-6
            @test enu_from_ecef.u ≈ enu_from_lla.u atol=1e-6
        end
        @testset "ECEF → NED direct" begin
            ref = LLA(0.0, 0.0, 0.0)
            point = LLA(0.001, 0.001, 100.0)
            ecef = ECEF(point, WGS84)
            ned_from_lla = NED(point, ref, WGS84)
            ned_from_ecef = NED(ecef, ref, WGS84)
            @test ned_from_ecef.n ≈ ned_from_lla.n atol=1e-6
            @test ned_from_ecef.e ≈ ned_from_lla.e atol=1e-6
            @test ned_from_ecef.d ≈ ned_from_lla.d atol=1e-6
        end
        @testset "NED → ECEF" begin
            ref = LLA(-74.006, 40.7128, 0.0)
            orig = LLA(-73.996, 40.7228, 50.0)
            ned = NED(orig, ref, WGS84)
            ecef_rt = ECEF(ned, ref, WGS84)
            ecef_orig = ECEF(orig, WGS84)
            @test ecef_rt.x ≈ ecef_orig.x atol=1e-4
            @test ecef_rt.y ≈ ecef_orig.y atol=1e-4
            @test ecef_rt.z ≈ ecef_orig.z atol=1e-4
        end
        @testset "AER → ECEF" begin
            ref = LLA(-74.006, 40.7128, 0.0)
            orig = LLA(-73.996, 40.7228, 50.0)
            aer = AER(orig, ref, WGS84)
            ecef_rt = ECEF(aer, ref, WGS84)
            ecef_orig = ECEF(orig, WGS84)
            @test ecef_rt.x ≈ ecef_orig.x atol=1e-4
            @test ecef_rt.y ≈ ecef_orig.y atol=1e-4
            @test ecef_rt.z ≈ ecef_orig.z atol=1e-4
        end
        @testset "AER from ECEF" begin
            ref = LLA(0.0, 0.0, 0.0)
            point = LLA(0.001, 0.0, 0.0)
            aer_from_lla = AER(point, ref, WGS84)
            ecef = ECEF(point, WGS84)
            aer_from_ecef = AER(ecef, ref, WGS84)
            @test aer_from_ecef.az ≈ aer_from_lla.az atol=1e-6
            @test aer_from_ecef.el ≈ aer_from_lla.el atol=1e-6
            @test aer_from_ecef.range ≈ aer_from_lla.range atol=1e-6
        end
        @testset "AER ↔ NED round-trip" begin
            aer = AER(45.0, 30.0, 1000.0)
            ned = NED(aer)
            aer2 = AER(ned)
            @test aer2.az ≈ aer.az atol=1e-10
            @test aer2.el ≈ aer.el atol=1e-10
            @test aer2.range ≈ aer.range atol=1e-10
        end
        @testset "LLARad → ECEF → LLARad round-trip" begin
            orig = LLARad(deg2rad(-74.006), deg2rad(40.7128), 10.0)
            ecef = ECEF(orig, WGS84)
            rt = LLARad(ecef, WGS84)
            @test rt.lon ≈ orig.lon atol=1e-10
            @test rt.lat ≈ orig.lat atol=1e-10
            @test rt.alt ≈ orig.alt atol=1e-6
        end
        @testset "LLA → ECEF → LLA round-trip (multiple locations)" begin
            locations = [
                LLA(139.6917, 35.6895, 0.0),    # Tokyo
                LLA(-43.1729, -22.9068, 0.0),    # Rio
                LLA(37.6173, 55.7558, 200.0),    # Moscow
                LLA(-180.0, 0.0, 0.0),           # Antimeridian
                LLA(0.0, -90.0, 0.0),            # South pole
            ]
            for orig in locations
                rt = LLA(ECEF(orig, WGS84), WGS84)
                @test rt.lon ≈ orig.lon atol=1e-8
                @test rt.lat ≈ orig.lat atol=1e-8
                @test rt.alt ≈ orig.alt atol=1e-4
            end
        end
        @testset "UTM with GRS80 datum" begin
            nyc = LLA(-74.006, 40.7128)
            utm_wgs84 = UTM(nyc, WGS84)
            utm_grs80 = UTM(nyc, GRS80)
            # GRS80 and WGS84 are very close, so UTM should be nearly identical
            @test utm_wgs84.zone == utm_grs80.zone
            @test utm_wgs84.easting ≈ utm_grs80.easting atol=0.01
        end
    end

    @testset "haversine" begin
        @testset "NYC to London" begin
            nyc = LLA(-74.006, 40.7128)
            london = LLA(-0.1278, 51.5074)
            d = haversine(nyc, london)
            @test d ≈ 5.570e6 rtol=0.01  # ~5570 km within 1%
        end
        @testset "same point is zero" begin
            p = LLA(10.0, 20.0)
            @test haversine(p, p) == 0.0
        end
        @testset "antipodal points" begin
            # North pole to south pole ≈ half circumference ≈ π * R
            np = LLA(0.0, 90.0)
            sp = LLA(0.0, -90.0)
            R = (2 * semi_major_axis(WGS84) + semi_minor_axis(WGS84)) / 3
            @test haversine(np, sp) ≈ π * R rtol=1e-10
        end
        @testset "symmetry" begin
            a = LLA(-122.4194, 37.7749)  # San Francisco
            b = LLA(139.6917, 35.6895)   # Tokyo
            @test haversine(a, b) ≈ haversine(b, a)
        end
        @testset "accepts AbstractFixedCoords" begin
            nyc_lla = LLA(-74.006, 40.7128)
            nyc_ecef = ECEF(nyc_lla, WGS84)
            london = LLA(-0.1278, 51.5074)
            @test haversine(nyc_ecef, london) ≈ haversine(nyc_lla, london) atol=1.0
        end
        @testset "explicit datum kwarg" begin
            a = LLA(-74.006, 40.7128)
            b = LLA(-0.1278, 51.5074)
            @test haversine(a, b; datum=GRS80) ≈ haversine(a, b) rtol=1e-6
        end
        @testset "short distance" begin
            # Two points ~11m apart
            a = LLA(0.0, 0.0)
            b = LLA(0.0001, 0.0)
            d = haversine(a, b)
            @test d ≈ 11.1 atol=1.0
        end
        @testset "accepts UTM" begin
            nyc_lla = LLA(-74.006, 40.7128)
            nyc_utm = UTM(nyc_lla, WGS84)
            london = LLA(-0.1278, 51.5074)
            @test haversine(nyc_utm, london) ≈ haversine(nyc_lla, london) atol=10.0
        end
        @testset "accepts LLARad" begin
            a_lla = LLA(-74.006, 40.7128)
            a_rad = LLARad(a_lla)
            b = LLA(-0.1278, 51.5074)
            @test haversine(a_rad, b) ≈ haversine(a_lla, b) atol=1.0
        end
    end

    @testset "vincenty_inverse" begin
        @testset "Flinders Peak to Buninyong" begin
            # Standard Vincenty test case
            fp = LLA(144.42486789, -37.95103342)
            bu = LLA(143.92649553, -37.65282114)
            inv = vincenty_inverse(fp, bu)
            @test inv.distance ≈ 54972.271 atol=0.1
            @test inv.azimuth1 ≈ 306.8686 atol=0.01
            @test inv.azimuth2 ≈ 127.1107 atol=0.1
        end
        @testset "same point" begin
            p = LLA(10.0, 20.0)
            inv = vincenty_inverse(p, p)
            @test inv.distance == 0.0
        end
        @testset "symmetry" begin
            a = LLA(-74.006, 40.7128)
            b = LLA(-0.1278, 51.5074)
            @test vincenty_inverse(a, b).distance ≈ vincenty_inverse(b, a).distance
        end
        @testset "agrees with haversine" begin
            a = LLA(-74.006, 40.7128)
            b = LLA(-0.1278, 51.5074)
            vh = haversine(a, b)
            vd = distance(a, b)
            @test vd ≈ vh rtol=0.005  # within 0.5%
        end
        @testset "accepts AbstractFixedCoords" begin
            a_ecef = ECEF(LLA(-74.006, 40.7128), WGS84)
            b_rad = LLARad(LLA(-0.1278, 51.5074))
            d = distance(a_ecef, b_rad)
            @test d > 5.5e6
        end
        @testset "near-antipodal fallback" begin
            # Nearly antipodal points may not converge; should fall back to haversine
            a = LLA(0.0, 89.999)
            b = LLA(180.0, -89.999)
            inv = vincenty_inverse(a, b)
            @test inv.distance > 0.0
        end
        @testset "due west azimuth" begin
            az = forward_azimuth(LLA(0.0, 0.0), LLA(-1.0, 0.0))
            @test az ≈ 270.0 atol=0.01
        end
        @testset "due south azimuth" begin
            az = forward_azimuth(LLA(0.0, 0.0), LLA(0.0, -1.0))
            @test az ≈ 180.0 atol=0.01
        end
        @testset "with non-default datum" begin
            a = LLA(-74.006, 40.7128)
            b = LLA(-0.1278, 51.5074)
            d_wgs84 = distance(a, b)
            d_grs80 = distance(a, b; datum=GRS80)
            @test d_grs80 ≈ d_wgs84 rtol=1e-6
        end
        @testset "accepts UTM" begin
            a = LLA(-74.006, 40.7128)
            b = LLA(-73.935, 40.730)
            a_utm = UTM(a, WGS84)
            d_lla = distance(a, b)
            d_utm = distance(a_utm, b)
            @test d_utm ≈ d_lla atol=1.0
        end
    end

    @testset "forward_azimuth" begin
        @testset "due east at equator" begin
            az = forward_azimuth(LLA(0.0, 0.0), LLA(1.0, 0.0))
            @test az ≈ 90.0 atol=0.01
        end
        @testset "due north" begin
            az = forward_azimuth(LLA(0.0, 0.0), LLA(0.0, 1.0))
            @test az ≈ 0.0 atol=0.01
        end
        @testset "bearing alias" begin
            a = LLA(0.0, 0.0)
            b = LLA(1.0, 1.0)
            @test bearing(a, b) == forward_azimuth(a, b)
        end
    end

    @testset "forward_geodesic" begin
        @testset "due east at equator" begin
            dest = forward_geodesic(LLA(0.0, 0.0), 90.0, 111_320.0)
            @test dest.lon ≈ 1.0 atol=0.01
            @test abs(dest.lat) < 0.01
        end
        @testset "due north at equator" begin
            dest = forward_geodesic(LLA(0.0, 0.0), 0.0, 110_574.0)
            @test abs(dest.lon) < 0.01
            @test dest.lat ≈ 1.0 atol=0.01
        end
        @testset "round-trip with vincenty_inverse" begin
            a = LLA(-74.006, 40.7128)
            b = LLA(-0.1278, 51.5074)
            inv = vincenty_inverse(a, b)
            dest = forward_geodesic(a, inv.azimuth1, inv.distance)
            @test dest.lon ≈ b.lon atol=1e-6
            @test dest.lat ≈ b.lat atol=1e-6
        end
        @testset "accepts AbstractFixedCoords" begin
            ecef = ECEF(LLA(0.0, 0.0), WGS84)
            dest = forward_geodesic(ecef, 0.0, 110_574.0)
            @test dest.lat ≈ 1.0 atol=0.01
        end
        @testset "due south at equator" begin
            dest = forward_geodesic(LLA(0.0, 0.0), 180.0, 110_574.0)
            @test abs(dest.lon) < 0.01
            @test dest.lat ≈ -1.0 atol=0.01
        end
        @testset "due west at equator" begin
            dest = forward_geodesic(LLA(0.0, 0.0), 270.0, 111_320.0)
            @test dest.lon ≈ -1.0 atol=0.01
            @test abs(dest.lat) < 0.01
        end
        @testset "zero distance returns origin" begin
            origin = LLA(-74.006, 40.7128)
            dest = forward_geodesic(origin, 45.0, 0.0)
            @test dest.lon ≈ origin.lon atol=1e-10
            @test dest.lat ≈ origin.lat atol=1e-10
        end
        @testset "long distance (NYC to London)" begin
            # NYC to London is ~5570 km at azimuth ~51°
            dest = forward_geodesic(LLA(-74.006, 40.7128), 51.2, 5.570e6)
            @test dest.lon ≈ -0.1278 atol=1.0
            @test dest.lat ≈ 51.5074 atol=1.0
        end
        @testset "with non-default datum" begin
            dest_wgs = forward_geodesic(LLA(0.0, 0.0), 90.0, 111_320.0; datum=WGS84)
            dest_grs = forward_geodesic(LLA(0.0, 0.0), 90.0, 111_320.0; datum=GRS80)
            @test dest_wgs.lon ≈ dest_grs.lon atol=1e-6
        end
        @testset "accepts UTM" begin
            origin_lla = LLA(0.0, 0.0)
            origin_utm = UTM(origin_lla, WGS84)
            dest_from_lla = forward_geodesic(origin_lla, 90.0, 111_320.0)
            dest_from_utm = forward_geodesic(origin_utm, 90.0, 111_320.0)
            @test dest_from_utm.lon ≈ dest_from_lla.lon atol=1e-6
        end
    end

    @testset "intermediate_point / midpoint" begin
        @testset "endpoints" begin
            a = LLA(0.0, 0.0)
            b = LLA(10.0, 0.0)
            p0 = intermediate_point(a, b, 0.0)
            p1 = intermediate_point(a, b, 1.0)
            @test p0.lon ≈ a.lon atol=1e-6
            @test p0.lat ≈ a.lat atol=1e-6
            @test p1.lon ≈ b.lon atol=1e-6
            @test p1.lat ≈ b.lat atol=1e-6
        end
        @testset "midpoint on equator" begin
            a = LLA(0.0, 0.0)
            b = LLA(10.0, 0.0)
            m = midpoint(a, b)
            @test m.lon ≈ 5.0 atol=0.01
            @test abs(m.lat) < 0.01
        end
        @testset "midpoint equidistant" begin
            a = LLA(-74.006, 40.7128)
            b = LLA(-0.1278, 51.5074)
            m = midpoint(a, b)
            d_am = distance(a, m)
            d_mb = distance(m, b)
            @test d_am ≈ d_mb rtol=0.001
        end
        @testset "quarter and three-quarter points" begin
            a = LLA(0.0, 0.0)
            b = LLA(10.0, 0.0)
            total = distance(a, b)
            q1 = intermediate_point(a, b, 0.25)
            q3 = intermediate_point(a, b, 0.75)
            @test distance(a, q1) ≈ total * 0.25 rtol=0.01
            @test distance(a, q3) ≈ total * 0.75 rtol=0.01
        end
        @testset "accepts AbstractFixedCoords" begin
            a_ecef = ECEF(LLA(0.0, 0.0), WGS84)
            b = LLA(10.0, 0.0)
            m = midpoint(a_ecef, b)
            @test m.lon ≈ 5.0 atol=0.01
        end
        @testset "midpoint with non-default datum" begin
            a = LLA(0.0, 0.0)
            b = LLA(10.0, 0.0)
            m_wgs = midpoint(a, b)
            m_grs = midpoint(a, b; datum=GRS80)
            @test m_wgs.lon ≈ m_grs.lon atol=1e-6
        end
    end

    @testset "BoundingBox" begin
        @testset "construction and show" begin
            bb = BoundingBox(-75.0, 40.0, -73.0, 41.0)
            @test bb.min_lon == -75.0
            @test bb.max_lat == 41.0
            s = sprint(show, bb)
            @test contains(s, "-75.0")
        end
        @testset "containment" begin
            bb = BoundingBox(-75.0, 40.0, -73.0, 41.0)
            @test in(LLA(-74.0, 40.5), bb)
            @test !in(LLA(-72.0, 40.5), bb)
            @test !in(LLA(-74.0, 42.0), bb)
        end
        @testset "containment with AbstractFixedCoords" begin
            bb = BoundingBox(-75.0, 40.0, -73.0, 41.0)
            ecef = ECEF(LLA(-74.0, 40.5), WGS84)
            @test in(ecef, bb)
        end
        @testset "boundingbox from points" begin
            pts = [LLA(-74.0, 40.0), LLA(-73.0, 41.0), LLA(-75.0, 40.5)]
            bb = boundingbox(pts)
            @test bb.min_lon == -75.0
            @test bb.max_lon == -73.0
            @test bb.min_lat == 40.0
            @test bb.max_lat == 41.0
        end
        @testset "antimeridian wrapping containment" begin
            # Box wrapping around 180°: min_lon > max_lon
            bb = BoundingBox(170.0, -10.0, -170.0, 10.0)
            @test in(LLA(175.0, 0.0), bb)    # east of antimeridian
            @test in(LLA(-175.0, 0.0), bb)   # west of antimeridian
            @test !in(LLA(0.0, 0.0), bb)     # far from antimeridian
        end
        @testset "boundary containment" begin
            bb = BoundingBox(-75.0, 40.0, -73.0, 41.0)
            @test in(LLA(-75.0, 40.0), bb)   # min corner
            @test in(LLA(-73.0, 41.0), bb)   # max corner
            @test in(LLA(-74.0, 40.0), bb)   # on min_lat edge
        end
        @testset "destination_boundingbox" begin
            bb = destination_boundingbox(LLA(0.0, 0.0), 100_000.0)
            @test bb.min_lon < 0.0
            @test bb.max_lon > 0.0
            @test bb.min_lat < 0.0
            @test bb.max_lat > 0.0
            # ~100km ≈ ~0.9 degrees at equator
            @test bb.max_lon ≈ 0.9 atol=0.1
        end
        @testset "destination_boundingbox at high latitude" begin
            bb = destination_boundingbox(LLA(0.0, 60.0), 100_000.0)
            @test bb.min_lon < 0.0
            @test bb.max_lon > 0.0
            # At 60° lat, longitude spread should be larger than at equator
            @test (bb.max_lon - bb.min_lon) > 1.0
        end
        @testset "destination_boundingbox with ECEF" begin
            bb_lla = destination_boundingbox(LLA(0.0, 0.0), 100_000.0)
            bb_ecef = destination_boundingbox(ECEF(LLA(0.0, 0.0), WGS84), 100_000.0)
            @test bb_ecef.min_lon ≈ bb_lla.min_lon atol=1e-6
            @test bb_ecef.max_lat ≈ bb_lla.max_lat atol=1e-6
        end
        @testset "boundingbox with mixed AbstractFixedCoords" begin
            pts = [LLA(-74.0, 40.0), ECEF(LLA(-73.0, 41.0), WGS84), LLA(-75.0, 40.5)]
            bb = boundingbox(pts)
            @test bb.min_lon ≈ -75.0 atol=1e-6
            @test bb.max_lat ≈ 41.0 atol=1e-6
        end
    end

    @testset "polygon_area" begin
        @testset "1°×1° square at equator" begin
            pts = [LLA(0, 0), LLA(1, 0), LLA(1, 1), LLA(0, 1)]
            area = polygon_area(pts)
            # ~12,300 km² for 1°×1° at equator
            @test area / 1e6 ≈ 12308.0 rtol=0.01
        end
        @testset "tiny triangle" begin
            pts = [LLA(0, 0), LLA(0.001, 0), LLA(0, 0.001)]
            area = polygon_area(pts)
            @test area > 0.0
            @test area < 1e7  # should be ~6150 m²
        end
        @testset "accepts AbstractFixedCoords" begin
            pts_lla = [LLA(0, 0), LLA(1, 0), LLA(1, 1), LLA(0, 1)]
            pts_ecef = [ECEF(p, WGS84) for p in pts_lla]
            @test polygon_area(pts_ecef) ≈ polygon_area(pts_lla) rtol=0.001
        end
        @testset "fewer than 3 points returns 0" begin
            @test polygon_area(LLA{Float64}[]) == 0.0
            @test polygon_area([LLA(0, 0)]) == 0.0
            @test polygon_area([LLA(0, 0), LLA(1, 0)]) == 0.0
        end
        @testset "reversed winding order gives same area" begin
            pts = [LLA(0, 0), LLA(1, 0), LLA(1, 1), LLA(0, 1)]
            pts_rev = reverse(pts)
            @test polygon_area(pts) ≈ polygon_area(pts_rev) rtol=0.001
        end
        @testset "larger polygon (10°×10° at equator)" begin
            pts = [LLA(0, 0), LLA(10, 0), LLA(10, 10), LLA(0, 10)]
            area = polygon_area(pts)
            # ~1,233,000 km² (roughly)
            @test area / 1e6 > 1.0e6  # > 1,000,000 km²
            @test area / 1e6 < 1.5e6  # < 1,500,000 km²
        end
        @testset "1°×1° square at high latitude (smaller)" begin
            pts_eq = [LLA(0, 0), LLA(1, 0), LLA(1, 1), LLA(0, 1)]
            pts_hi = [LLA(0, 60), LLA(1, 60), LLA(1, 61), LLA(0, 61)]
            area_eq = polygon_area(pts_eq)
            area_hi = polygon_area(pts_hi)
            @test area_hi < area_eq  # Higher latitude → smaller area
        end
        @testset "with non-default datum" begin
            pts = [LLA(0, 0), LLA(1, 0), LLA(1, 1), LLA(0, 1)]
            area_wgs = polygon_area(pts)
            area_grs = polygon_area(pts; datum=GRS80)
            @test area_wgs ≈ area_grs rtol=1e-6
        end
    end
end
