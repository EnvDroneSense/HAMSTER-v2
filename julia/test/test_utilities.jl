const REFDIR = joinpath(@__DIR__, "reference")

# Reference data: manifest.txt lists "name dtype dims" in numpy order;
# reading the bytes with reversed dims keeps the memory layout.
function load_reference()
    types = Dict("float32" => Float32, "float64" => Float64, "int64" => Int64)
    ref = Dict{String,Array}()
    for line in eachline(joinpath(REFDIR, "manifest.txt"))
        isempty(strip(line)) && continue
        name, dtype, dims... = split(line)
        a = Array{types[dtype]}(undef, reverse(parse.(Int, dims))...)
        read!(joinpath(REFDIR, name * ".bin"), a)
        ref[name] = a
    end
    return ref
end

@testset "utilities" begin
    ref = load_reference()
    nx, ny, mass = Int(ref["grid_dims"][1]), Int(ref["grid_dims"][2]), ref["grid_dims"][3]
    xi, yi = ref["grid_x"] .+ 1, ref["grid_y"] .+ 1

    lon, lat = midpoints(ref["lona"], ref["lonb"], ref["lata"], ref["latb"])
    dlon = abs.(lon .- ref["midp_lon"])
    @test maximum(min.(dlon, 360 .- dlon)) < 0.01
    @test maximum(abs.(lat .- ref["midp_lat"])) < 0.01

    @test only_theta(ref["qv"], ref["temp"], ref["rho"]) ≈ ref["theta"] rtol = 1e-6

    @test ongrid(nx, ny, xi, yi, ref["ongrid_f64"], mass) == ref["ongrid_f64_out"]
    @test ongrid1d(nx, ny, xi[:, 1], yi[:, 1], ref["ongrid1d_f64"], mass) ==
          (ref["ongrid1d_f64_out1"], ref["ongrid1d_f64_out2"])

    # parcels outside the grid are skipped
    @test sum(ongrid1d(2, 2, [1, 99, -5], [1, 1, 1], [1.0, 1.0, 1.0], 1.0)[1]) == 1
end
