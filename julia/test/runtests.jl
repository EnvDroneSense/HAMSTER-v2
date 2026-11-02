using Test
using HAMSTER

@testset "HAMSTER" begin
    @test HAMSTER.version() isa VersionNumber
    include("test_utilities.jl")
end
