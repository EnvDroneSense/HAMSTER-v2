using Test
using HAMSTER

@testset "HAMSTER" begin
    @test HAMSTER.version() isa VersionNumber
end
