using Aqua
using StandardHaloModelsPP
using Test

@testset "Aqua" begin
    Aqua.test_all(StandardHaloModelsPP)
end
