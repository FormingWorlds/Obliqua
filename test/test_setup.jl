using Test
using Obliqua

@testset "Obliqua tidal model import check" begin

    @testset "solid-tide models" begin
        solid_funcs = [
            :run_solid0d, 
            :run_solid1d, 
            :run_solid1d_relax, 
            :run_solid1d_mush, 
            :run_solid1d_mush_relax
        ]
        
        for func in solid_funcs
            @test isdefined(Obliqua, func)
        end
    end

    @testset "fluid-tide models" begin
        fluid_funcs = [
            :run_fluid0d, 
            :run_fluid1d
        ]
        
        for func in fluid_funcs
            @test isdefined(Obliqua, func)
        end
    end

    @testset "mushy-tide models" begin
        mushy_funcs = [
            :run_interp
        ]
        
        for func in mushy_funcs
            @test isdefined(Obliqua, func)
        end
    end
end