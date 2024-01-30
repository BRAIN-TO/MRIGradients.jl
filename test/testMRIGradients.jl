function test_grads()

    # Placeholder for now
    @test true

end

function test_package()

    @testset "Package" begin
        
        test_grads()

    end

end

test_package()
