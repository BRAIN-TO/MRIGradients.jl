function test_grads()

    # Placeholder for now
    girf_applier_k1 = load_object("data/girf_app_k1.jld2")

    time_vector = collect(0:1000) ./ 100
    testSinusoid = map(x -> trianglewave1(x), time_vector) # options are squarewave1 or sawtoothwave1 as well

    # correct test data to increase frequency artificially, and also set the first gradient sample to start at 0
    testSinusoid[1] = 0
    time_vector = time_vector ./ 1000

    ## Filter test data with GIRF
    correctedGradients = apply_girf(girf_applier_k1, testSinusoid, time_vector, time_vector, 1)

    @test abs.(correctedGradients[end] - -0.0258168) < 1e-6

end

function test_package()

    @testset "Package" begin

        test_grads()

    end

end

test_package()