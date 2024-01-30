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

    nodes = gradients_to_nodes(ones(2,1) .* correctedGradients')
    gradients = nodes_to_gradients(nodes)

    @test sqrt.(sum(abs2,ones(2,1).*correctedGradients' - gradients)) < 1e-9

    nodes = gradients_to_nodes(ones(3,1) .* correctedGradients')
    gradients = nodes_to_gradients(nodes)
    
    @test sqrt.(sum(abs2,ones(3,1).*correctedGradients' - gradients)) < 1e-9

end

function test_display_girf()

    girf_applier_k1 = load_object("data/girf_app_k1.jld2")
    displayGirf(girf_applier_k1.essential)
    @test true

end

function test_girf_applier()

    girf_applier_k1 = load_object("data/girf_app_k1.jld2")

    girf_applier_2 = GirfApplier(deepcopy(girf_applier_k1.essential),deepcopy(girf_applier_k1.gamma))

    @test girf_applier_2.gamma == girf_applier_k1.gamma
    @test girf_applier_2.essential.df == girf_applier_k1.essential.df

end

function test_package()

    @testset "Package" begin

        test_grads()
        test_girf_applier()
        test_display_girf()

    end

end

test_package()