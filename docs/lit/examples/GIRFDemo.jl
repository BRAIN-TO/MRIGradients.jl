## GIRF Demo file showing basic loading of GIRF

using Waveforms, MRIGradients

## Demonstration:

## Load saved GIRF into variable of type GIRFEssential by calling the loadGirf function which takes the degree of the GIRF to load as the argument

# 0th degree (k0)
girf_k0 = loadGirf(0, 1)
setIdentifier!(girf_k0, "k0 GIRF")

# first degree (k1)
girf_k1 = loadGirf(1, 1)
setIdentifier!(girf_k1, "k1 GIRF")

## Create GirfApplier objects

girfApplier_k0 = GirfApplier(girf_k0, 42577478)
girfApplier_k1 = GirfApplier(girf_k1, 42577478)

## Generate Test Data

time_vector = collect(0:1000) ./ 100
testSinusoid = map(x -> trianglewave1(x), time_vector) # options are squarewave1 or sawtoothwave1 as well

# correct test data to increase frequency artificially, and also set the first gradient sample to start at 0
testSinusoid[1] = 0
time_vector = time_vector ./ 1000

## Filter test data with GIRF
correctedGradients = apply_girf(girfApplier_k1, testSinusoid, time_vector, time_vector, 1)