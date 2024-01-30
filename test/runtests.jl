using MRIGradients
using Test
using Scratch
using Statistics
using JLD2

global const tmpdir = @get_scratch!("temp")
@info "For the output of tests, please refer to the directory $tmpdir"

include("testMRIGradients.jl")