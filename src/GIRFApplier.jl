## Mutable struct definition for GirfApplier type: following convention of JV in MATLAB repo
# composite type of GirfEssential and a field for gamma, the gyromagnetic ratio
mutable struct GirfApplier

    # PROPERTIES (fields)

    # GirfEssential type containing majority of GIRF data
    essential::GirfEssential

    # Gamma (gyromagnetic ratio)
    gamma::Float64

end

"""
    apply_girf(g::GirfApplier, gradient_in, time_in, time_out, direction)

    Function for predicting the gradient by applying the GIRF
    This is done in one dimension only. Repeat this function for multiple dimensions
# Arguments
* `g::GirfApplier` - GirfApplier containing GIRF data
* `gradient_in` - input gradient waveform corresponding to time_in vector
* `time_in` - time vector corresponding to input gradient waveform
* `time_out` - output time vector desired (decimation)
* `direction` - direction of the GIRF to use (x,y,z etc...)

# Outputs:
* `grad_OUT_resampled` - Vector, gradient waveform with GIRF applied, in the time points of `time_out`.
"""
function apply_girf(g::GirfApplier, gradient_in, time_in, time_out, direction)

    # GET GIRF DATA FROM THE GIRF DATATYPE
    freq_GIRF = g.essential.freq
    GIRF = g.essential.girf[:, direction]

    # check dimensions
    ns_GIRF = length(GIRF)
    df_GIRF = abs.(freq_GIRF[2] - freq_GIRF[1])

    # make copies so that original variables and objects are not mutated
    gradient_in_original = deepcopy(gradient_in)
    time_in_original = deepcopy(time_in)

    t_or1, t_or2, t_or3 = time2freq(time_in_original)

    # Set length of prediction
    T_GIRF = 1 ./ df_GIRF # Signal Length for nyquist fs = 2B and T = 1/2B
    T_O = max(time_in[end], time_out[end]) - min(time_in[1], time_out[1])
    T_ZF = min(T_GIRF, T_O) # Zero pad to the shortest time of the two

    # prepare input in time
    dt_in = abs.(time_in[2] - time_in[1])
    nZF = Int(ceil(T_ZF ./ dt_in)) # test this out

    # Do Zero filling
    gradient_in = vcat(zeros(nZF), gradient_in, zeros(nZF))
    time_in = vcat(time_in[1] .+ dt_in * (-nZF:-1), time_in, time_in[end] .+ dt_in * (1:nZF))

    # get updated frequency vector corresponding to Zero padded time vector
    f_in, df, f_max = time2freq(time_in)

    # Get frequency domain gradient
    IN = fftshift(fft(gradient_in)) .* dt_in
    ns_in = length(time_in)

    # interpolate GIRF onto the input grid
    GIRF_real_spline = Spline1D(freq_GIRF, real.(GIRF), w = ones(length(freq_GIRF)), k = 2, bc = "zero")
    GIRF_imaginary_spline = Spline1D(freq_GIRF, imag.(GIRF), w = ones(length(freq_GIRF)), k = 2, bc = "zero")

    # recombine interpolated reals and imaginary parts
    GIRF_ip = GIRF_real_spline(f_in) .+ 1im .* GIRF_imaginary_spline(f_in)

    ## filter the output to avoid aliasing
    l_girf_ip = length(GIRF_ip)
    npad = Int(floor(l_girf_ip * 0.7))

    windowFunction = fftshift(tukey(l_girf_ip - npad, 0.1; padding = npad, zerophase = true))
    GIRF_ip = windowFunction .* GIRF_ip

    # Convolution in Time (mult in freq) to get output Gradient train in frequency domain
    OUT = GIRF_ip .* IN

    grad_OUT = real.(ifft(ifftshift(OUT)) ./ dt_in)

    grad_OUT_spline = Spline1D(time_in, grad_OUT, w = ones(length(time_in)), k = 2, bc = "zero")
    grad_OUT_resampled = grad_OUT_spline(time_out)

    return grad_OUT_resampled

end