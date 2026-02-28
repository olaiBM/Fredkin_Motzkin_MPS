using ITensors, LinearAlgebra
using ITensorMPS
using LsqFit

function dense_mps(qn_mps::MPS)
    """
    Returns a dense version of the MPS
    """
     dense_tensors = [dense(A) for A in qn_mps]
     return MPS(dense_tensors)
end


function expect_Sz_analytic_F(r::Int)
    """
    The analytic expression of the expectation value of Sz for the Fredkin chain
    """
    tot = 0
    if r == 1
         tot = 0.5
    else      
         part_1 = 1/(sqrt(2*pi)*sqrt(r-1))
         part_2 = (-1)^(r+1)/(4*sqrt(pi) * (2 * (r-1))^(1.5))
         tot = part_1 + part_2
    end
   return tot
end



function expect_Sz_analytic_M(r::Int)
    """
    The analytic expression of the expectation value of Sz for the Motzkin chain
    """
    tot = 0
    if r == 1
     ### This was calculated from the intergals in Appendix C at r = 1
         K = 8*sqrt(6)/9
         I_y_and_const = K * (1/(12 * (2/3)^(3/2) * pi)) 
         tot = 2 * pi * I_y_and_const
    else     
         ### This was calculated by doing the integrals in Appendix C, keeping all terms not exponentially small in N
         K = 8*sqrt(6)/9
         I_1_contribution = ((3*K)/(4*sqrt(2*pi)*sqrt(r-1)))   #

         I_3_y_and_const = K * (1/(12 * (2/3)^(3/2) * pi)) 
         I_3_x_pi2 = (0.5 * (1/3)^(r-1) * ((1/(r+1))*(1 - (1-pi)^(r+1)) + (1/(r))*(1 - (1-pi)^(r))))
         I_3_x_pi = (-1/3)^(r-1) * (1/(r-1)^(3/2)) * sqrt(pi/2)
         tot = I_1_contribution + I_3_y_and_const * (I_3_x_pi2 + I_3_x_pi)
    end
   return tot       
end




function compute_decay_envelope(corr_values::Vector)
     """
     Computes decay behaviour by setting the value of corr_values 
     at r = i+0.5 as the average of corr_values[i] and corr_values[i+1]
     """
     n = length(corr_values)
 
     # Initialize output arrays
     envelope = Float64[]
     r_vals = Float64[]
 
     # First value is fixed (boundary)
     push!(envelope, corr_values[1])
     push!(r_vals, 1.0)
 
     # Average zig-zag pairs and assign them r = i + 0.5
     for i in 2:2:(n-2)
         avg = (corr_values[i] + corr_values[i+1]) / 2
         push!(envelope, avg)
         push!(r_vals, i + 0.5)
     end
 
     return r_vals, envelope
 end


### Define unzip function 
unzip(pairs) = map(x -> getindex.(pairs, x), (1, 2))

function fit_line(x, y) # Utility function for least squares fitting of y = ax + b
    A = hcat(x, ones(length(x)))
    return A \ y
end



### Choose which values to use in line fit (for N = 2000) for the various qs (found by plotting and testing)

F_idx_qs = Dict{Float64, Tuple{Int, Int}}()
F_idx_qs[0.9995] = (10, 400)
F_idx_qs[0.999] = (10, 400)
F_idx_qs[0.9985] = (10, 400)
F_idx_qs[0.998] = (10, 400)
F_idx_qs[0.9975] = (10, 350)
F_idx_qs[0.997] = (10, 300)
F_idx_qs[0.9965] = (10, 275)
F_idx_qs[0.996] = (10, 275)
F_idx_qs[0.9955] = (10, 250)
F_idx_qs[0.995] = (10, 225)
F_idx_qs[0.994] = (10, 225)
F_idx_qs[0.993] = (10, 200)
F_idx_qs[0.992] = (10, 150)
F_idx_qs[0.991] = (10, 125)
F_idx_qs[0.99] = (10, 125)

M_idx_qs = Dict{Float64, Tuple{Int, Int}}()
M_idx_qs[0.9995] = (10, 400)
M_idx_qs[0.999] = (10, 400)
M_idx_qs[0.9985] = (10, 300)
M_idx_qs[0.998] = (10, 250)
M_idx_qs[0.9975] = (10, 225)
M_idx_qs[0.997] = (10, 225)
M_idx_qs[0.9965] = (10, 200)
M_idx_qs[0.996] = (10, 175)
M_idx_qs[0.9955] = (10, 125)
M_idx_qs[0.995] = (10, 125)
M_idx_qs[0.994] = (10, 125)
M_idx_qs[0.993] = (10, 110)
M_idx_qs[0.992] = (10, 110)
M_idx_qs[0.991] = (10, 110)
M_idx_qs[0.99] = (10, 110)


function compute_correlation_length_MPS(N::Int, q_vals::Vector{Float64}, type::String)
    """
    Computes the correlation lengths xi from the MPS expectation values
    """
    # slope_vals = Float64[]
    xis = Float64[]
    slope_vals = Float64[]
    idx_qs = Dict{Float64, Tuple{Int, Int}}()
    if type == "Fredkin"
        idx_qs = F_idx_qs
    elseif type == "Motzkin"
        idx_qs = M_idx_qs
    else
        error("Type can only be Fredkin or Motzkin!")
    end

    # Loop over q values
    for q in q_vals
        # Load data
        start_idx_q,end_idx_q = idx_qs[q]
        prefix = type * "_Sz"
        @load datadir("sims", savename(prefix, (N=N, q=q), "jld2"; digits = 4)) Sz_expect_data
        filtered_data = []
        # Process data

        r_vals_decay, avg_exp_data= compute_decay_envelope(Sz_expect_data)
        filtered_data = [(r, z) for (r, z) in zip(r_vals_decay, avg_exp_data) if z > 0]

        filtered_r, filtered_z = unzip(filtered_data)
        log_z = log10.(filtered_z)

        # Slice
        x_slice = filtered_r[start_idx_q:end_idx_q]
        y_slice = log_z[start_idx_q:end_idx_q]

        # Fit line
        coeffs = fit_line(x_slice, y_slice)
        slope, _ = coeffs
        xi = -(1/slope)
        # Store
        push!(slope_vals, slope)
        push!(xis, xi)
    end
    return xis, slope_vals
end



function load_and_filter_ql1_data(N::Int, q_vals::Vector{Float64}, type::String)
    """
    Loads and filters data used to to line fits for the q<1 case (for semilog plot)
    """
    q_filtered_data_dict = Dict{Float64, Tuple}()
    for q in q_vals
        prefix = type * "_Sz"
        @load datadir("sims", savename(prefix, (N=N, q=q), "jld2"; digits = 4)) Sz_expect_data
        r_vals, envelope = compute_decay_envelope(Sz_expect_data) # Average between even and odd sites (raw data is rugged)
        filtered_data = [(r, z) for (r, z) in zip(r_vals, envelope) if z > 0]
        r_filtered, expectation_filtered = unzip(filtered_data)
        q_filtered_data_dict[q] = (r_filtered, expectation_filtered)
    end
    return q_filtered_data_dict
end



function return_slope_and_intercept_ql1(q_filtered_data_dict::Dict, q_start_and_stop_index::Dict)
    """
    Returns the intercept and slope for line fits in the q<1 cases (line fit in the semilog plot)
    """
    q_intercept_slope_dict = Dict{Float64, Vector}()
    for q in keys(q_filtered_data_dict)
        log_data = log.(q_filtered_data_dict[q][2])
        r_vals = q_filtered_data_dict[q][1]

        # Determine window for line fit
        log_data = log_data[q_start_and_stop_index[q][1]:q_start_and_stop_index[q][2]] 
        r_vals = r_vals[q_start_and_stop_index[q][1]:q_start_and_stop_index[q][2]]
        q_intercept_slope_dict[q] = fit_line(r_vals, log_data)
    end

    return q_intercept_slope_dict
end

function linear_fit(r::Int, slope::Float64, intercept::Float64)
    return intercept + r*slope
end


### Fermi - Dirac models used for curce fitting
function fermi_dirac_model_F(x, p)
    xi = p[1]
    r0 = N / 2
    return (1.0 ./ (exp.((x .- r0) ./ xi) .+ 1.0)).- 0.5
end

function fermi_dirac_model_M(x, p)
    xi = p[1]
    r0 = N / 2
    return (2.0 ./ (exp.((x .- r0) ./ xi) .+ 1.0)).- 1.0
end


function fermi_dirac_fit(r_vals::Vector, qg1_data_dict::Dict, type::String)
    initial_guess = [10.0]
    q_fit_dict = Dict{Float64, Float64}()
    for q in keys(qg1_data_dict)
        if type == "Fredkin"
            q_fit_dict[q] = curve_fit(fermi_dirac_model_F, r_vals, qg1_data_dict[q], initial_guess).param[1]
        elseif type == "Motzkin"
            q_fit_dict[q] = curve_fit(fermi_dirac_model_M, r_vals, qg1_data_dict[q], initial_guess).param[1]
        else
            error("Type can only be Fredkin or Motzkin!")
        end

    end
    return q_fit_dict
end

function load_q_data(N::Int, q_vals::Vector{Float64}, type::String)
    """
    Loads and filters data used to to line fits for the q<1 case (for semilog plot)
    """
    q_data_dict = Dict{Float64, Vector}()
    for q in q_vals
        prefix = type * "_Sz"
        @load datadir("sims", savename(prefix, (N=N, q=q), "jld2"; digits = 4)) Sz_expect_data
        q_data_dict[q] = Sz_expect_data
    end
    return q_data_dict
end
