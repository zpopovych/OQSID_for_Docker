include("LiPoSID.jl")

using QuantumOptics
basis = NLevelBasis(2)
using LinearAlgebra

using HDF5
using DynamicPolynomials

using Dates

using Statistics

using TSSOS

function lindblad_rhs(ρ, H, J::Array)
    """
    Right hand side of the Lindblad master equation with multiple disipators
    """
   
    Σ = sum([ ( Jⱼ * ρ * Jⱼ' - (Jⱼ' * Jⱼ  * ρ + ρ * Jⱼ' * Jⱼ)/2 ) for Jⱼ in J ])
    
    return -im * (H * ρ - ρ * H) + Σ 
    
end

# Define polynomial variables
@polyvar ϵ α β r p

Hˢʸᵐᵇ = [ ϵ            α + im* β
          α - im* β   -ϵ         ]/2

J₁ˢʸᵐᵇ = [ 0   r
           0   0. ]

J₂ˢʸᵐᵇ = p * [ 0    1
               1    0. ]

J₃ˢʸᵐᵇ = p * [ 0   -im
               im   0  ]

J₄ˢʸᵐᵇ = p * [ 1    0
               0   -1. ]


Jˢʸᵐᵇ = [ J₁ˢʸᵐᵇ, J₂ˢʸᵐᵇ, J₃ˢʸᵐᵇ, J₄ˢʸᵐᵇ]

function simpson_obj(ρ, t, H, J)
    
    obj = 0
    for i in 3:length(ρ)
        obj += LiPoSID.frobenius_norm2(
            ρ[i] - ρ[i-2] - (t[i]-t[i-1])lindblad_rhs((ρ[i-2] + 4ρ[i-1] + ρ[i])/3, H, J)
        )
    end
    obj = sum(real(coef) * mon for (coef, mon) in zip(coefficients(obj), monomials(obj)))
    return obj
end

function lindblad_GEXY_obj(ρᵍᵉˣʸ, tᵍᵉˣʸ, Hˢʸᵐᵇ, Jˢʸᵐᵇ)

    ρᵍ, ρᵉ, ρˣ, ρʸ = ρᵍᵉˣʸ
    tᵍ, tᵉ, tˣ, tʸ = tᵍᵉˣʸ

    polyG = simpson_obj(ρᵍ, tᵍ, Hˢʸᵐᵇ, Jˢʸᵐᵇ)
    polyE = simpson_obj(ρᵉ, tᵉ, Hˢʸᵐᵇ, Jˢʸᵐᵇ)
    polyX = simpson_obj(ρˣ, tˣ, Hˢʸᵐᵇ, Jˢʸᵐᵇ)
    polyY = simpson_obj(ρʸ, tʸ, Hˢʸᵐᵇ, Jˢʸᵐᵇ)

    polyGEXY = polyG + polyE + polyX + polyY

    return polyGEXY
end

constraints = [ϵ, p] #r, e, ϕ, p]

# Here we set parameters what data to use for training and testing

ρᵍ₀ = [ 1 0
        0 0 ]    # state to measure initial distance from

dodeca_10_states = ["D"*string(n) for n=1:10];
basis_states = ["B"*string(n) for n=1:4];

train_states = basis_states 
test_states = dodeca_10_states

all_states = vcat(train_states, test_states)

function min_cs_tssos(p, constrs)

    coeffs = coefficients(p)
    reg_coef = 0

    pop =[p+reg_coef*sum(variables(p).^2), constrs...] ./ maximum(abs.(coeffs))

    d = maxdegree(p)
    
    # Initial optimization step
    opt, sol, data = cs_tssos_first(pop, variables(pop), d; solution=true, QUIET=true)
    ref_sol, flag = TSSOS.refine_sol(opt, sol, data; QUIET=true)
    prev_opt, prev_sol, prev_data = opt, sol, data 

    # Check if the solution needs further refinement
    if flag != 0
        while ~isnothing(sol) && flag != 0
            prev_opt, prev_sol, prev_data = opt, sol, data
            opt, sol, data = cs_tssos_higher!(data; solution=true, QUIET=true) 
        end
        ref_sol, flag = TSSOS.refine_sol(prev_opt, prev_sol, prev_data; QUIET=true)
    end

    solution = variables(p) => ref_sol

    if flag == 0 
        status_name = "GLOBAL"
    else
        status_name = "LOCAL/FAIL"
    end

    return solution, status_name

end

println(" MULTI DURATION SYSTEM IDENTIFICATION w CONSTRAINED TSSOS and LINDBLAD 4ops Frobenius objective QO simulation")

γ = [ "0.079477",  "0.25133", "0.79477", "2.5133", "7.9477", "25.133", "79.477", "251.33"]

date_and_time_string = string(Dates.format(now(), "yyyy-u-dd_at_HH-MM"))

evol_data_file_name = "../data/ALL_GAMMAS_B4_D10.h5"

tests_dir = "../results/"

relative_threshold = 1e-15

rltrs = string(convert(Int, floor(log10(relative_threshold))))

tests_data_file_name = "MULT_LINDBLAD_CONSTR_TSSOS_threshold_1e"*rltrs*"_FROB_QO.h5"#*date_and_time_string * ".h5"

FminGammas = []
FmedianGammas = []
FmeanGammas = []
Epsilons = []
CoefRanges = []

for γᵢ in γ
    println("γ =  "*γᵢ)

    # Save results into HDF5
    h5open(tests_dir*tests_data_file_name,"cw") do fid
        γ_group = create_group(fid, γᵢ) # Create gamma coupling group
    end
    
    # Read time series data for the current γ
    tᵍᵉˣʸ, ρᵍᵉˣʸ = LiPoSID.read_GEXY_timeevolution(evol_data_file_name, γᵢ)
    tᵍ, tᵉ, tˣ, tʸ = tᵍᵉˣʸ
    
    # Find the end of the shortest time series
    t_end = minimum([tᵍ[end], tᵉ[end], tˣ[end], tʸ[end]])

    # Dynamically generate the duration steps
    max_steps = Int(floor(t_end * 2 * parse(Float64, γᵢ)))  # Maximum steps based on t_end and γᵢ
    duration_steps = [i/(2*parse(Float64, γᵢ)) for i in 1:max_steps]  # Generates dynamic steps up to t_end
    
    # Ensure steps do not exceed the shortest time series
    duration_steps = duration_steps[duration_steps .<= t_end]

    for t_duration in duration_steps
        println("Processing time duration: ", t_duration)

        elapsed_time = @timed begin
            # Masking the time series up to t_duration
            t_masked = [t[t .<= t_duration] for t in [tᵍ, tᵉ, tˣ, tʸ]]
            min_length = minimum(length.(t_masked))  # Ensure all time series are truncated to the same length
            t_masked = [t[1:min_length] for t in t_masked]

            ρ_masked = [ρ[1:min_length] for ρ in ρᵍᵉˣʸ] # Truncate ρᵍᵉˣʸ accordingly

            # Objective function construction and minimization
            polyGEXYfull = lindblad_GEXY_obj(ρ_masked, t_masked, Hˢʸᵐᵇ, Jˢʸᵐᵇ)
            polyGEXY = LiPoSID.filter_odd_terms_by_relative_threshold(polyGEXYfull, relative_threshold)
            sol, status = min_cs_tssos(polyGEXY, constraints)
        end

        @show minimum(abs.(coefficients(polyGEXY)))
        @show maximum(abs.(coefficients(polyGEXY)))
        push!(CoefRanges, LiPoSID.coefficient_range(polyGEXYfull))

        print(" status:", status)
        print(" runtime :", elapsed_time.time)

        Hˢⁱᵈ = convert.(ComplexF64, subs(Hˢʸᵐᵇ, sol))
        Jˢⁱᵈ = [convert.(ComplexF64,subs(Jᵢˢʸᵐᵇ, sol)) for Jᵢˢʸᵐᵇ in Jˢʸᵐᵇ]

        epsilon = subs(ϵ, sol)
        push!(Epsilons, epsilon)

        # Save results into HDF5
        h5open(tests_dir*tests_data_file_name,"cw") do fid
            γ_group = open_group(fid, γᵢ) # Create gamma coupling group
            t_group = create_group(γ_group, string(t_duration)) # Create time duration group
            t_group["epsilon"] = convert(Float64, epsilon)
            t_group["H"] = convert.(ComplexF64, Hˢⁱᵈ)
            t_group["J1"] = convert.(ComplexF64, Jˢⁱᵈ[1])
            t_group["J2"] = convert.(ComplexF64, Jˢⁱᵈ[2])
            t_group["J3"] = convert.(ComplexF64, Jˢⁱᵈ[3])
            t_group["J4"] = convert.(ComplexF64, Jˢⁱᵈ[4])
            t_group["status"] = status
            t_group["runtime"] = elapsed_time.time
        end

        println()

        FminStates = []
        FmedianStates = []
        FmeanStates = []

        print("Processing states: ")

        for state in test_states # Loop over initial states
            print(state*" ")
            start_time = time()

            tₛ, ρₛ = LiPoSID.read_timeevolution(evol_data_file_name, state, γᵢ)
            ρₛ = convert(Vector{Matrix{ComplexF64}}, ρₛ)
            ρᵗˢᵗ = [DenseOperator(basis, Hermitian(ρₜ)) for ρₜ in ρₛ]
            tᵗˢᵗ = convert(Vector{Float64}, tₛ)

            ρₒ = DenseOperator(basis, ρₛ[1])
            dt = tᵗˢᵗ[2] - tᵗˢᵗ[1]
            tᵉⁿᵈ = tᵗˢᵗ[end]

            effective_Lindblad = Jˢⁱᵈ
            effective_Lindblad_ops = [DenseOperator(basis, j) for j in effective_Lindblad]

            tout, ρ_t_kossak = timeevolution.master(convert.(Float64, tᵗˢᵗ), ρₒ, DenseOperator(basis, Hˢⁱᵈ), effective_Lindblad_ops)
            ρˢⁱᵈ = [ρₜ.data for ρₜ in ρ_t_kossak]

            F = LiPoSID.fidelity_series(basis, ρₛ, ρˢⁱᵈ)

            FminState = minimum(F)
            FmedianState = mean(F)
            FmeanState = mean(F)

            h5open(tests_dir*tests_data_file_name,"cw") do fid
                γ_group = open_group(fid, γᵢ)
                t_group = open_group(γ_group, string(t_duration))
                init_state_group = create_group(t_group, state)
                init_state_group["Fidelity"] = convert.(Float64, F)
                init_state_group["Fmin"] = convert(Float64, FminState)
                init_state_group["Fmedian"] = convert(Float64, FmedianState)
                init_state_group["Fmean"] = convert(Float64, FmeanState)

                init_state_group["tr_dist_grnd"] = LiPoSID.TrDist(ρₛ[1], ρᵍ₀)
                init_state_group["time"] = tᵗˢᵗ
            end

            push!(FminStates, FminState)
            push!(FmedianStates, FmedianState)
            push!(FmeanStates, FmeanState)
        end

        F_mean_value = mean(FmeanStates)
        F_median_value = median(FmedianStates)
        F_min_value = minimum(FminStates)

        push!(FminGammas, F_min_value)
        push!(FmedianGammas, F_median_value)
        push!(FmeanGammas, F_mean_value)

        println()

        println("Median fidelity for γ = "*γᵢ*" and time duration "*string(t_duration)*": ", F_median_value)

    end
end
    
    
println(tests_data_file_name)