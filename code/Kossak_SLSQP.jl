using JuMP
using NLopt
using LinearAlgebra

using HDF5

using QuantumOptics
basis = NLevelBasis(2)

include("LiPoSID.jl")

using Statistics

using Dates

function det2x2(m)
    real(m[1,1]*m[2,2] - m[1,2]*m[2,1])
end

function infidelity_norm(ρ, σ)
    abs(1 - real(tr(ρ * σ)) - 2*sqrt(abs(det(ρ)*det(σ))))#^2
end

function fidelity(ρ, σ)
    real(tr(ρ * σ)) + 2*sqrt(abs(det(ρ)*det(σ)))#^2
end

function frobenius_norm2(m)
    return abs(tr(m * m'))
end

# Define the Pauli matrices
σˣ = [0 1; 1 0]
σʸ = [0 im; -im 0]
σᶻ = [1 0; 0 -1]

# Define the basis elements
fᴷ₁ = σˣ / 2
fᴷ₂ = σʸ / 2
fᴷ₃ = σᶻ / 2

# Check orthogonality and normalization
@assert tr(fᴷ₁ * fᴷ₂) ≈ 0
@assert tr(fᴷ₁ * fᴷ₃) ≈ 0
@assert tr(fᴷ₂ * fᴷ₃) ≈ 0
@assert tr(fᴷ₁ * fᴷ₁) ≈ 1 / 2
@assert tr(fᴷ₂ * fᴷ₂) ≈ 1 / 2
@assert tr(fᴷ₃ * fᴷ₃) ≈ 1 / 2

fᴼᴺᴮ = [fᴷ₁, fᴷ₂, fᴷ₃]

# Function to calculate Dc
function DcKossak(ρ, t, H, C)
    U = (H * ρ - ρ * H) / im
    D = sum(C .* [2 * fᵢ * ρ * fⱼ' - ρ * fⱼ' * fᵢ - fⱼ' * fᵢ * ρ for fᵢ in fᴼᴺᴮ, fⱼ in fᴼᴺᴮ]) / 2
    return U + D
end

# Define functions to construct H and C from parameters
function construct_H(ϵ, h_Re, h_Im)
    return [
             ϵ                     h_Re + im * h_Im
             h_Re - im * h_Im     -ϵ
           ] / 2
end

function construct_C(γ, a)
    return [
        -γ[1] + γ[2] + γ[3]  -im * a[3]            im * a[2]
        im * a[3]            γ[1] - γ[2] + γ[3]   -im * a[1]
        -im * a[2]           im * a[1]             γ[1] + γ[2] - γ[3]
    ]
end

# Objective function

function kossak_obj(ϵ, h_Re, h_Im, γ, a, ρ, t) #Simpson 
    H = construct_H(ϵ, h_Re, h_Im)
    C = construct_C(γ, a)
    obj = 0.0
    for i in 3:length(ρ)
        ρ1 = ρ[i]
        ρ2 = ρ[i - 2] + 
        (t[i] - t[i - 1]) * 
        (DcKossak(ρ[i], t[i], H, C) +
         4 * DcKossak(ρ[i - 1], t[i - 1], H, C) 
         + DcKossak(ρ[i - 2], t[i - 2], H, C)) / 3
        #obj += infidelity_norm(ρ1,ρ2)
        #obj += fidelity(ρ1,ρ2)
        obj += frobenius_norm2(ρ1-ρ2)
    end
    return obj
end

# Define the objective function wrapper
function kossak_objectiveGEXY(ϵ, h_Re, h_Im, γ1, γ2, γ3, a1, a2, a3)

    ρᵍ, ρᵉ, ρˣ, ρʸ = ρᵍᵉˣʸ
    tᵍ, tᵉ, tˣ, tʸ = tᵍᵉˣʸ

    γ = [γ1, γ2, γ3]
    a = [a1, a2, a3]
    objGEXY = kossak_obj(ϵ, h_Re, h_Im, γ, a, ρᵍ, tᵍ) + 
              kossak_obj(ϵ, h_Re, h_Im, γ, a, ρᵉ, tᵉ) + 
              kossak_obj(ϵ, h_Re, h_Im, γ, a, ρˣ, tˣ) + 
              kossak_obj(ϵ, h_Re, h_Im, γ, a, ρʸ, tʸ)
    return objGEXY
end

# Define a function to create and solve the model with given initial guesses and a general upper bound
function solve_kossak_model(initial_guess, upper_bound)

    #upper_bound = 260
    
    model = Model(optimizer_with_attributes(NLopt.Optimizer, "algorithm" => :LD_SLSQP))

    # Register the element-wise power operator
    register(model, :.^, 2, (x, y) -> x.^y, autodiff=true)

    # Register the custom objective function
    register(model, :kossak_objectiveGEXY, 9, kossak_objectiveGEXY, autodiff = true)

    # Define lower and upper bounds
    lower_bound = -upper_bound

    # Define variables with general upper and lower bounds and initial guesses
    @variable(model, 0 <= ϵ <= upper_bound, start = initial_guess[1])
    @variable(model, lower_bound <= h_Re <= upper_bound, start = initial_guess[2])
    @variable(model, lower_bound <= h_Im <= upper_bound, start = initial_guess[3])
    @variable(model, lower_bound <= γ[i=1:3] <= upper_bound)
    set_start_value.(γ, initial_guess[4:6])
    @variable(model, lower_bound <= a[i=1:3] <= upper_bound)
    set_start_value.(a, initial_guess[7:9])

    # Define κ terms using anonymous construction with bounds
    @variable(model, lower_bound <= κ[1:3] <= upper_bound)

    @constraint(model, κ[1] == -γ[1] + γ[2] + γ[3])
    @constraint(model, κ[2] == γ[1] - γ[2] + γ[3])
    @constraint(model, κ[3] == γ[1] + γ[2] - γ[3])

    # Define the real part of Cˢʸᵐᵇ matrix using fewer constraints
    @variable(model, C_Re[1:3, 1:3], lower_bound = lower_bound, upper_bound = upper_bound)
    @constraint(model, [i=1:3], C_Re[i, i] == κ[i])
    @constraint(model, [i=1:3, j=1:3; i != j], C_Re[i, j] == 0)

    # Define the imaginary part of Cˢʸᵐᵇ matrix using fewer constraints
    @variable(model, C_Im[1:3, 1:3], lower_bound = lower_bound, upper_bound = upper_bound)
    @constraint(model, C_Im[1, 1] == 0)
    @constraint(model, C_Im[2, 2] == 0)
    @constraint(model, C_Im[3, 3] == 0)
    @constraint(model, C_Im[1, 2] == -a[3])
    @constraint(model, C_Im[1, 3] == a[2])
    @constraint(model, C_Im[2, 1] == a[3])
    @constraint(model, C_Im[2, 3] == -a[1])
    @constraint(model, C_Im[3, 1] == -a[2])
    @constraint(model, C_Im[3, 2] == a[1])

    # Define the κ-related constraints
    @constraint(model, sum(κ) >= 0)
    @constraint(model, [i=1:3], κ[i] >= 0)

    # Nonlinear constraints
    @NLconstraint(model, κ[1] * κ[2] + κ[3] * κ[1] + κ[2] * κ[3] - (a[1]^2 + a[2]^2 + a[3]^2) >= 0)
    @NLconstraint(model, κ[1] * κ[2] * κ[3] - (κ[1] * a[1]^2 + κ[2] * a[2]^2 + κ[3] * a[3]^2) >= 0)

    # Objective function
    @NLobjective(model, Min, kossak_objectiveGEXY(ϵ, h_Re, h_Im, γ[1], γ[2], γ[3], a[1], a[2], a[3]))

    # Solve the model
    JuMP.optimize!(model)

    # Retrieve and print results
    objective_value = JuMP.objective_value(model)
    ϵ_value = value(ϵ)
    h_Re_value = value.(h_Re)
    h_Im_value = value.(h_Im)
    γ_value = value.(γ)
    a_value = value.(a)
    C_Re_value = value.(C_Re)
    C_Im_value = value.(C_Im)

    #println("Initial Guess: ", initial_guess)
    println("Objective Value: ", objective_value)
    println("ϵ: ", ϵ_value)
    #println("γ: ", γ_value)
    #println("a: ", a_value)
    #println("C_Re: ", C_Re_value)
    #println("C_Im: ", C_Im_value)

    return objective_value, ϵ_value, h_Re_value, h_Im_value, γ_value, a_value
end

function test_and_save_D10(Hˢⁱᵈ, Jˢⁱᵈ, γᵢ, output_file_name)

    γ = ["0.079477", "0.25133", "0.79477", "2.5133", "7.9477", "25.133", "79.477", "251.33"]
    
    tests_data_file_name = res_dir * "BENCHMARK_TEST_" * output_file_name * ".h5"

    ρᵍ₀ = [ 1 0.
            0 0 ]    # state to measure initial distance from

    dodeca_10_states = ["D"*string(n) for n=1:10];

    basis_states = ["B"*string(n) for n=1:4];

    train_states = basis_states 
    test_states = dodeca_10_states;

    FminStates = []
    FmedianStates = []
    FmeanStates = []

    h5open(tests_data_file_name, "cw") do fid
        γ_group = create_group(fid, γᵢ) # create coupling group "gamma_" *
        #γ_group["H"] = convert.(ComplexF64, Hˢⁱᵈ)
        #γ_group["C"] = convert.(ComplexF64, Cˢⁱᵈ)
    end

    for state in test_states # loop over initial states
        
        print(state*" ")

        start_time = time()

        tₛ, ρₛ = LiPoSID.read_timeevolution(evol_data_file_name, state, γᵢ)
        ρₛ = convert(Vector{Matrix{ComplexF64}}, ρₛ)
        
        ρᵗˢᵗ = [DenseOperator(basis,Hermitian(ρₜ)) for ρₜ in ρₛ]
        tᵗˢᵗ = convert.(Float64, tₛ)

        ρₒ = DenseOperator(basis,ρₛ[1])
        dt = tᵗˢᵗ[2] - tᵗˢᵗ[1]
        tᵉⁿᵈ = tᵗˢᵗ[end]

        #print("effective_Lindblad_ops for Kossakowski")       
        effective_Lindblad_ops = [DenseOperator(basis,j) for j in Jˢⁱᵈ]

        #print("Simulating Kossakowski")
        tout, ρ_t_kossak = timeevolution.master(tᵗˢᵗ, ρₒ, DenseOperator(basis, Hˢⁱᵈ), effective_Lindblad_ops)
        ρˢⁱᵈ  = [ρₜ.data for ρₜ in ρ_t_kossak]

        #print("Calculating Fidelity")

        #F = LiPoSID.fidelity_series(basis, [ρₜ.data for ρₜ in ρˢⁱᵐ], ρˢⁱᵈ)
        F = LiPoSID.fidelity_series(basis, ρₛ, ρˢⁱᵈ)

        h5open(tests_data_file_name, "cw") do fid
            #γ_group = create_group(fid, "gamma_" * γᵢ) # open coupling group
            γ_group = open_group(fid, γᵢ) # open coupling group "gamma_" *
            init_state_group = create_group(γ_group, state) # create initial state group
            init_state_group["Fidelity"] = convert.(Float64, F)
        end
        
        FminState = minimum(F)
        FmedianState = median(F)
        FmeanState = mean(F)
        
        push!(FminStates, FminState)
        push!(FmedianStates, FmedianState)
        push!(FmeanStates, FmeanState)

    end

    # Calculate the mean
    F_mean_value = mean(FmeanStates)

    # Calculate the median
    F_median_value = median(FmedianStates)

    # Calculate the min
    F_min_value = minimum(FminStates)

    println()
    println("Mimimal fidelity for "*γᵢ*": ", F_min_value)
    println("Median fidelity for "*γᵢ*": ", F_median_value)

    return(F_min_value, F_median_value)
    
end

# Define the density matrix evolution and time points

ρᵍ₀ = [ 1 0.
        0 0 ]    # state to measure initial distance from

dodeca_10_states = ["D"*string(n) for n=1:10];

basis_states = ["B"*string(n) for n=1:4];

data_dir = "../data/"

train_states = basis_states 
test_states = dodeca_10_states

date_and_time_string = string(Dates.format(now(), "yyyy-u-dd_at_HH-MM"))

evol_data_file_name = data_dir * "ALL_GAMMAS_B4_D10.h5"

res_dir = "../results/"

model_name = "KossakConstrFrob_Jump_NLOPT_LD_SLSQP"

models_data_file_name = model_name #* date_and_time_string 

γ = [ "0.079477",  "0.25133", "0.79477", "2.5133", "7.9477", "25.133", "79.477", "251.33"]

global tᵍᵉˣʸ
global ρᵍᵉˣʸ

for γᵢ in  γ

    println("γ=",γᵢ)

    tρ_GEXY = LiPoSID.read_GEXY_timeevolution(evol_data_file_name, γᵢ)

    global tᵍᵉˣʸ = tρ_GEXY[1]
    global ρᵍᵉˣʸ = tρ_GEXY[2] 

    γᶠ = parse(ComplexF64, γᵢ)

    # Define a general upper bound for all variables
    upper_bound = 260

    #Smart guesses for pure relaxation
    smart_guess = zeros(9)
    smart_guess[1] = 25.133 # ϵ 
    smart_guess[2] = 0.0
    smart_guess[3] = 0.0
    smart_guess[4] = γᶠ/2  # γ1
    smart_guess[5] = γᶠ/2  # γ2
    smart_guess[6] = γᶠ    # γ3

    elapsed_time = @timed begin

        object_val, ϵ_value, h_Re_value, h_Im_value, γ_value, a_value = solve_kossak_model(smart_guess, upper_bound)
    
    end

    # Substitute the optimized parameters back into H and C
    optimized_H = construct_H(ϵ_value, h_Re_value, h_Im_value)
    optimized_C = construct_C(γ_value, a_value)

    Hˢⁱᵈ = convert.(ComplexF64,optimized_H)
    Cˢⁱᵈ = convert.(ComplexF64,optimized_C)

    effective_Lindblad = LiPoSID.get_lindblad_operators(Cˢⁱᵈ, fᴼᴺᴮ)

    h5open(res_dir * "MODEL_" * models_data_file_name * ".h5", "cw") do fid
        γ_group = create_group(fid, γᵢ) # create coupling group "gamma_" *
        γ_group["H"] = convert.(ComplexF64, Hˢⁱᵈ)
        γ_group["C"] = convert.(ComplexF64, Cˢⁱᵈ)
        γ_group["runtime"] = elapsed_time.time
    end

    #test_D10(Hˢⁱᵈ, effective_Lindblad, γᵢ)

    test_and_save_D10(Hˢⁱᵈ, effective_Lindblad, γᵢ, models_data_file_name)
    
end