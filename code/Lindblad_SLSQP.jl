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
function DcLindblad(ρ, t, H, J)
    U = (H * ρ - ρ * H) / im
    D = sum([ ( Jⱼ * ρ * Jⱼ' - (Jⱼ' * Jⱼ  * ρ + ρ * Jⱼ' * Jⱼ)/2 ) for Jⱼ in J ])
    return U + D
end

# Define functions to construct H and C from parameters
function construct_H(ϵ, h_Re, h_Im)
    return [
             ϵ                     h_Re + im * h_Im
             h_Re - im * h_Im     -ϵ
           ] / 2
end

function construct_J1(r)
    return [
         0   r
         0   0    
    ]
end

function construct_J2(p)
    return [
         0   p
         p   0    
    ]
end

function construct_J3(p)
    return [
         0     -im*p
         im*p   0    
    ]
end

function construct_J4(p)
    return [
         p    0
         0   -p   
    ]
end

# Objective function

function lindblad_obj(ϵ, h_Re, h_Im, r, p, ρ, t) #Simpson 
    H = construct_H(ϵ, h_Re, h_Im)
    J1 = construct_J1(r)
    J2 = construct_J2(p)
    J3 = construct_J3(p)
    J4 = construct_J4(p)

    J = [J1, J2, J3, J4]
    obj = 0.0
    for i in 3:length(ρ)
        ρ1 = ρ[i]
        ρ2 = ρ[i - 2] + 
        (t[i] - t[i - 1]) * 
        (DcLindblad(ρ[i], t[i], H, J1) +
         4 * DcLindblad(ρ[i - 1], t[i - 1], H, J) 
         + DcLindblad(ρ[i - 2], t[i - 2], H, J)) / 3
        #obj += infidelity_norm(ρ1,ρ2)
        #obj += fidelity(ρ1,ρ2)
        obj += frobenius_norm2(ρ1-ρ2)
    end
    return obj
end

# Define the objective function wrapper
function lindblad_objectiveGEXY(ϵ, h_Re, h_Im, r, p)

    ρᵍ, ρᵉ, ρˣ, ρʸ = ρᵍᵉˣʸ
    tᵍ, tᵉ, tˣ, tʸ = tᵍᵉˣʸ

    objGEXY = lindblad_obj(ϵ, h_Re, h_Im, r, p, ρᵍ, tᵍ) + 
              lindblad_obj(ϵ, h_Re, h_Im, r, p, ρᵉ, tᵉ) + 
              lindblad_obj(ϵ, h_Re, h_Im, r, p, ρˣ, tˣ) + 
              lindblad_obj(ϵ, h_Re, h_Im, r, p, ρʸ, tʸ)
    return objGEXY
end


# Define a function to create and solve the model with given initial guesses and a general upper bound
function solve_lindblad_model(initial_guess, upper_bound) #, γᵢ)

    #upper_bound = 260
    
    model = Model(optimizer_with_attributes(NLopt.Optimizer, "algorithm" => :LD_SLSQP))

    # Register the element-wise power operator
    register(model, :.^, 2, (x, y) -> x.^y, autodiff=true)

    # Register the custom objective function
    register(model, :lindblad_objectiveGEXY, 5, lindblad_objectiveGEXY, autodiff = true)

    # Define lower and upper bounds
    lower_bound = -upper_bound

    # Define variables with general upper and lower bounds and initial guesses
    @variable(model, 0 <= ϵ <= upper_bound, start = initial_guess[1])
    @variable(model, lower_bound <= h_Re <= upper_bound, start = initial_guess[2])
    @variable(model, lower_bound <= h_Im <= upper_bound, start = initial_guess[3])
    @variable(model, lower_bound <= r <= upper_bound, start = initial_guess[4])
    @variable(model, lower_bound <= p <= upper_bound, start = initial_guess[5])

    # Objective function
    @NLobjective(model, Min, lindblad_objectiveGEXY(ϵ, h_Re, h_Im, r, p))

    # Solve the model
    JuMP.optimize!(model)

    # Retrieve and print results
    objective_value = JuMP.objective_value(model)
    ϵ_value = value(ϵ)
    h_Re_value = value(h_Re)
    h_Im_value = value(h_Im)
    r_value = value(r)
    p_value = value(p)


    #println("Initial Guess: ", initial_guess)
    println("Objective Value: ", objective_value)
    println("ϵ: ", ϵ_value)
    println("r: ", r_value)
    println("p: ", p_value)

    return objective_value, ϵ_value, h_Re_value, h_Im_value, r_value, p_value
end

# Define the density matrix evolution and time points

ρᵍ₀ = [ 1 0.
        0 0 ]    # state to measure initial distance from

dodeca_10_states = ["D"*string(n) for n=1:10];

basis_states = ["B"*string(n) for n=1:4];

train_states = basis_states 
test_states = dodeca_10_states;

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

using Dates, Random

date_and_time_string = string(Dates.format(now(), "yyyy-u-dd_at_HH-MM"))

data_dir = "../data/"

evol_data_file_name = data_dir * "ALL_GAMMAS_B4_D10.h5"

res_dir = "../results/"

model_name = "LindbladFrob_Jump_NLOPT_LD_SLSQP"

models_data_file_name = model_name

γ = [ "0.079477",  "0.25133", "0.79477", "2.5133", "7.9477", "25.133", "79.477", "251.33"]

global tᵍᵉˣʸ
global ρᵍᵉˣʸ

for γᵢ in  γ
    println("γ=",γᵢ)

    γᶠ = parse(ComplexF64, γᵢ)

    tρ_GEXY = LiPoSID.read_GEXY_timeevolution(evol_data_file_name, γᵢ)

    global tᵍᵉˣʸ = tρ_GEXY[1]
    global ρᵍᵉˣʸ = tρ_GEXY[2] 

    # Define a general upper bound for all variables
    upper_bound = 260

    #Smart guesses for pure relaxation
    smart_guess = zeros(5)
    smart_guess[1] = 25.133 # ϵ 
    smart_guess[2] = 0.0
    smart_guess[3] = 0.0
    smart_guess[4] = sqrt(γᶠ)  # r
    smart_guess[5] = 0.0  # p


    elapsed_time = @timed begin

        object_val, ϵ_value, h_Re_value, h_Im_value, r_value, p_value = solve_lindblad_model(smart_guess, upper_bound)
    
    end

    # Substitute the optimized parameters back into H and C
    optimized_H = construct_H(ϵ_value, h_Re_value, h_Im_value)
    optimized_J1 = construct_J1(r_value)
    optimized_J2 = construct_J2(p_value)
    optimized_J3 = construct_J3(p_value)
    optimized_J4 = construct_J4(p_value)

    Hˢⁱᵈ = convert.(ComplexF64,optimized_H)
    J1ˢⁱᵈ = convert.(ComplexF64,optimized_J1)
    J2ˢⁱᵈ = convert.(ComplexF64,optimized_J2)
    J3ˢⁱᵈ = convert.(ComplexF64,optimized_J3)
    J4ˢⁱᵈ = convert.(ComplexF64,optimized_J4)

    effective_Lindblad = [J1ˢⁱᵈ, J2ˢⁱᵈ, J3ˢⁱᵈ, J4ˢⁱᵈ]

    h5open(res_dir * "MODEL_" * models_data_file_name * ".h5", "cw") do fid
        γ_group = create_group(fid, γᵢ) # create coupling group "gamma_" *
        γ_group["H"] = convert.(ComplexF64, Hˢⁱᵈ)
        γ_group["J1"] = convert.(ComplexF64, J1ˢⁱᵈ)
        γ_group["J2"] = convert.(ComplexF64, J2ˢⁱᵈ)
        γ_group["J3"] = convert.(ComplexF64, J3ˢⁱᵈ)
        γ_group["J4"] = convert.(ComplexF64, J4ˢⁱᵈ)
        γ_group["runtime"] = elapsed_time.time
    end

    #test_D10(Hˢⁱᵈ, effective_Lindblad, γᵢ)

    test_and_save_D10(Hˢⁱᵈ, effective_Lindblad, γᵢ, models_data_file_name)
    
end