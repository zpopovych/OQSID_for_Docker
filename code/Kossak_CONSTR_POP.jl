include("LiPoSID.jl")

using QuantumOptics
basis = NLevelBasis(2)
using LinearAlgebra

using HDF5
using DynamicPolynomials

using Dates

using Statistics

using TSSOS

σˣ = [ 0 1 
       1 0 ]

σʸ = [ 0.   im*1
      -im*1 0    ]

σᶻ = [ 1.  0
       0  -1 ] 

fᴷ₁ = σˣ/2; fᴷ₂ = σʸ/2; fᴷ₃ = σᶻ/2

@assert tr(σˣ/2*σʸ/2) == tr(σˣ/2*σᶻ/2) ==  tr(σʸ/2*σᶻ/2) ≈ 0
@assert tr(σˣ/2*σˣ/2) == tr(σʸ/2*σʸ/2) == tr(σᶻ/2*σᶻ/2) ≈ 1/2

fᴷᴼᴺᴮ = [fᴷ₁, fᴷ₂, fᴷ₃]

@polyvar ϵ h₁ h₂ 

Hˢʸᵐᵇ = [ ϵ           h₁+im*h₂
           h₁-im*h₂   -ϵ            ] / 2

Hᴷˢʸᵐᵇ = h₁ * fᴷ₁ + h₂ * fᴷ₂  + ϵ * fᴷ₃ 

@assert tr(Hᴷˢʸᵐᵇ) == 0
@assert Hᴷˢʸᵐᵇ == Hˢʸᵐᵇ

@polyvar κ[1:3]; κ₁,κ₂,κ₃ = κ[1],κ[2],κ[3]
@polyvar a[1:3]; a₁,a₂,a₃ = a[1],a[2],a[3]

Cˢʸᵐᵇ = [  κ[1]     -im*a[3]     im*a[2]
           im*a[3]   κ[2]       -im*a[1] 
          -im*a[2]   im*a[1]     κ[3]     ]

constr1 = κ₁ + κ₂ + κ₃ 
constr2 = κ₁*κ₂ + κ₃*κ₁ + κ₂*κ₃ - - a₁^2 - a₂^2 - a₃^2
constr3 = κ₁*κ₂*κ₃ - κ₁*a₁^2 - κ₂*a₂^2 - κ₃*a₃^2

constraints = [ ϵ, κ₁, κ₂, κ₃, constr1, constr2, constr3] 

function kossak_obj(ρ, t, Hˢʸᵐᵇ, Cˢʸᵐᵇ, Fᴼᴺᴮ)

    function Dc(ρ, t)
        U = (Hˢʸᵐᵇ*ρ - ρ*Hˢʸᵐᵇ)/im 
        D = sum(Cˢʸᵐᵇ .* [2*fᵢ*ρ*fⱼ' - ρ*fⱼ'*fᵢ - fⱼ'*fᵢ*ρ  for fᵢ in Fᴼᴺᴮ, fⱼ in Fᴼᴺᴮ])/2
        return U + D
    end 

    obj = 0
    for i in 3:length(ρ)
        obj += LiPoSID.frobenius_norm2(
            ρ[i] - ρ[i-2] - (t[i]-t[i-1])*(Dc(ρ[i], t[i])+
            4*Dc(ρ[i-1], t[i-1])+Dc(ρ[i-2], t[i-2]))/3
        )
    end

    if isempty(monomials(obj))
        obj = 0. 
    else
        obj = sum(real(coef) * mon for (coef, mon) in zip(coefficients(obj), monomials(obj)))
    end

    return obj

end

function kossak_GEXY_obj(ρᵍᵉˣʸ, tᵍᵉˣʸ, H0ˢʸᵐᵇ, Cˢʸᵐᵇ, fᴼᴺᴮ)

    ρᵍ, ρᵉ, ρˣ, ρʸ = ρᵍᵉˣʸ
    tᵍ, tᵉ, tˣ, tʸ = tᵍᵉˣʸ

    polyG = kossak_obj(ρᵍ, tᵍ, H0ˢʸᵐᵇ, Cˢʸᵐᵇ, fᴼᴺᴮ)
    polyE = kossak_obj(ρᵉ, tᵉ, H0ˢʸᵐᵇ, Cˢʸᵐᵇ, fᴼᴺᴮ)
    polyX = kossak_obj(ρˣ, tˣ, H0ˢʸᵐᵇ, Cˢʸᵐᵇ, fᴼᴺᴮ)
    polyY = kossak_obj(ρʸ, tʸ, H0ˢʸᵐᵇ, Cˢʸᵐᵇ, fᴼᴺᴮ)

    polyGEXY = polyG + polyE + polyX + polyY

    return polyGEXY
end

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

println(" SYSTEM IDENTIFICATION w CONSTRAINED TSSOS and KOSSAKOWSKI Frobenius objective QO simulation")

γ = [ "0.079477",  "0.25133", "0.79477", "2.5133", "7.9477", "25.133", "79.477", "251.33"]

date_and_time_string =  string(Dates.format(now(), "yyyy-u-dd_at_HH-MM"))

data_dir = "../data/"

evol_data_file_name = data_dir * "ALL_GAMMAS_B4_D10.h5"

relative_threshold = 1e-15

rltrs = string(convert(Int, floor(log10(relative_threshold))))

tests_dir = "../results/"

tests_data_file_name = "KOSSAK_CONSTR_TSSOS_treshold_1e"*rltrs*"_FROB_QO.h5" #*date_and_time_string * ".h5"

FminGammas = []
FmedianGammas = []
FmeanGammas = []
Epsilons = []
CoefRanges = []

for γᵢ in γ

    println("γ =  "*γᵢ)

    tᵍᵉˣʸ , ρᵍᵉˣʸ  = LiPoSID.read_GEXY_timeevolution(evol_data_file_name, γᵢ)

    elapsed_time = @timed begin

        polyGEXYfull = kossak_GEXY_obj(ρᵍᵉˣʸ, tᵍᵉˣʸ, Hˢʸᵐᵇ, Cˢʸᵐᵇ, fᴷᴼᴺᴮ)     
        #polyGEXY = polyGEXYfull
        polyGEXY = LiPoSID.filter_odd_terms_by_relative_threshold(polyGEXYfull, relative_threshold)
        sol, status = min_cs_tssos(polyGEXY, constraints)

    end

    @show minimum(abs.(coefficients(polyGEXY)))
    @show maximum(abs.(coefficients(polyGEXY)))
    push!(CoefRanges, LiPoSID.coefficient_range(polyGEXYfull))

    print(" status:", status)
    print(" runtime :", elapsed_time.time)

    Hˢⁱᵈ = convert.(ComplexF64, subs(Hˢʸᵐᵇ, sol))
    Cˢⁱᵈ = convert.(ComplexF64,subs(Cˢʸᵐᵇ, sol))
    epsilon = subs(ϵ, sol)

    push!(Epsilons, epsilon)
    
    h5open(tests_dir*tests_data_file_name,"cw") do fid
        γ_group = create_group(fid, γᵢ) # create gamma coupling group   
        γ_group["epsilon"] = convert(Float64, epsilon)
        γ_group["H"] = convert.(ComplexF64, Hˢⁱᵈ)
        γ_group["C"] = convert.(ComplexF64, Cˢⁱᵈ)
        γ_group["status"] = status
        γ_group["runtime"] = elapsed_time.time

    end

    println()

    FminStates = []
    FmedianStates = []
    FmeanStates = []

    for state in test_states # loop over initial states
        
        print(state*" ")

        start_time = time()

        tₛ, ρₛ = LiPoSID.read_timeevolution(evol_data_file_name, state, γᵢ)
        ρₛ = convert(Vector{Matrix{ComplexF64}}, ρₛ)
        #bᵗˢᵗ = LiPoSID.bloch(ρₛ)
        ρᵗˢᵗ = [DenseOperator(basis,Hermitian(ρₜ)) for ρₜ in ρₛ]
        tᵗˢᵗ = convert(Vector{Float64}, tₛ)

        ρₒ = DenseOperator(basis, ρₛ[1])
        dt = tᵗˢᵗ[2] - tᵗˢᵗ[1]
        tᵉⁿᵈ = tᵗˢᵗ[end]

        #print("effective_Lindblad_ops for Kossakowski")

        effective_Lindblad = LiPoSID.get_lindblad_operators(convert.(ComplexF64, Cˢⁱᵈ), fᴷᴼᴺᴮ)
        effective_Lindblad_ops = [DenseOperator(basis,j) for j in effective_Lindblad]

        #print("Simulating Kossakowski")

        tout, ρ_t_kossak = timeevolution.master(convert.(Float64, tᵗˢᵗ), ρₒ, DenseOperator(basis, Hˢⁱᵈ), effective_Lindblad_ops)
        ρˢⁱᵈ  = [ρₜ.data for ρₜ in ρ_t_kossak]

        #bˢⁱᵈ = LiPoSID.bloch(ρˢⁱᵈ)
        
        # Calculating fidelity
        
        F = LiPoSID.fidelity_series(basis, ρₛ, ρˢⁱᵈ)


        h5open(tests_dir*tests_data_file_name,"cw") do fid
            γ_group = open_group(fid, γᵢ) # open gamma coupling group
            init_state_group = create_group(γ_group, state) # create initial state group
            init_state_group["Fidelity"] = convert.(Float64, F)
            #init_state_group["bloch_exact"] = convert.(Float64, bᵗˢᵗ)
            #init_state_group["bloch_sid"] = convert.(Float64, bˢⁱᵈ)
            init_state_group["tr_dist_grnd"] = LiPoSID.TrDist(ρₛ[1], ρᵍ₀)
            init_state_group["time"] = tᵗˢᵗ
        end
        
        FminState = minimum(F)
        FmedianState = mean(F)
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

    push!(FminGammas, F_min_value)
    push!(FmedianGammas, F_median_value)
    push!(FmeanGammas, F_mean_value)

    println("Median fidelity for "*γᵢ*": ", F_median_value)

end

h5open(tests_dir*tests_data_file_name,"cw") do fid
    fid["F_min"] = convert.(Float64,FminGammas)
    fid["F_median"] = convert.(Float64,FmedianGammas)
    fid["F_mean"] = convert.(Float64,FmeanGammas)
    fid["Energy"] = convert.(Float64,Epsilons)
    fid["CoefsRanges"] = convert.(Float64,CoefRanges)

end


println(tests_data_file_name)