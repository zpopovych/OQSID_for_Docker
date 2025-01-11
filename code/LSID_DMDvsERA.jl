include("LiPoSID.jl")
using QuantumOptics
basis = NLevelBasis(2)
using LinearAlgebra

using Dates
using HDF5

function DMD_SVD(Y, r, Δt)
    
    X₋ = Y[:,1:end-1]
    X₊ = Y[:,2:end]
    U, Σ, V = svd(X₋)
    
    Uʳ = U[:, 1:r] #12 x 4
    Σʳ = diagm(Σ[1:r])
    Vʳ = V[:, 1:r]
    Ã = Uʳ' * X₊ * Vʳ / Σʳ
    Λ, W = eigen(Ã)
    Φ = X₊ * Vʳ / Σʳ * W
    Ω = log.(Λ)/Δt
    x₁ = X₋[:,1]
    b₁ = Φ \ x₁
    
    return Φ, Ω, b₁, Ã

end    

function DMD_reconstruct(Φ, Ω, b₁, Δt, steps)

    Yᵈᵐᵈ = hcat([real.(Φ * (b₁ .* exp.(Ω * (i * Δt)))) for i in 0:steps]...)
    
    return Yᵈᵐᵈ

end


function DMDvsERA_sb_basis(γᵢ, n)
    
    γᵢ = string(γᵢ)

    tᵍᵉˣʸ , ρᵍᵉˣʸ  = LiPoSID.read_GEXY_timeevolution(evol_data_file_name, γᵢ)

    ρᵍ, ρᵉ, ρˣ, ρʸ = ρᵍᵉˣʸ
    tᵍ, tᵉ, tˣ, tʸ = tᵍᵉˣʸ
    
    lᵉ = length(ρᵉ); lᵍ = length(ρᵍ); lˣ = length(ρˣ); lʸ = length(ρʸ)
    lᵐᵃˣ = min(lᵉ, lᵍ,  lˣ, lʸ)  #choose time limit by shortest series

    tᵉᶠ = convert.(Float64, tᵉ); tᵍᶠ = convert.(Float64, tᵍ); 
    tˣᶠ = convert.(Float64, tˣ); tʸᶠ = convert.(Float64, tʸ); 

    bᵉ = LiPoSID.bloch(ρᵉ[1:lᵐᵃˣ])
    bᵍ = LiPoSID.bloch(ρᵍ[1:lᵐᵃˣ])
    bˣ = LiPoSID.bloch(ρˣ[1:lᵐᵃˣ])
    bʸ = LiPoSID.bloch(ρʸ[1:lᵐᵃˣ])

    Y = [bᵉ; bᵍ; bˣ; bʸ]

    t = convert.(Float64, tᵉᶠ[1:lᵐᵃˣ])
    Δt = t[2]-t[1]

    # DMD (Dynamic mode decomposition)
    
    Φ, Ω, b₁, Aᴰᴹᴰ = DMD_SVD(Y, n, Δt)

    Aᴰᴹᴰc = log(Aᴰᴹᴰ)/Δt
    Λᴰᴹᴰ, Wᴰᴹᴰ = eigen(Aᴰᴹᴰc)

    Yᴰᴹᴰ = DMD_reconstruct(Φ, Ω, b₁, Δt, length(t))

    bᵉᴰᴹᴰ = Yᴰᴹᴰ[1:3,:]
    bᵍᴰᴹᴰ = Yᴰᴹᴰ[4:6,:]
    bˣᴰᴹᴰ = Yᴰᴹᴰ[7:9,:]
    bʸᴰᴹᴰ = Yᴰᴹᴰ[10:12,:]

    ρᵉᴰᴹᴰ = LiPoSID.rho_series_from_bloch(bᵉᴰᴹᴰ)
    ρᵍᴰᴹᴰ = LiPoSID.rho_series_from_bloch(bᵍᴰᴹᴰ)
    ρˣᴰᴹᴰ = LiPoSID.rho_series_from_bloch(bˣᴰᴹᴰ)
    ρʸᴰᴹᴰ = LiPoSID.rho_series_from_bloch(bʸᴰᴹᴰ)

    ρᴰᴹᴰ = [ρᵉᴰᴹᴰ, ρᵍᴰᴹᴰ, ρˣᴰᴹᴰ, ρʸᴰᴹᴰ]

    # ERA (Eigenvalue Realization Algorithm)

    Aᴱᴿᴬ, Cᴱᴿᴬ, x₀ᴱᴿᴬ, Σᴱᴿᴬ = LiPoSID.lsid_n_ACx0Σ(Y, Δt, n) 

    Aᴱᴿᴬc = log(Aᴱᴿᴬ)/Δt

    Λᴱᴿᴬ, Wᴱᴿᴬ = eigen(Aᴱᴿᴬc)

    Yᴱᴿᴬ = LiPoSID.propagate_LTI(Aᴱᴿᴬ, Cᴱᴿᴬ, x₀ᴱᴿᴬ, n, length(t))

    bᵉᴱᴿᴬ = Yᴱᴿᴬ[1:3,:]
    bᵍᴱᴿᴬ = Yᴱᴿᴬ[4:6,:]
    bˣᴱᴿᴬ = Yᴱᴿᴬ[7:9,:]
    bʸᴱᴿᴬ = Yᴱᴿᴬ[10:12,:]

    ρᵉᴱᴿᴬ = LiPoSID.rho_series_from_bloch(bᵉᴱᴿᴬ)
    ρᵍᴱᴿᴬ = LiPoSID.rho_series_from_bloch(bᵍᴱᴿᴬ)
    ρˣᴱᴿᴬ = LiPoSID.rho_series_from_bloch(bˣᴱᴿᴬ)
    ρʸᴱᴿᴬ = LiPoSID.rho_series_from_bloch(bʸᴱᴿᴬ)

    ρᴱᴿᴬ = [ρᵉᴱᴿᴬ, ρᵍᴱᴿᴬ, ρˣᴱᴿᴬ, ρʸᴱᴿᴬ]
    
    return ρᴱᴿᴬ, ρᴰᴹᴰ, Λᴱᴿᴬ, Λᴰᴹᴰ, t[1:lᵐᵃˣ]
    
end

function propagate_rho_O1XY(ρ₀, ρᵉᵍˣʸ, steps)

    hcat(vec[ρᵢ] for ρᵢ in ρᵉᵍˣʸ)

    kᵉᵍˣʸ = hcat([vec(ρᵢ[1]) for ρᵢ in ρᵉᵍˣʸ]...)\vec(ρ₀)

    kᵉ, kᵍ, kˣ, kʸ = kᵉᵍˣʸ              
    ρᵉ, ρᵍ, ρˣ, ρʸ = ρᵉᵍˣʸ

    ρ = kᵉ * ρᵉ + kᵍ * ρᵍ + kˣ * ρˣ + kʸ * ρʸ

    return ρ
end 

parentdir = "../"
data_dir = parentdir*"data/"
println(data_dir)

models_dir = parentdir*"results/"
tests_dir = parentdir*"results/"

dodeca_files = ["D"*string(n) for n=1:10];
basis_files = ["B"*string(n) for n=1:4];

all_files = vcat(dodeca_files, basis_files)
test_states = dodeca_files;

evol_data_file_name = parentdir * "data/ALL_GAMMAS_B4_D10.h5"

date_and_time_string =  string(Dates.format(now(), "yyyy-u-dd_at_HH-MM"))

tests_data_file_name = "DMDvsERA_rank5_SB_trn4_tst10.h5" #_"*date_and_time_string * ".h5"

γ = [ "0.079477",  "0.25133", "0.79477", "2.5133", "7.9477", "25.133",  "79.477", "251.33"]

println("Coupling levels to be avaluated γ ∈ ", γ)

n = 5

for γᵢ in  γ

    ρᴱᴿᴬ, ρᴰᴹᴰ, Λᴱᴿᴬ, Λᴰᴹᴰ, t = DMDvsERA_sb_basis(γᵢ, n)

    ρᵉᴰᴹᴰ, ρᵍᴰᴹᴰ, ρˣᴰᴹᴰ, ρʸᴰᴹᴰ = ρᴰᴹᴰ
    ρᵉᴱᴿᴬ, ρᵍᴱᴿᴬ, ρˣᴱᴿᴬ, ρʸᴱᴿᴬ = ρᴱᴿᴬ

    h5open(tests_dir*tests_data_file_name,"cw") do fid
        γ_group = create_group(fid, γᵢ) # create coupling group
        γ_group["Re_eigvals_dmd_sb"] = convert.(Float64, real.(Λᴰᴹᴰ))
        γ_group["Re_eigvals_era_sb"] = convert.(Float64, real.(Λᴱᴿᴬ))
    end 


    for state in test_states # loop over initial states

        tᵗˢᵗ, ρᵗˢᵗ = LiPoSID.read_timeevolution(evol_data_file_name, state, γᵢ)
        
        if length(tᵗˢᵗ) > 1200 end_tst = 1200 else end_tst = length(tᵗˢᵗ) end
            
        ρᵗˢᵗ = convert(Vector{Matrix{ComplexF64}}, ρᵗˢᵗ[1:end_tst])

        steps = min(end_tst, length(t))

        ρᵗˢᵗᴱᴿᴬ =  propagate_rho_O1XY(ρᵗˢᵗ[1], ρᴱᴿᴬ, steps)
        ρᵗˢᵗᴰᴹᴰ =  propagate_rho_O1XY(ρᵗˢᵗ[1], ρᴰᴹᴰ, steps)

        Fᴱᴿᴬ = LiPoSID.fidelity_series(basis, ρᵗˢᵗᴱᴿᴬ[1:steps], ρᵗˢᵗ[1:steps])
        Fᴰᴹᴰ = LiPoSID.fidelity_series(basis, ρᵗˢᵗᴰᴹᴰ[1:steps], ρᵗˢᵗ[1:steps])
        
        h5open(tests_dir*tests_data_file_name,"cw") do fid
            γ_group = open_group(fid, γᵢ) # create coupling group
            state_group = create_group(γ_group, state) # create coupling group

            state_group["F_dmd_sb"] = convert.(Float64, Fᴰᴹᴰ)
            state_group["F_era_sb"] = convert.(Float64, Fᴱᴿᴬ)
            state_group["time"] = convert.(Float64, tᵗˢᵗ[1:steps]) 
            
        end

    end
    
end

println(tests_data_file_name * " done.")