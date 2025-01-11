using LinearAlgebra
using Combinatorics

using HDF5

include("LiPoSID.jl")

function NonMarkovianity(file_name, states, γᵢ)
    N=[]
    for (i,j)  in combinations(1:length(states), 2)
        t₁, ρs₁ = LiPoSID.read_timeevolution(file_name, states[i],γᵢ)
        t₂, ρs₂ = LiPoSID.read_timeevolution(file_name, states[j],γᵢ)
        dD = diff([LiPoSID.TrDist(ρ₁, ρ₂) for (ρ₁, ρ₂) in zip(ρs₁, ρs₂)])
        append!(N, sum(dD[dD.>0]))
    end
    maximum(N)
end  

file_name = "../data/ALL_GAMMAS_B4_D10.h5"

dodeca_10_states = ["D"*string(n) for n=1:10];
basis_states = ["B"*string(n) for n=1:4];

all_states = vcat(basis_states, dodeca_10_states);

γ = [ "0.079477",  "0.25133", "0.79477", "2.5133", "7.9477", "25.133", "79.477", "251.33"]

N = [NonMarkovianity(file_name, all_states, γᵢ) for γᵢ in γ]

h5open("../results/NonMark.h5", "cw") do fid
    fid["N"] = N
end 