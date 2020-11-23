module MultiLJ

using JuLIP
import JuLIP: JMat, rnn
import JuLIP.Potentials: cutoff, evaluate!, evaluate_d!, MPairPotential

export MLJ

struct MLJ{T} <: MPairPotential
   Z2idx::Vector{Int} 
   rcut::Float64
   V::Matrix{T}
end

function MLJ(Z, ϵ0, σ0; rcutfact = 2.7)
   @assert length(Z) == length(unique(Z))
   # create a mapping from atomic numbers to indices
   Z2idx = zeros(Int, maximum(Z))
   Z2idx[Z] = 1:length(Z)
   # generate the potential
   rcut = maximum(σ0) * rcutfact
   V = [ lennardjones(; rcut = :auto, σ = σ0[a,b], e0 = ϵ0[a, b])
         for a = 1:length(Z), b = 1:length(Z) ]
   return MLJ(Z2idx, rcut, V)
end

cutoff(V::MLJ) = V.rcut

function evaluate!(tmp, V::MLJ{T}, r12::AbstractFloat, z1::Integer, z0::Integer) where {T}
    return V.V[z1, z0](r12)
end

function evaluate_d!(tmp, V::MLJ{T}, r12, z1, z0) where {T}
    return (@D V.V[z1, z0](r12))
end


#function precon!(tmp, V::HertzModel{T}, r12, R12, z1, z0, innerstab) where {T}
#    r1, r2 = V.radii[z1], V.radii[z0]
#    return r1 + r2 < r12 ? zero(JMat{T}) : precon!(tmp, V.fpair, r12, R12, innerstab)
#end

end