"Initial empty mass distribution."
struct EmptyMassDistribution{N} <: InitialMassDistribution{N} end

Base.show(io::IO, ::Type{EmptyMassDistribution{N}}) where {N} =
    print(io, "EmptyMassDistribution")

(::EmptyMassDistribution{N})() where {N} = MassDistributionData{N}(1.0, [], [])
