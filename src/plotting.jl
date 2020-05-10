export sphplot

@recipe f(::Type{Vector{SVectorF{1}}}, array::Vector{SVectorF{1}}) =
    [x[] for x in array]

@userplot SPHPlot
@recipe function f(p::SPHPlot)
    x, y = p.args
    seriestype := :scatter
    legend := false
    x, y
end
