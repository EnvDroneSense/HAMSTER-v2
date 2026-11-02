"Change between the two time steps (columns) of `x`; changes of at most 1e-4 in size are set to 0."
function delta(x)
    d = x[:, 2] .- x[:, 1]
    return @. ifelse(abs(d) <= 1.0f-4, zero(d), d)
end
