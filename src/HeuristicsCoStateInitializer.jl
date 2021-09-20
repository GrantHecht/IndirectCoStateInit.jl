
abstract type HeuristicsCoStateInitializer{HOT <: Heuristics.Optimizer} end

function initialize!(hcsi::HeuristicsCoStateInitializer)
    optimize!(hcsi.ho, hcsi.hoOpts)
    hcsi.λh .= hcsi.ho.results.xbest 
    return hcsi.ho.results
end

function GetInitializedCostates(hcsi::HeuristicsCoStateInitializer)
    return hcsi.λh
end

function SetInitializedCostates!(hsci::HeuristicsCoStateInitializer, λh::AbstractVector)
    hsci.λh .= λh 
    return nothing
end