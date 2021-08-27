
abstract type HeuristicsCoStateInitializer{HOT <: Heuristics.Optimizer} end

function initialize!(hcsi::HeuristicsCoStateInitializer)
    optimize!(hcsi.ho, hcsi.hoOpts)
    return hcsi.ho.result
end