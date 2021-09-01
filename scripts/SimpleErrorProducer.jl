using IndirectTrajOpt
using DifferentialEquations
using StaticArrays
using LinearAlgebra

# Initialize BVP function
ps = initCR3BPIndirectParams("Low Thrust 10 CR3BP")
tspan = (0.0, 8.6404*24*3600/ps.crp.TU)

# Error producing initial conditions
y0 =   @SVector [-0.0194885115,
                 -0.0160334798,
                 0.0,
                 8.9188819237,
                 -4.081936888,
                 0.0,
                 1.0,
                 -18.20557373360214,
                 -16.683730836929307,
                 40.0,
                 0.04174822756326783,
                 -0.005163201648349852,
                 -0.03929039276763477,
                 0.12683528760241147]

# Set thrust type
cSc = ps.sp.isp*9.81*ps.crp.TU / (ps.crp.LU*1000.0)
λv = norm(view(y0, 11:13))
S = IndirectTrajOpt.computeS(y0, λv, cSc)
if S > ps.ϵ; ps.utype = 0 
elseif S < -ps.ϵ; ps.utype = 2
else; ps.utype = 1; end

cb = VectorContinuousCallback(
            IndirectTrajOpt.cr3bpEomsCondition,
            IndirectTrajOpt.cr3bpEomsAffect!,
            IndirectTrajOpt.cr3bpEomsAffect!, 4;
            idxs = nothing,
            rootfind = DiffEqBase.LeftRootFind,
            interp_points = 10,
            abstol = 1e-14,
            reltol = 0.0,
            save_positions = (true, true))
ff = ODEFunction{false}(IndirectTrajOpt.cr3bpEomIndirect)
prob = ODEProblem(ff, y0, tspan, ps; callback=cb)

sol = solve(prob, Vern9(),reltol=1e-14,abstol=1e-14)