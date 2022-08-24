
# FFSCoStateInitializer : (F)orward (S)ingle (S)hooting CoState Initializer
# This struct consaints information and has associated methods for 
# initializing the initial time co-state variables for solving 
# indirect trajectory optimization problems using a forward single
# shooting based approach.
#
# Note: Currently, implementation assumes the trajectory optimization 
# problem is for a single spacecraft with a 6-element state vector + mass,
# and therefore, 7 co-state variables. Other assumptions also made in cost
# function(s) that should likely be addressed

mutable struct FSSCoStateInitializer{HOT,HOOT} <: HeuristicsCoStateInitializer{HOT}
    # Heuristic optimizer 
    ho::HOT 

    # Heuristic optimizer options 
    hoOpts::HOOT

    # Initialized co-state vector
    λh::Vector{Float64}
end

function FSSCoStateInitializer(bvpFunc, tspan, ICS, FCS; 
                               costFunc = :WSS, optimizer = :PSO, numParticles = 100, numSwarms = 4,
                               initMethod = :Uniform,
                               UBs = [100, 100, 100, 50, 50, 50, 50], 
                               LBs = [-100, -100, -100, -50, -50, -50, -50],
                               iUBs = nothing, 
                               iLBs = nothing,
                               weights = nothing,
                               display = true,
                               displayInterval = 1,
                               maxIters = 1000,
                               funcTol = 1e-6,
                               MFD = 2e-6,
                               minNeighborhoodFraction = 0.25,
                               maxStallIters = 25,
                               maxStallTime = 500,
                               maxTime = 1800,
                               useParallel = true
                               )

    # Check size of ICS and FCS 
    if length(ICS) != 7 && length(FCS) != 7
        throw(ArgumentError("Initial and final boundary condition vectors must be of length 7."))
    end

    # Generate cost function and heuristic optimization problem
    if costFunc == :WSS # Weighted Sum of Squares
        
        # Check length of weights vector
        if weights === nothing
            weights = [10, 10, 10, 1, 1, 1, 1]
        elseif length(weights) != 7
            throw(ArgumentError("Weights vector must be of length 7 to use weighted sum of squares cost function."))
        end

        # Instantiate problem
        prob = Problem(x->cfFSSWSS(x, y->bvpFunc(y,tspan), ICS, FCS, weights), LBs, UBs)

        # Test cost function
        try
        xtest = @SVector [UBs[1] - LBs[1], UBs[2] - LBs[2], UBs[3] - LBs[3], 
                          UBs[4] - LBs[4], UBs[5] - LBs[5], UBs[6] - LBs[6], UBs[7] - LBs[7]]
        prob.f(xtest)
        catch
            throw(ErrorException("Cost function test failed! Likely cause of this error is selecting the incorrect initialization integrator for Weighted Sum of Squares cost without integral cost (:WSS)."))
        end

    elseif costFunc == :WSSWC

        # Check length of weights vector
        if weights === nothing 
            weights = [10, 10, 10, 1, 1, 1, 1, 10]
        elseif length(weights) != 8
            throw(ArgumentError("Weights vector must be of length 8 to use weighted sum of squares w/ integral cost function."))
        end

        # Instantiate problem
        prob = Problem(x->cfFSSWSSWC(x, y->bvpFunc(y,tspan), ICS, FCS, weights), LBs, UBs)

        # Test cost function
        try
        xtest = @SVector [UBs[1] - LBs[1], UBs[2] - LBs[2], UBs[3] - LBs[3], 
                          UBs[4] - LBs[4], UBs[5] - LBs[5], UBs[6] - LBs[6], UBs[7] - LBs[7]]
        prob.f(xtest)
        catch
            throw(ErrorException("Cost function test failed! Likely cause of this error is selecting the incorrect initialization integrator for Weighted Sum of Squares cost with integral cost (:WSSWC)."))
        end
    else
        throw(ArgumentError("Cost function type not implemented."))
    end

    # Initialize optimizer and options
    opts = Options(;display = display, displayInterval = displayInterval,
                    maxIters = maxIters, useParallel = useParallel,
                    maxStallIters = maxStallIters, maxStallTime = maxStallTime,
                    maxTime = maxTime, iUB= iUBs, iLB = iLBs, funcTol = funcTol)

    if optimizer == :PSO 
        ho   = PSO(prob; numParticles = numParticles, initMethod = initMethod, minNeighborFrac = minNeighborhoodFraction)
    elseif optimizer == :MS_PSO 
        ho   = MS_PSO(prob; numParticlesPerSwarm = numParticles, numSwarms = numSwarms, initMethod = initMethod)
    elseif optimizer == :MPSO
        ho   = MPSO(prob; numParticles = numParticles, MFD = MFD)
    else
        throw(ArgumentError("Optimizer is not implemented"))
    end

    # Instantiate FFS initializer 
    FSSCoStateInitializer{typeof(ho), typeof(opts)}(ho, opts, Vector{Float64}(undef, 7))
end

# Cost functions 
function cfFSSWSS(x, bvpFunc, ICS, FCS, ws)
    # Initialize state/co-state vector
    y0 = @SVector [ICS[1], ICS[2], ICS[3], ICS[4], ICS[5], ICS[6], ICS[7],
                   x[1], x[2], x[3], x[4], x[5], x[6], x[7]]

    # Evaluate bvp function
    yf, timeToFinalTime = bvpFunc(y0)

    # Compute cost
    cost = ws[1]*(yf[1] - FCS[1])^2 + 
           ws[2]*(yf[2] - FCS[2])^2 +
           ws[3]*(yf[3] - FCS[3])^2 +
           ws[4]*(yf[4]  - FCS[4])^2 +
           ws[5]*(yf[5]  - FCS[5])^2 +
           ws[6]*(yf[6]  - FCS[6])^2 +
           ws[7]*(yf[14] - FCS[7])^2 + 
           100.0*timeToFinalTime

    # Shift cost if terminated early
    return cost 
end

function cfFSSWSSWC(x, bvpFunc, ICS, FCS, ws)
    # Initialize state/co-state vector
    y0 = @SVector [ICS[1], ICS[2], ICS[3], ICS[4], ICS[5], ICS[6], ICS[7],
                   x[1], x[2], x[3], x[4], x[5], x[6], x[7], 0.0]

    # Evaluate bvp function
    yf, timeToFinalTime = bvpFunc(y0)

    # Compute cost
    cost = ws[1]*(yf[1]  - FCS[1])^2 + 
           ws[2]*(yf[2]  - FCS[2])^2 +
           ws[3]*(yf[3]  - FCS[3])^2 +
           ws[4]*(yf[4]  - FCS[4])^2 +
           ws[5]*(yf[5]  - FCS[5])^2 +
           ws[6]*(yf[6]  - FCS[6])^2 +
           ws[7]*(yf[14] - FCS[7])^2 + 
           ws[8]*yf[15]^2 +
           timeToFinalTime

    # Shift cost if terminated early
    return cost 
end