
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
end

function FSSCoStateInitializer(bvpFunc, ICS, FCS; 
                               costFunc = :WSS, optimizer = :PSO, numParticles = 100,
                               UBs = [40, 40, 40, 2, 2, 2, 2], 
                               LBs = [-40, -40, -40, -2, -2, -2, -2])

    # Check size of ICS and FCS 
    if length(ICS) != 7 && length(FCS) != 7
        throw(ArgumentError("Initial and final boundary condition vectors must be of length 7."))
    end

    # Generate cost function
    if costFunc == :WSS
        cf(x) = cfFSSWSS(x, bvpFunc, ICS, FCS)
    else
        throw(ArgumentError("Cost function type not implemented."))
    end

    # Initialize optimizer and options
    opts = Options(;useParallel = true)
    if optimizer == :PSO 
        prob = Problem(cf, LBs, UBs)
        ho   = PSO(prob; numParticles = numParticles)
    end

    # Instantiate FFS initializer 
    FSSCoStateInitializer{typeof(ho), typeof(opts)}(ho, opts)
end

# Cost functions 
function cfFSSWSS(x, bvpFunc, ICS, FCS)

    # Initialize state/co-state vector
    y0 = @SVector [ICS[1], ICS[2], ICS[3], ICS[4], ICS[5], ICS[6], ICS[7],
                   x[1], x[2], x[3], x[4], x[5], x[6], x[7]]

    # Evaluate bvp function
    yf, timeToFinalTime = bvpFunc(y0)

    # Compute cost
    cost = 10.0*(yf[1] - FCS[1])^2 + 
           10.0*(yf[2] - FCS[2])^2 +
           10.0*(yf[3] - FCS[3])^2 +
           1.0*(yf[4]  - FCS[4])^2 +
           1.0*(yf[5]  - FCS[5])^2 +
           1.0*(yf[6]  - FCS[6])^2 +
           1.0*(yf[14] - FCS[7])^2 + 
           10.0*timeToFinalTime

    # Shift cost if terminated early

    return cost 
end