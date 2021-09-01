using IndirectTrajOpt
using IndirectCoStateInit

# Initialize BVP function
ps = initCR3BPIndirectParams("Low Thrust CR3BP")
tspan = (0.0, 8.6404*24*3600/ps.crp.TU)
bvpFunc(y0) = cr3bpOptIntegrate(y0, tspan, ps, copyParams = true, termCallbacks = true)

# Set initial and final conditions 
ics = [-0.0194885115, -0.0160334798, 0.0,
       8.9188819237, -4.081936888, 0.0, 1.0]
fcs = [0.8233851820, 0.0, -0.02225563,
       0.0, 0.134184103, 0.0, 0.0]

# Initialize FFS Initializer 
csInitializer = FSSCoStateInitializer(bvpFunc, ics, fcs; 
                                      numParticles = 1000,
                                      displayInterval = 5);

# Initialize co-states 
initialize!(csInitializer)