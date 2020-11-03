import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()

include("src/vine.jl")
include("src/dynamics.jl")
include("src/visualize.jl")

# Box maze environment
ipm = 39.3701/1000 # inches per mm
b1 = createRect(5.5/ipm,-5/ipm,4/ipm,7/ipm)
b2 = createRect(4/ipm,-17/ipm,7/ipm,5/ipm)
b3 = createRect(13.5/ipm,-17/ipm,8/ipm,6/ipm)
b4 = createRect(13.5/ipm,-30/ipm,4/ipm,5/ipm)
b5 = createRect(20.5/ipm,-30/ipm,4/ipm,5/ipm)
b6 = createRect(13.5/ipm,-30/ipm,10/ipm,1/ipm)
objects = Env([b1 b2 b3 b4 b5 b6])

# Vine model
links = 25 # number of pairs of bodies
r = 3.0 # distance from CoM to body end (mm)
vine = create_vine(links, r, stiffness=30000., damping=10.,
                    m_b=.002, # mass
                    J_b=5, # inertia
                    θ0 = -52 * pi / 180, # initial heading
                    diam = 15.0, # tube diameter in mm
                    Δt = 1/90, # timestep
                    objects = objects) # environment

# Vine dimensions
nq = vine.nq
nb = vine.nb
nu = vine.nu

# Set initial state
x = vine.θ0 * ones(nq)
x[1:2] = [vine.r * cos(vine.θ0); vine.r * sin(vine.θ0)]
for i = 1:nb-1
    x[3*i+1:3*i+2] =
        (2 * i + 1) * [vine.r * cos(vine.θ0); vine.r * sin(vine.θ0)]
end

q0 = SVector{nq}(x)
x0 = SVector{2*nq}([x; zeros(nq)])

# Set up trajectory matrices
N = round(Int, 3/vine.Δt) # 3-second trajectory
Z = zeros(2 * nq, N)
U = 205/links * ones(nu, N - 1)

# Generate and visualize trajectory
quick_rollout!(vine, Z, x0, U)
Z_meters = change_units_Z(vine, Z)
visualize!(vine, Z_meters, showIntermediate=true)
