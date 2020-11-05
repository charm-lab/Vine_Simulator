import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
using Plots
using JLD2

# Load model
include("src/vine.jl")
include("src/dynamics.jl")

# Load fitting problem data
@load joinpath(@__DIR__,"fitting_problem_data.jld2") vine A U n_opt z0 x_vid y_vid
T = length(U[1,:]) # time horizon

# Load optimization functions
include("src/fitting/ipopt.jl")
include("src/fitting/objective.jl")
include("src/fitting/constraints.jl")

# Set up problem
N = length(z0)# number of decision variables
M = 2*vine.nq*(T-1) + (vine.nc + 1 + vine.nu)*T # number of constraints
prob = ProblemIpopt(N, M) # set up optimization problem for Ipopt

# Solve
z_sol = solve(z0, prob, tol=1., d_tol=1., c_tol=1.0e-3, max_iter = 1)

# Extract trajectory
Z_sol = reshape(z_sol[1:2*vine.nq*T], 2*vine.nq, T)
x_kin_sol = [kinematics(Z_sol[1:vine.nq,t], vine.d, A)[1] for t=1:T]
y_kin_sol = [kinematics(Z_sol[1:vine.nq,t], vine.d, A)[2] for t=1:T]

# Compare tip trajectories
T_vid = length(x_vid)
T = length(x_kin_sol)
vid_t = [i for i=1:3:T]
kin_t = [i for i=1:T]
vid_x = [x_vid[i][end] for i=1:T_vid]
kin_x = [x_kin_sol[i][end] for i=1:22]
plot(vid_t, vid_x)
plot!(kin_t, kin_x, xlabel="Timestep", ylabel="x-coordinate (mm)")
