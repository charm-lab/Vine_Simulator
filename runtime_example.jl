import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
using BenchmarkTools
using Plots

# Load model
include("src/vine.jl")
include("src/dynamics.jl")
include("src/visualize.jl")

# Compute runtime of simulator and visualize trajectory
function run_sim(links, N=400, growth_rate=40; vis=false)
    global vine

    # Create vine model
    println("Setting up model with $(2*links) rigid bodies...")
    d = 100/links # link radius in mm
    vine = create_vine(links, d)
    nq = vine.nq

    # Set starting position
    q = zeros(nq)
    q[1:3:nq] .= [d*i for i=1:2:2*vine.nb]

    # Set initial state
    q = SVector{nq}(q)
    v = @SVector zeros(nq)
    z0 = [q;v]

    # Set up matrices for simulation
    Z = zeros(2 * vine.nq, N)
    U = growth_rate/links * ones(vine.nu, N - 1)

    # Simulate trajectory
    b = @benchmark quick_rollout!($vine, $Z, $z0, $U) samples=1 evals=1
    println("Trajectory computed in ", round(b.times[1]/1e9, digits=2), " seconds")

    # Visualize
    if vis
        println("Visualizing trajectory...")
        Z_meters = change_units_Z(vine, Z) # convert units to meters
        visualize!(vine, Z_meters, showIntermediate=true)
    end

    return b.times[1]/(N*1e6)
end

# Batch run
b_all = []
num_links = collect(10:5:50)
for nl in num_links
    push!(b_all, run_sim(nl))
end

# Single run with visualization
run_sim(15, vis=true)

# Plot runtime scaling
scatter(2*num_links, b_all, xlabel = "Number of Rigid Bodies in Model", ylabel="Average Runtime per Timestep (ms)", label="Simulator")
plot!(2*[num_links[1], num_links[end]], [10, 10], label="Real-time")
