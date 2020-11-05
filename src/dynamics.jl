using ForwardDiff
include("joints.jl")
include("contact.jl")
include("growth.jl")
include("solver.jl")

function dynamics(m::Vine,x, u)
    nb = m.nb
    nc = m.nc
    nq = m.nq
    q = view(x,1:nq)
    qd = view(x,nq+1:2*nq)
    J = m.J
    λ = m.λ
    F = m.Fext
    MI = m.MInv
    r = m.r
    Δt = m.Δt
    L = m.L
    n = m.n

    # assemble vector of angles and angular velocities for calculating fExt
    θall = [i == 3 ?  q[i] - m.θ0 : q[i] - q[i-3] for i=3:6:length(q)]
    θdall = [i == 3 ? qd[i] : qd[i] - qd[i-3] for i=3:6:length(qd)]

    # external forces
    F.= -m.R * m.k * [θall; θdall]

    # joint constraints
    c!(m.c,q)
    ForwardDiff.jacobian!(J, c!, ones(nc), q)

    # contact constraints
    calcLΦ!(m,q,m.env.objects)

    # growth constraints
    g!(m.g,x)
    ForwardDiff.jacobian!(m.G, g!, ones(length(m.g)), x)

    # solve convex problem
    vkp1 = solveQP!(m, qd, u)
    qnext = q + Δt*vkp1

    return [qnext; vkp1]
end

function quick_rollout!(m, Z::Array, x0, U)
    Z[:,1] = x0
    for k = 2:length(Z[1,:])
        Z[:,k] = dynamics(m, SVector{2*m.nq}(Z[:,k-1]), SVector{m.nu}(U[:,k-1]))
    end
end
