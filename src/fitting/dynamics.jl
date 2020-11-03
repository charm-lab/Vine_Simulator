function compute_jacobians(m::Vine, x::Vector, u::Vector, stiffness::Number, damping::Number)
    nc = m.nc
    nq = m.nq
    q = view(x,1:nq)
    qd = view(x,nq+1:2*nq)

    # assemble vector of angles and angular velocities for calculating fExt
    θall = [i == 3 ?  q[i] - m.θ0 : q[i] - q[i-3] for i=3:6:length(q)]
    θdall = [i == 3 ? qd[i] : qd[i] - qd[i-3] for i=3:6:length(qd)]

    # external forces
    F = -m.R * (stiffness*θall + damping*θdall)

    # joint constraints
    cq = zeros(eltype(x), m.nc)
    J = zeros(eltype(x), m.nc, m.nq)
    C!(cq, q)
    ForwardDiff.jacobian!(J, C!, ones(eltype(q), nc), q)

    # contact constraints
    L, Φ = calcL(m,q,m.env.objects)

    # growth constraints
    g = zeros(eltype(x), m.nu)
    G = zeros(eltype(x), size(m.G))
    W!(g,[q;qd])
    ForwardDiff.jacobian!(G, W!, ones(eltype(q), length(m.g)), [q;qd])

    return J, cq, L, Φ, G, g, F
end

function compute_impulses(m::Vine, x::Vector, u::Vector, stiffness::Number, damping::Number, normal_force::Number, M, MInv)
    Δt = m.Δt
    nq = m.nq
    qk = x[1:nq]
    vk = x[1+nq:2*nq]

    J, cq, L_all, Φ_all, G, g, F = compute_jacobians(m, x, u, stiffness, damping)

    # growth
    g_coeff = G[:,1:nq]*m.Δt + G[:,nq+1:2*nq]
    g_con = g - u - G[:,nq+1:2*nq]*vk

    # extract last row of L and Φ
    L = L_all[end,:]'
    Φ = Φ_all[end]

    impulses = F*Δt + L'*normal_force

    J̃ = [J;g_coeff]
    J̃M = [J;g_coeff;L]*MInv
    A = J̃M*J̃'
    b = -([J*vk + cq/Δt;g_coeff*vk + g_con;L*vk + Φ/Δt] + J̃M*impulses)
    λ̃ = A\b

    return J̃'*λ̃  + impulses
end
