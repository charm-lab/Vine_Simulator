using JuMP
using OSQP
using LinearAlgebra

function solveQP!(m::Vine, v, u, verbose = false)
    M = m.M
    MInv = m.MInv
    τ = m.Fext
    Δt = m.Δt
    J = m.J
    c = m.c
    L = m.L
    Φ = m.Φ
    nb = m.nb
    nq = m.nq
    G = m.G
    g = m.g
    g_con = deepcopy(g)
    g_con .= g - u - G[:,nq+1:2*nq]*v
    g_coeff = G[:,1:nq]*Δt + G[:,nq+1:2*nq]

    model = JuMP.Model(OSQP.Optimizer)
    set_silent(model)
    @variable(model, vkp1[1:3*nb])
    @objective(model, Min, .5*vkp1'*M*vkp1 - vkp1'*(M*v+τ*Δt))
    @constraint(model, sdf, Φ/Δt + L*vkp1 .>= 0)
    @constraint(model, joint, c/Δt + J*vkp1 .== 0)
    @constraint(model, growth, g_con + g_coeff*vkp1 .== 0)

    JuMP.optimize!(model)
    m.λ .= dual.(joint)
    m.n .= dual.(sdf)
    m.w .= dual.(growth)

    if verbose
        println("Objective value: ", JuMP.objective_value(model))
    end

    return JuMP.value.(vkp1)
end
