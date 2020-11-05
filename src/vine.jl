using LinearAlgebra
using StaticArrays
include("../src/objects.jl")

mutable struct Vine{T,S,D,V}
    θ0::T       # initial heading
    diam::T     # tube diameter
    d::T        # distance from body CoM to body endpoint
    k::S        # stiffness and damping matrix

    M::D        # mass matrix
    MInv::D     # inverted mass matrix
    R::S        # mapping from torques to maximal coordinates
    Fext::V     # bending torques
    J::S        # mapping from joint impulses to maximal coordinates
    c::V        # joint constraint c(q)
    λ::V        # joint constraint impulse
    L::S        # mapping from contact impulses to maximal coordiantes
    Φ::V        # signed distance function Φ(q)
    n::V        # contact impulses
    G::S        # mapping from growth impulses to maximal coordinates
    g::V        # growth rate of prismatic joints g(q)
    w::V        # growth impulses w

    Δt::T       # length of timestep

    nb::Int     # number of rigid bodies
    nq::Int     # size of configuration vector
    nc::Int     # size of joint constraint
    nu::Int     # size of growth constraint
    ne::Int     # number of objects in environment
    nΦ::Int     # size of contact constraint

    env::Env    # list of objects in environment
end

function create_vine(links, d; stiffness=500000, damping=30,
                    m_b=.001, # mass
                    J_b=200, # inertia
                    θ0 = 0., # initial heading
                    diam = 24.0, # tube diameter in mm
                    Δt = .01, # timestep
                    objects = Env([Circ(470,25,60)]) # environment
                    )
    # Dimensions
    nb = 2*links # number of rigid bodies
    nq = 3*nb # configuration dimension
    nc = 2*nb # number of primary constraints
    nu = links # control dim
    ne = length(objects.objects) # num objects in environment
    nΦ = links

    # Stiffness and damping matrices
    k = Matrix([Diagonal(stiffness*ones(links)) Diagonal(damping*ones(links))]) # rotaional spring

    # Mass matrix
    M = Diagonal(SVector{3*nb}( repeat([m_b; m_b; J_b]; outer = [nb])))
    MInv = Diagonal(SVector{3*nb}( repeat([1/m_b; 1/m_b; 1/J_b]; outer = [nb])))

    # Preallocate matrices
    Fext = zeros(nq)
    J = zeros(nc,nq)
    c = zeros(nc)
    λ = zeros(nc)
    L = zeros(nΦ*ne,nq)
    Φ = zeros(nΦ*ne)
    n = zeros(nΦ*ne)
    G = zeros(nu,2*nq)
    g = zeros(links)
    w = zeros(links)

    # Mapping torque to maximal coordinates
    R = zeros(6*nu,nu)
    R[3,1] = 1
    for i = 1:links-1
        R[6*i,i+1] = -1
        R[6*i+3,i+1] = 1
    end

    return Vine(θ0,diam,d,k,M,MInv,R,Fext,J,c,λ,L,Φ,n,G,g,w,Δt,nb,nq,nc,nu,ne,nΦ,objects)
end
