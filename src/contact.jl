include("distances.jl")

# Compute L and Φ in place
function calcLΦ!(m::Vine, q, objects)
    L = m.L
    Φ = m.Φ
    for i=1:m.ne
        calcLΦ!(L, Φ, m, q, objects[i], i)
    end
    return L, Φ
end

# Return L and Φ
function calcLΦ(m::Vine, q, objects)
    L = zeros(eltype(q), m.nΦ*m.ne, m.nq)
    Φ = zeros(eltype(q), m.nΦ*m.ne)
    for i=1:m.ne
        calcLΦ!(L, Φ, m, q, objects[i], i)
    end
    return L, Φ
end

# Compute L and Φ for a circular object
function calcLΦ!(L_all, Φ_all, m::Vine, q, o::Circ, j::Int)
    nΦ = m.nΦ
    Φ = view(Φ_all, (j-1)*nΦ+1:j*nΦ)
    L = view(L_all, (j-1)*nΦ+1:j*nΦ, :)
    d = m.d
    for i=1:nΦ
        x = q[6*i-2]
        y = q[6*i-1]
        Θ = q[6*i]
        dist = sqrt((x + d*cos(Θ) - o.x)^2 + (y + d*sin(Θ) - o.y)^2)
        Φ[i] = dist - (o.r+m.diam/2)
        nhatx = (x + d*cos(Θ) - o.x)/dist
        nhaty = (y + d*sin(Θ) - o.y)/dist

        L[i,6*i-2] = nhatx
        L[i,6*i-1] = nhaty
        L[i,6*i] = -nhatx*d*sin(Θ) + nhaty*d*cos(Θ)
    end
end

# Compute L and Φ for a rectangle object
function calcLΦ!(L_all, Φ_all, m::Vine, q, o::Rect, j::Int)
    d = m.d
    nΦ = m.nΦ
    Φ = view(Φ_all, (j-1)*nΦ+1:j*nΦ)
    L = view(L_all, (j-1)*nΦ+1:j*nΦ,:)
    for i=1:nΦ
        x = q[6*i-2]
        y = q[6*i-1]
        Θ = q[6*i]

        Φ[i],nhatx,nhaty = distToPoly(o.polygon,[x+d*cos(Θ);y+d*sin(Θ)])
        Φ[i] = Φ[i] - m.diam/2

        L[i,6*i-2] = nhatx
        L[i,6*i-1] = nhaty
        L[i,6*i] = -nhatx*d*sin(Θ) + nhaty*d*cos(Θ)
    end
end

# Compute L and Φ for an angled rectangle object
function calcLΦ!(L_all, Φ_all, m::Vine, q, o::AngledRect, j::Int)
    d = m.d
    nΦ = m.nΦ
    Φ = view(Φ_all,(j-1)*nΦ+1:j*nΦ)
    L = view(L_all,(j-1)*nΦ+1:j*nΦ,:)
    for i=1:nΦ
        x = q[6*i-2]
        y = q[6*i-1]
        Θ = q[6*i]

        Φ[i],nhatx,nhaty = distToPoly(o.polygon,[x+d*cos(Θ);y+d*sin(Θ)])
        Φ[i] = Φ[i] - m.diam/2
        # println(Φ[i])
        L[i,6*i-2] = nhatx
        L[i,6*i-1] = nhaty
        L[i,6*i] = -nhatx*d*sin(Θ) + nhaty*d*cos(Θ)
    end
end
