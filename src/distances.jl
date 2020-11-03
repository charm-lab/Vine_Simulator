using LinearAlgebra

function distToSeg(p,a,b,normal)
    n = b - a;

    pa = a - p;
    c = dot(n, pa);
    # Closest point is a
    if ( c > 0.0 )
        return norm(pa)
    end

    bp = p - b;
    d = dot(n, bp);
    # Closest point is b
    if ( d > 0.0 )
        return norm(bp)
    end

    return dot(bp, normal)
end

function distToPoly(polygon,p)
    numsides = length(polygon[1,:])
    distances = zeros(eltype(p), numsides)

    for i = 1:numsides
        distances[i] = distToSeg(p, polygon[1:2,i], polygon[3:4,i], polygon[5:6,i])
    end

    if sum(distances .>= 0) == 0 # all values negative
        Φ = maximum(distances)
    else
        Φ = minimum(filter((x) -> x >= 0, distances))
    end
    nhatx = polygon[5,findfirst(x->x==Φ, distances)]
    nhaty = polygon[6,findfirst(x->x==Φ, distances)]
    return (Φ,nhatx,nhaty)
end
