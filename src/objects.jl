using LinearAlgebra, Rotations

struct Circ{T}
    x::T
    y::T
    r::T
end

struct Rect{M,T}
    polygon::M
    corner_x::T
    corner_y::T
    width::T
    height::T
end

struct AngledRect{M,T}
    polygon::M
    cen_x::T
    cen_y::T
    theta::T
    width::T
    height::T
end

struct Env{O}
    objects::O
end

# Create Angled Rectangle by specifying the center, angle, width, and height
function createAngledRect(cen_x,cen_y,theta,w,h)
    c = [cen_x;cen_y]
    v1 = [cos(-theta); sin(-theta)]
    v2 = [cos(-theta-pi/2); sin(-theta-pi/2)]

    p1 = c + .5*(h*v1 + w*v2)
    p2 = c + .5*(h*v1 - w*v2)
    p3 = c + .5*(-h*v1 - w*v2)
    p4 = c + .5*(-h*v1 + w*v2)

    polygon = [p1[1] p2[1] p3[1] p4[1];
               p1[2] p2[2] p3[2] p4[2];
               p2[1] p3[1] p4[1] p1[1];
               p2[2] p3[2] p4[2] p1[2];
               v1    -v2   -v1   v2]

    return AngledRect(polygon,cen_x,cen_y,theta,w,h)
end

# Create Angled Rectangle by specifying the midline endpoints and width
function createAngledRectFromEndpoints(x1,y1,x2,y2,w)
    cen_x = .5*(x1+x2)
    cen_y = .5*(y1+y2)
    theta = -atan((y2-y1)/(x2-x1))
    h = norm([x1;y1]-[x2;y2])
    return createAngledRect(cen_x,cen_y,theta,w,h)
end

# Create Rectangle by specifying the lower left corner, width, and height
function createRect(x1,y1,w,h)
    x2 = x1
    y2 = y1 + h
    x3 = x1 + w
    y3 = y1 + h
    x4 = x1 + w
    y4 = y1

    polygon = [x1 x2 x3 x4;
               y1 y2 y3 y4;
               x2 x3 x4 x1;
               y2 y3 y4 y1;
               -1. 0. 1. 0.;
               0. 1. 0. -1.]

    return Rect(polygon,x1,y1,w,h)
end

# Create Rectangle by specifying the four corners
# Input: p1 is bottom left corner, p2/p3/p4 goes counterclockwise
function createRectFromPoints(p1,p2,p3,p4)
    n1 = reverse((p1-p2)/norm(p1-p2))
    n2 = reverse((p2-p3)/norm(p2-p3))
    n3 = reverse((p3-p4)/norm(p3-p4))
    n4 = reverse((p4-p1)/norm(p4-p1))

    polygon = [p1[1] p2[1] p3[1] p4[1];
               p1[2] p2[2] p3[2] p4[2];
               p2[1] p3[1] p4[1] p1[1];
               p2[2] p3[2] p4[2] p1[2];
               n1    n2    n3    n4]

    return Rect(polygon,p1[1],p1[2],p3[1]-p1[1],p2[2]-p1[2])
end
