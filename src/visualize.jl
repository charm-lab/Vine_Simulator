using MeshCat, GeometryBasics, CoordinateTransformations, Colors
using GeometryBasics: HyperRectangle, HyperSphere, Vec, Point, Mesh

function setobstacle!(vis,name,o::Circ)
    G = RGB(0, 1, 0)
    setobject!(vis[name], HyperSphere(GeometryBasics.Point(o.x/1000,o.y/1000,0.),o.r/1000), MeshPhongMaterial(color=G))
end

function setobstacle!(vis,name,o::Rect)
    G = RGB(0, 1, 0)
    setobject!(vis[name], HyperRectangle(GeometryBasics.Vec(o.corner_x/1000,o.corner_y/1000,0.), GeometryBasics.Vec(o.width/1000,o.height/1000,.01)), MeshPhongMaterial(color=G))
end

function setobstacle!(vis,name,o::AngledRect)
    G = RGB(0, 1, 0)
    setobject!(vis[name], HyperRectangle(GeometryBasics.Vec(-o.height/2000,-o.width/2000,0.), GeometryBasics.Vec(o.height/1000,o.width/1000,.01)), MeshPhongMaterial(color=G))
    settransform!(vis[name], compose(Translation(o.cen_x/1000, o.cen_y/1000, 0.0),LinearMap(AngleAxis(-o.theta, 0, 0, 1))))
end

function change_units_Z(vine, Z; scale = .001)
    _, N = size(Z)
    return Z .* repeat([scale;scale;1],outer = (2*vine.nb,N))
end

function visualize!(m::Vine,Z,withEnv=true;showIntermediate = false)
    diam = m.diam/1000 # tube diameter
    r = m.r/1000
    links = m.nΦ
    Δt = m.Δt
    N = length(Z[1,:])

    vis = Visualizer()
    open(vis)

    num = 6 # num intermediate points

    # Define colors for convenience
    R = RGB(1, 0, 0.)
    B = RGB(0, 0, 1.)

    # Visualize objects in environment
    if withEnv == true
        for i=1:length(m.env.objects)
            setobstacle!(vis,"obj$i",m.env.objects[i])
        end
    end

    # Create points of robot
    for i = 1:links
        # Proximal endpoint of link
        if i == 1
            setobject!(vis["prox$i"], HyperSphere(GeometryBasics.Point(0.,0.,0.),diam/2), MeshPhongMaterial(color=B))
        else
            setobject!(vis["prox$i"], HyperSphere(GeometryBasics.Point(0.,0.,0.),diam/2), MeshPhongMaterial(color=R))
        end

        # Intermediate points
        if showIntermediate == true
            for j = 1:num
                setobject!(vis["p$i$j"], HyperSphere(GeometryBasics.Point(0.,0.,0.),diam/2), MeshPhongMaterial(color=B))
            end
        end

        # Distal endpoint of link
        setobject!(vis["dist$i"], HyperSphere(GeometryBasics.Point(0.,0.,0.),diam/2), MeshPhongMaterial(color=R))
    end

    # Generate animation
    anim = MeshCat.Animation(Int(1/Δt))
    for k = 1:N
        atframe(anim, k-1) do
            for i = 1:links
                # Proximal endpoint of link
                x1 = Z[6*i-5,k]
                y1 = Z[6*i-4,k]
                θ1 = Z[6*i-3,k]
                X1 = x1 - r*cos(θ1)
                Y1 = y1 - r*sin(θ1)

                # Distal endpoint of link
                x2 = Z[6*i-2,k]
                y2 = Z[6*i-1,k]
                θ2 = Z[6*i,k]
                X2 = x2 + r*cos(θ2)
                Y2 = y2 + r*sin(θ2)

                # Set point locations
                settransform!(vis["prox$i"], Translation(X1, Y1, 0))
                if showIntermediate == true
                    for j = 1:num
                        ratio = j/(num+1)
                        settransform!(vis["p$i$j"], Translation(ratio*X1 + (1-ratio)*X2, ratio*Y1 + (1-ratio)*Y2, 0))
                    end
                end
                settransform!(vis["dist$i"], Translation(X2, Y2, 0))
            end
        end
    end

    setanimation!(vis, anim)
end
