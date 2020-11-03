function C!(cq,q)
    nc = length(cq)

    # base pin
    cq[1] = q[1] - vine.r*cos(q[3])
    cq[2] = q[2] - vine.r*sin(q[3])

    # pin elements
    for i=2:Int(nc/4)
        x1 = q[6*i-8]
        y1 = q[6*i-7]
        Θ1 = q[6*i-6]
        x2 = q[6*i-5]
        y2 = q[6*i-4]
        Θ2 = q[6*i-3]
        cq[4*i-3] = x2-x1 - vine.r*(cos(Θ1)+cos(Θ2))
        cq[4*i-2] = y2-y1 - vine.r*(sin(Θ1)+sin(Θ2))
    end

    # prismatic elements
    for i=1:Int(nc/4)
        x1 = q[6*i-5]
        y1 = q[6*i-4]
        Θ1 = q[6*i-3]
        x2 = q[6*i-2]
        y2 = q[6*i-1]
        Θ2 = q[6*i]
        cq[4*i-1:4*i] .= [-sin(Θ1) cos(Θ1); -sin(Θ2) cos(Θ2)]*[x2-x1;y2-y1]
    end
end
