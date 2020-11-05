function c!(c,q)
    nc = length(c)

    # base pin
    c[1] = q[1] - vine.d*cos(q[3])
    c[2] = q[2] - vine.d*sin(q[3])

    # pin elements
    for i=2:Int(nc/4)
        x1 = q[6*i-8]
        y1 = q[6*i-7]
        Θ1 = q[6*i-6]
        x2 = q[6*i-5]
        y2 = q[6*i-4]
        Θ2 = q[6*i-3]
        c[4*i-3] = x2-x1 - vine.d*(cos(Θ1)+cos(Θ2))
        c[4*i-2] = y2-y1 - vine.d*(sin(Θ1)+sin(Θ2))
    end

    # prismatic elements
    for i=1:Int(nc/4)
        x1 = q[6*i-5]
        y1 = q[6*i-4]
        Θ1 = q[6*i-3]
        x2 = q[6*i-2]
        y2 = q[6*i-1]
        Θ2 = q[6*i]
        c[4*i-1:4*i] .= [-sin(Θ1) cos(Θ1); -sin(Θ2) cos(Θ2)]*[x2-x1;y2-y1]
    end
end
