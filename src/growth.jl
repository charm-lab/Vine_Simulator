function W!(g,qv)
    links = length(g)
    nq = Int(length(qv)/2)
    q = qv[1:nq]
    v = qv[nq+1:2*nq]
    for i = 0:links-1
        x1 = q[6*i+1]
        y1 = q[6*i+2]
        x2 = q[6*i+4]
        y2 = q[6*i+5]
        vx1 = v[6*i+1]
        vy1 = v[6*i+2]
        vx2 = v[6*i+4]
        vy2 = v[6*i+5]
        g[i+1] = ((x2-x1)*(vx2-vx1) + (y2-y1)*(vy2-vy1))/sqrt((x2-x1)^2+(y2-y1)^2);
    end
end
