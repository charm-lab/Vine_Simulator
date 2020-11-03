include("kinematics.jl")

# indices for convenience
q_idx = [(t-1)*(2*vine.nq) .+ (1:vine.nq) for t = 1:T]
v_idx = [(t-1)*(2*vine.nq) + vine.nq .+ (1:vine.nq) for t = 1:T]
r_idx = [2*vine.nq*T + (t-1)*vine.nq .+ (1:vine.nq) for t = 1:T-1]
s_idx = [2*vine.nq*T + vine.nq*(T-1) + (t-1)*vine.nq .+ (1:vine.nq) for t = 1:T-1]

# hyperparameters
αx = .5
αy = .05
βr = 5e4
βs = 5e4

# objective function
function obj(z)
    # unpack decision variables
    q = [z[q_idx[t]] for t = 1:T]
    r = [z[r_idx[t]] for t = 1:T-1]
    s = [z[s_idx[t]] for t = 1:T-1]

    vid_idx = 1
    cost = 0.0
    for t = 1:T
        if mod(t, 3) == 1
            xcoords, ycoords = kinematics(q[t], vine.r, A)
            cost += αx*(xcoords - x_vid[vid_idx])'*(xcoords - x_vid[vid_idx])
            cost += αy*(ycoords - y_vid[vid_idx])'*(ycoords - y_vid[vid_idx])
            # println((xcoords - x_vid[vid_idx])'*(xcoords - x_vid[vid_idx]))
            vid_idx += 1
        end

        t==T && continue
        cost += βr*r[t]'*r[t] + βs*s[t]'*s[t]
    end

    return cost
end
