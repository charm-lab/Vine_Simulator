using Plots, LinearAlgebra

function create_A(sv, sm)
    A = zeros(sv, sm)
    for i = 0:sv-2
        z = i * sm/(sv-1) + 1
        row = i + 1
        col = floor(Int, z)
        A[row, col] = z - col
        if col > 1
            A[row, col-1] = 1 - z + col
        end
        A[end, end] = 1
    end
    return A
end

# input: state vector
# output: xy vector
function kinematics(q, d, A)
    nq = length(q)
    x = [q[i] for i = 1:3:nq]
    y = [q[i] for i = 2:3:nq]
    θ = [q[i] for i = 3:3:nq]

    x_pin = x + d*cos.(θ)
    y_pin = y + d*sin.(θ)

    return A*x_pin, A*y_pin
end
