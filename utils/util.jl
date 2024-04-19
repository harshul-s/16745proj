function rk4(model, ode, x, u, dt)
    # rk4 
    k1 = dt * ode(model, x, u)
    k2 = dt * ode(model, x + k1 / 2, u)
    k3 = dt * ode(model, x + k2 / 2, u)
    k4 = dt * ode(model, x + k3, u)
    x + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
end


function quadrotor_dynamics(model::NamedTuple, x, u)
    # quadrotor dynamics with an MRP for attitude
    # and velocity in the world frame (not body frame)

    r = x[1:3]     # position in world frame 
    q = normalize(Quaternion(x[4:7]))     # quaterion in body frame 
    v = x[8:10]     # n_p_b (MRP) attitude 
    ω = x[11:13]   # angular velocity 

    mass = model.mass
    J = model.J
    gravity = model.gravity
    L = model.L
    kf = model.kf
    km = model.km

    w1 = u[1]
    w2 = u[2]
    w3 = u[3]
    w4 = u[4]

    F1 = max(0, kf * w1)
    F2 = max(0, kf * w2)
    F3 = max(0, kf * w3)
    F4 = max(0, kf * w4)
    F = [0.0, 0.0, F1 + F2 + F3 + F4] #total rotor force in body frame

    M1 = km * w1
    M2 = km * w2
    M3 = km * w3
    M4 = km * w4
    τ = [L * (F2 - F4), L * (F3 - F1), (M1 - M2 + M3 - M4)] #total rotor torque in body frame

    f = mass * gravity + Q * F # forces in world frame

    # this is xdot 
    [
        v
        SVector(0.5 * q * Quaternion(zero(x[1]), ω))
        g + (1 / m) * (q * F + u[5])
        J \ (τ - cross(ω, J * ω))
    ]
end


function combined_dynamics(params, x, u)
    # dynamics for three planar quadrotors, assuming the state is stacked
    # in the following manner: x = [x1;x2;x3]

    N = params.NUM_QUAD
    nx = params.nx
    xdots = [zeros(nx) for i in 1:N]

    for i in 1:N
        x_i_start = (i-1)*13 + 1
        x_i_end = i*13
        u_i_start = (i-1)*5 + 1
        u_i_end = i*5
        x_i = x[x_i_start:x_i_end]
        u_i = u[u_i_start:u_i_end]
        xdots[i] = quadrotor_dynamics(params, x_i, u_i)
    end

    # return stacked dynamics 
    return xdots
end