function rk4(model,ode,x,u,dt)
    # rk4 
    k1 = dt*ode(model,x, u)
    k2 = dt*ode(model,x + k1/2, u)
    k3 = dt*ode(model,x + k2/2, u)
    k4 = dt*ode(model,x + k3, u)
    x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
end


function single_quad_dynamics(params, x, u)
    # planar quadrotor dynamics for a single quadrotor 

    # unpack state
    px, pz, θ, vx, vz, ω = x

    xdot = [
        vx,
        vz,
        ω,
        (1 / params.mass) * (u[1] + u[2]) * sin(θ),
        (1 / params.mass) * (u[1] + u[2]) * cos(θ) - params.g,
        (params.ℓ / (2 * params.J)) * (u[2] - u[1])
    ]

    return xdot
end


function combined_dynamics(params, x, u)
    # dynamics for three planar quadrotors, assuming the state is stacked
    # in the following manner: x = [x1;x2;x3]

    # NOTE: you would only need to use this if you chose option 2 where 
    # you optimize over all three trajectories simultaneously 

    # quadrotor 1 
    x1 = x[1:6]
    u1 = u[1:2]
    xdot1 = single_quad_dynamics(params, x1, u1)

    # quadrotor 2
    x2 = x[(1:6).+6]
    u2 = u[(1:2).+2]
    xdot2 = single_quad_dynamics(params, x2, u2)

    # quadrotor 3
    x3 = x[(1:6).+12]
    u3 = u[(1:2).+4]
    xdot3 = single_quad_dynamics(params, x3, u3)

    # return stacked dynamics 
    return [xdot1; xdot2; xdot3]
end