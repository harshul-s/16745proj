include("quaternions.jl")

quad_params = (
    mass = 0.5,
    J = SMatrix{3,3}(Diagonal([0.0023, 0.0023, 0.004])),
    Jinv=SMatrix{3,3}(Diagonal(1.0./[0.0023, 0.0023, 0.004])),
    gravity = SVector(0,0,-9.81),
    motor_dist = 0.1750,
    kf = 1.0,
    km = 0.0245
)

function load_dynamics(x,u)
    ẋ = zeros(6)
    ẋ[1:3] = x[4:6]
    ẋ[4:6] = u[1:3]
    ẋ[6] -= 9.81 # gravity
    return ẋ
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
    L = model.motor_dist
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

function combined_quadrotor_dynamics(params, x, u)
    # dynamics for three planar quadrotors, assuming the state is stacked
    # in the following manner: x = [x1;x2;x3]

    N = params.NUM_QUAD
    # nx = params.nx
    idx = params.idx
    xdots = [zeros(params.n_lift) for i=1:N]

    for i=1:N
        xi = x[(i-1)*13 + 1:13*i]
        ui = u[(i-1)*5 + 1:5*i]
        xdots[i] = quadrotor_dynamics(params.quad_params, xi, ui)
    end
    # return stacked dynamics 
    return xdots
end


function combined_load_dynamics(params, x, u)
    # u1;u2;u3;u4;u5 --> quadrotors
    # ul1,ul2,ul3,ul4 --> 
    #[u1;u2;u3;u4;u5 u1;u2;u3;u4;u5 u1;u2;u3;u4;u5 u1;u2;u3;u4;u5 ul1,ul2,ul3,ul4]
    idx = params.idx
    L = params.NUM_QUAD
    load_params = params.load
    load_mass = load_params.m
    
    #Get 3D Positions
    r_inds = [(1:3) .+ i for i in (0:13:13*num_lift)]
    r = x[r_inds]
    r_quad = r[1:L]
    r_load = r[end]

    #Get U 
    u_inds = [(1:4) .+ i for i in (0:5:5*num_lift-1)]
    u_quad = u[u_inds]
    u_load = u[(1:L) .+ 5*L]
    s_quad = u[5:5:5*L]

    dir = [(r_load - r)/norm(r_load - r) for r in r_quad]
    lift_control = [[u_quad[i]; s_quad[i]*dir[i]] for i = 1:L]
    u_slack_load = -1.0*sum(dir .* u_load)

    load_inds = (1:6) .+ 13*num_lift
    xdots = load_dynamics(x[load_inds], u_slack_load / load_mass)
    return xdots
end 

function combined_system_dynamics(params,x̄,ū)
    return [combined_quadrotor_dynamics(params,x̄,ū) ; combined_load_dynamics(params,x̄,ū)]
end


# function combined_dynamics(params, x, u)
#     # dynamics for three planar quadrotors, assuming the state is stacked
#     # in the following manner: x = [x1;x2;x3]

#     N = params.NUM_QUAD
#     nx = params.nx
#     xdots = [zeros(nx) for i in 1:N]

#     for i in 1:N
#         x_i_start = (i-1)*13 + 1
#         x_i_end = i*13
#         u_i_start = (i-1)*5 + 1
#         u_i_end = i*5
#         x_i = x[x_i_start:x_i_end]
#         u_i = u[u_i_start:u_i_end]
#         xdots[i] = quadrotor_dynamics(params, x_i, u_i)
#     end

#     # return stacked dynamics 
#     return xdots
# end