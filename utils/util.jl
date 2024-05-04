
function rk4(model, ode, x, u, dt)
    # rk4 
    k1 = dt * ode(model, x, u)
    # @show size(k1)
    k2 = dt * ode(model, x + k1 / 2, u)
    k3 = dt * ode(model, x + k2 / 2, u)
    k4 = dt * ode(model, x + k3, u)
    x + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
end
function create_idx(nx, nu, N)
    # This function creates some useful indexing tools for Z 
    # x_i = Z[idx.x[i]]
    # u_i = Z[idx.u[i]]

    # Feel free to use/not use anything here.
    # our Z vector is [x0, u0, x1, u1, …, xN]
    nz = (N - 1) * nu + N * nx # length of Z 
    x = [(i - 1) * (nx + nu) .+ (1:nx) for i = 1:N]
    u = [(i - 1) * (nx + nu) .+ ((nx+1):(nx+nu)) for i = 1:(N-1)]

    # constraint indexing for the (N-1) dynamics constraints when stacked up
    c = [(i - 1) * (nx) .+ (1:nx) for i = 1:(N-1)]
    nc = (N - 1) * nx # (N-1)*nx 

    return (nx=nx, nu=nu, N=N, nz=nz, nc=nc, x=x, u=u, c=c)
end

# function get_obstic