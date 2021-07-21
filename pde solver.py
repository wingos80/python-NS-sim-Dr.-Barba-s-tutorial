import numpy as np, timeit
from matplotlib import pyplot as plt, cm
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

toc = timeit.default_timer()

xmax, ymax = 4, 2

##variable declarations
nx = 360
ny = 180
nt = 10
nit = 50
c = 1
dx = xmax / (nx - 1)
dy = ymax / (ny - 1)
x = np.linspace(0, xmax, nx)
y = np.linspace(0, ymax, ny)
X, Y = np.meshgrid(x, y)

##physical variables
rho = 1
nu = .1
F = 0
dt = .0003

# Initialization
u = np.zeros((ny,nx))*-1
un = np.zeros((ny,nx))

v = np.zeros((ny,nx))
vn = np.zeros((ny,nx))

p  = np.zeros((ny, nx))
b = np.zeros((ny, nx))


def plot2D(x, y, p):
    fig = plt.figure(figsize=(11, 7), dpi=100)
    ax = fig.gca(projection='3d')
    X, Y = np.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.viridis,
            linewidth=0, antialiased=False)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

    plt.show()

def build_up_b(rho, dx, dy, u, v):
    b = np.zeros_like(u)
    b[1:-1, 1:-1] = rho*(((u[1:-1, 2:] - u[1:-1, :-2]) / (2 * dx)+(v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * dy))/dt
                   -(u[1:-1, 2:] - u[1:-1, :-2]) / (2 * dx)**2
                   -2*(u[1:-1, 2:] - u[1:-1, :-2]) / (2 * dy)*(v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * dx)
                   -(v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * dy)**2)

    # # one sided ns @ x = 2
    # b[1:-1, -1] = rho*(1/dt*((u[1:-1,-1]-u[1:-1, -2])/dx + (v[2:,-1] - v[:-2,-1])/(2*dy))-
    #                    ((u[1:-1, -1]-u[1:-1,-2])/dx)**2-2*(u[2:,-1]-u[:-2,-1])*(u[1:-1,-1]-u[1:-1,-2])/(2*dy*dx))

    # b[1:-1, -1] =0

    # b[1:-1, 0] = 1
    # # Periodic BC Pressure @ x = 2
    # b[1:-1, -1] = (rho * (1 / dt * ((u[1:-1, 0] - u[1:-1,-2]) / (2 * dx) +
    #                                 (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) -
    #                       ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx))**2 -
    #                       2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) *
    #                            (v[1:-1, 0] - v[1:-1, -2]) / (2 * dx)) -
    #                       ((v[2:, -1] - v[0:-2, -1]) / (2 * dy))**2))
    #
    # # Periodic BC Pressure @ x = 0
    # b[1:-1, 0] = (rho * (1 / dt * ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx) +
    #                                (v[2:, 0] - v[0:-2, 0]) / (2 * dy)) -
    #                      ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx))**2 -
    #                      2 * ((u[2:, 0] - u[0:-2, 0]) / (2 * dy) *
    #                           (v[1:-1, 1] - v[1:-1, -1]) / (2 * dx))-
    #                      ((v[2:, 0] - v[0:-2, 0]) / (2 * dy))**2))
    return b

def pressure_poisson_periodic(p, dx, dy):
    pn = np.empty_like(p)
    # print(b[1,1])
    for q in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] = ((pn[1:-1, 2:]+pn[1:-1, :-2])*dy**2 + (pn[2:, 1:-1]+pn[:-2,1:-1])*dx**2)/(2*(dx**2+dy**2))\
                  -(dx**2*dy**2)/(2*(dx**2+dy**2))*b[1:-1,1:-1]


        # one sided ns @ x = 2
        p[1:-1, -1] = p[1:-1, -2]

        p[1:-1, 0] = 1

        # object at xmax/2
        p[:int(ny / 3), int(nx / 2)-1] = p[:int(ny / 3), int(nx / 2)]
        p[:int(ny / 3), int(nx / 2) +1] = p[:int(ny / 3), int(nx / 2)]
        p[int(ny / 3)+1, int(nx / 2)] = p[int(ny / 3), int(nx / 2)]


        # # Periodic BC Pressure @ x = 2
        # p[1:-1, -1] = (((pn[1:-1, 0] + pn[1:-1, -2])* dy**2 +
        #                 (pn[2:, -1] + pn[0:-2, -1]) * dx**2) /
        #                (2 * (dx**2 + dy**2)) -
        #                dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, -1])
        #
        # # Periodic BC Pressure @ x = 0
        # p[1:-1, 0] = (((pn[1:-1, 1] + pn[1:-1, -1])* dy**2 +
        #                (pn[2:, 0] + pn[0:-2, 0]) * dx**2) /
        #               (2 * (dx**2 + dy**2)) -
        #               dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 0])
        # Wall boundary conditions, pressure
        p[-1, :] = p[-2, :]  # dp/dy = 0 at y = 2
        p[0, :] = p[1, :]  # dp/dy = 0 at y = 0

    return p


udiff = 1
stepcount = 0
running = 1
fix, ax = plt.subplots()

while running:
    plt.clf()
    un = u.copy()
    vn = v.copy()
    b = build_up_b(rho, dx, dy, u, v)
    p = pressure_poisson_periodic(p, dx, dy)

    test = (vn[1:-1,1:-1]-vn[1:-1,0:-2])/(dx)

    u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                     un[1:-1, 1:-1] * dt / dx *
                    (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy *
                    (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                     dt / (2 * rho * dx) *
                    (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                     nu * (dt / dx**2 *
                    (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                     dt / dy**2 *
                    (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])) + F * dt)

    v[1:-1, 1:-1] = vn[1:-1,1:-1] - un[1:-1,1:-1]*dt*test - vn[1:-1,1:-1]*dt*(vn[1:-1,1:-1]-vn[0:-2,1:-1])/(dy) \
                    -dt*(p[2:,1:-1]-p[0:-2,1:-1])/(rho*2*dy)\
                    +nu*dt*((vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,:-2])/(dx**2)+(vn[0:-2,1:-1]-2*vn[1:-1,1:-1]+vn[2:,1:-1])/(dy**2))


    # open ended tunnel @ x = 2
    u[:, -1] = 1
    u[:, 0] = 1
    v[:, -1] = 0
    v[:, 0] = 0

    #object at x = 1
    u[:int(ny / 3), int(nx / 2)-1] = 0
    v[:int(ny / 3), int(nx / 2)-1] = 0

    u[:int(ny / 3), int(nx / 2)+1] = 0
    v[:int(ny / 3), int(nx / 2)+1] = 0

    u[int(ny / 3)+1, int(nx / 2)] = 0
    v[int(ny / 3)+1, int(nx / 2)] = 0


    # # Periodic BC u @ x = 2
    # u[1:-1, -1] = (un[1:-1, -1] - un[1:-1, -1] * dt / dx *
    #                (un[1:-1, -1] - un[1:-1, -2]) -
    #                vn[1:-1, -1] * dt / dy *
    #                (un[1:-1, -1] - un[0:-2, -1]) -
    #                dt / (2 * rho * dx) *
    #                (p[1:-1, 0] - p[1:-1, -2]) +
    #                nu * (dt / dx ** 2 *
    #                      (un[1:-1, 0] - 2 * un[1:-1, -1] + un[1:-1, -2]) +
    #                      dt / dy ** 2 *
    #                      (un[2:, -1] - 2 * un[1:-1, -1] + un[0:-2, -1])) + F * dt)
    #
    # # Periodic BC u @ x = 0
    # u[1:-1, 0] = (un[1:-1, 0] - un[1:-1, 0] * dt / dx *
    #               (un[1:-1, 0] - un[1:-1, -1]) -
    #               vn[1:-1, 0] * dt / dy *
    #               (un[1:-1, 0] - un[0:-2, 0]) -
    #               dt / (2 * rho * dx) *
    #               (p[1:-1, 1] - p[1:-1, -1]) +
    #               nu * (dt / dx ** 2 *
    #                     (un[1:-1, 1] - 2 * un[1:-1, 0] + un[1:-1, -1]) +
    #                     dt / dy ** 2 *
    #                     (un[2:, 0] - 2 * un[1:-1, 0] + un[0:-2, 0])) + F * dt)
    #
    # # Periodic BC v @ x = 2
    # v[1:-1, -1] = (vn[1:-1, -1] - un[1:-1, -1] * dt / dx *
    #                (vn[1:-1, -1] - vn[1:-1, -2]) -
    #                vn[1:-1, -1] * dt / dy *
    #                (vn[1:-1, -1] - vn[0:-2, -1]) -
    #                dt / (2 * rho * dy) *
    #                (p[2:, -1] - p[0:-2, -1]) +
    #                nu * (dt / dx ** 2 *
    #                      (vn[1:-1, 0] - 2 * vn[1:-1, -1] + vn[1:-1, -2]) +
    #                      dt / dy ** 2 *
    #                      (vn[2:, -1] - 2 * vn[1:-1, -1] + vn[0:-2, -1])))
    #
    # # Periodic BC v @ x = 0
    # v[1:-1, 0] = (vn[1:-1, 0] - un[1:-1, 0] * dt / dx *
    #               (vn[1:-1, 0] - vn[1:-1, -1]) -
    #               vn[1:-1, 0] * dt / dy *
    #               (vn[1:-1, 0] - vn[0:-2, 0]) -
    #               dt / (2 * rho * dy) *
    #               (p[2:, 0] - p[0:-2, 0]) +
    #               nu * (dt / dx ** 2 *
    #                     (vn[1:-1, 1] - 2 * vn[1:-1, 0] + vn[1:-1, -1]) +
    #                     dt / dy ** 2 *
    #                     (vn[2:, 0] - 2 * vn[1:-1, 0] + vn[0:-2, 0])))

    # Wall BC: u,v = 0 @ y = 0,2

    u[0, :] = 0
    u[-1, :] = 0
    v[0, :] = 0
    v[-1, :] = 0

    udiff = (np.sum(u) - np.sum(un)) / np.sum(u)

    plt.quiver(X[::6, ::6], Y[::6, ::6], u[::6, ::6], v[::6, ::6])
    plt.pause(0.001)
    stepcount += 1


plt.show(block=True)

tic = timeit.default_timer()

print(tic-toc, '---', stepcount)
