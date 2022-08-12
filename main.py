import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.animation as animation


def rossler_system(x, y, z, dx, dy, dz, b, a=0.2, c=5.7):
    x_dot = - y - z
    y_dot = x + a * y
    z_dot = b + z * (x - c)

    dx_dot = - dy - dz
    dy_dot = dx + a * dy
    dz_dot = dz * x + z * dx - c * dz

    return x_dot, y_dot, z_dot, dx_dot, dy_dot, dz_dot


db = 0.001  # parameter step size
b = np.arange(0.1, 2, db)  # parameter range
b_choose = 0.5000000000000003
dt = 0.01  # time step
t = np.arange(1, 100, dt)  # time range

# initialize solution arrays
xs = np.empty(len(t) + 1)
ys = np.empty(len(t) + 1)
zs = np.empty(len(t) + 1)

dxs = np.empty(len(t) + 1)
dys = np.empty(len(t) + 1)
dzs = np.empty(len(t) + 1)

# initial values x0,y0,z0 for the system
xs[0], ys[0], zs[0] = (0.5, 0.3, -1.0)
dxs[0], dys[0], dzs[0] = (1.0, 0, 0)

# Save the plot points coordinates and plot the with a single call to plt.plot
# instead of plotting them one at a time, as it's much more efficient
b_maxes = []
z_maxes = []
b_mins = []
z_mins = []

x_3d = np.empty(len(t) + 1)
y_3d = np.empty(len(t) + 1)
z_3d = np.empty(len(t) + 1)

lamb_lyap = []
b_lyap = []

for B in b:
    # Print something to show everything is running
    # print(f"{R=:.2f}")
    K = 0
    lamb = 0.0

    for i in range(len(t)):
        # approximate numerical solutions to system
        x_dot, y_dot, z_dot, dx_dot, dy_dot, dz_dot = rossler_system(xs[i],
                                                                     ys[i], zs[i], dxs[i], dys[i], dzs[i], B)
        xs[i + 1] = xs[i] + (x_dot * dt)
        ys[i + 1] = ys[i] + (y_dot * dt)
        zs[i + 1] = zs[i] + (z_dot * dt)
        dxs[i + 1] = dxs[i] + (dx_dot * dt)
        dys[i + 1] = dys[i] + (dy_dot * dt)
        dzs[i + 1] = dzs[i] + (dz_dot * dt)

        K += 1
        dist = np.sqrt(dxs[i] * dxs[i] + dys[i] * dys[i] + dzs[i] * dzs[i])
        lamb += np.log(dist / 1.0)

    lamb /= K

    lamb_lyap.append(lamb)
    b_lyap.append(B)

    if (B == b_choose):
        print(B)
        x_3d = np.copy(xs)
        y_3d = np.copy(ys)
        z_3d = np.copy(zs)

    # calculate and save the peak values of the z solution
    for i in range(1, len(zs) - 1):
        # save the local maxima
        if zs[i - 1] < zs[i] and zs[i] > zs[i + 1]:
            b_maxes.append(B)
            z_maxes.append(zs[i])
        elif zs[i - 1] > zs[i] and zs[i] < zs[i + 1]:
            b_mins.append(B)
            z_mins.append(zs[i])

    # "use final values from one run as initial conditions for the next to stay near the attractor"
    xs[0], ys[0], zs[0] = xs[i], ys[i], zs[i]
    dxs[0], dys[0], dzs[0] = 1, 0, 0

# -- Plots

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')

ax.plot(x_3d, y_3d, z_3d, zdir='z',color='black',linewidth=0.5)

ax.set_xlabel('$x$', fontsize='16.5', horizontalalignment='center')
ax.set_ylabel('$y$', fontsize='16.5', horizontalalignment='center')
ax.set_zlabel('$z$', fontsize='16.5', horizontalalignment='center')
plt.show()

def rotate(angle):
    ax.view_init(azim=angle)


print("Making animation")
rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
rot_animation.save('rotation.gif', dpi=80, writer='imagemagick')


def x_y(xs, ys):
    fig = plt.figure()

    a = plt.plot(xs, ys, color="black", linewidth=0.5, label=' a = 0.2\n c = 5.7\n b= 0.5')

    plt.xlabel('x', fontsize='16.5', horizontalalignment='center')
    plt.ylabel('y', fontsize='16.5', horizontalalignment='center')
    plt.legend()

    return fig


def x_z(xs, zs):
    fig = plt.figure()

    a = plt.plot(xs, zs, color="black", linewidth=0.5, label=' a = 0.2\n c = 5.7\n b= 0.5')

    plt.xlabel('x', fontsize='16.5', horizontalalignment='center')
    plt.ylabel('z', fontsize='16.5', horizontalalignment='center')
    plt.legend()

    return fig


def y_z(ys, zs):
    fig = plt.figure()

    a = plt.plot(ys, zs, color="black", linewidth=0.5, label=' a = 0.2\n c = 5.7\n b= 0.5')

    plt.xlabel('y', fontsize='16.5', horizontalalignment='center')
    plt.ylabel('z', fontsize='16.5', horizontalalignment='center')
    plt.legend()

    return fig


def z_b(z_maxes, b_maxes):
    fig = plt.figure()

    a = plt.scatter(b_maxes, z_maxes, color="black", s=1, alpha=1)

    plt.xlabel('b', fontsize='16.5', horizontalalignment='center')
    plt.ylabel('z', fontsize='16.5', horizontalalignment='center')

    return fig


def lamb_b(lamb_lyap, b_lyap, b):
    fig = plt.figure()

    a = plt.plot(b_lyap, lamb_lyap, color="black", label=
    'Expoente de Lyapunov')
    b = plt.plot([min(b), max(b)], [0, 0], color="red", label='Função: y = 0')

    plt.xlabel('b', fontsize='16.5', horizontalalignment='center')
    plt.ylabel('z', fontsize='16.5', horizontalalignment='center')
    plt.legend()

    return fig


fig1 = x_y(x_3d, y_3d)
fig2 = x_z(x_3d, z_3d)
fig3 = y_z(y_3d, z_3d)
fig4 = z_b(z_maxes, b_maxes)
fig5 = lamb_b(lamb_lyap, b_lyap, b)

pp = PdfPages('rossler.pdf')
pp.savefig(fig1)
pp.savefig(fig2)
pp.savefig(fig3)
pp.savefig(fig4)
pp.savefig(fig5)
pp.close()