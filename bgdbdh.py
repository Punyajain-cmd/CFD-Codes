# import matplotlib.pyplot as plt
# import numpy as np
# from tqdm import tqdm

# # Constants
# N_POINTS = 51
# DOMAIN_SIZE = 1.0
# N_ITERATIONS = 1000
# TIME_STEP_LENGTH = 0.001
# KINEMATIC_VISCOSITY = 0.1
# DENSITY = 1.0
# HORIZONTAL_VELOCITY_TOP = 1.0

# N_PRESSURE_POISSON_ITERATIONS = 50
# STABILITY_SAFETY_FACTOR = 0.5

# def main():
#     # Discretization
#     element_length = DOMAIN_SIZE / (N_POINTS - 1)
#     x = np.linspace(0.0, DOMAIN_SIZE, N_POINTS)
#     y = np.linspace(0.0, DOMAIN_SIZE, N_POINTS)
#     X, Y = np.meshgrid(x, y)

#     # Initialize arrays
#     u_prev = np.zeros_like(X)
#     v_prev = np.zeros_like(X)
#     p_prev = np.zeros_like(X)

#     def central_difference_x(f):
#         diff = np.zeros_like(f)
#         diff[1:-1, 1:-1] = (
#             f[1:-1, 2:] - f[1:-1, :-2]
#         ) / (2 * element_length)
#         return diff

#     def central_difference_y(f):
#         diff = np.zeros_like(f)
#         diff[1:-1, 1:-1] = (
#             f[2:, 1:-1] - f[:-2, 1:-1]
#         ) / (2 * element_length)
#         return diff

#     def laplace(f):
#         diff = np.zeros_like(f)
#         diff[1:-1, 1:-1] = (
#             f[1:-1, :-2]
#             + f[:-2, 1:-1]
#             + f[1:-1, 2:]
#             + f[2:, 1:-1]
#             - 4 * f[1:-1, 1:-1]
#         ) / (element_length**2)
#         return diff

#     # Check time step stability
#     maximum_possible_time_step_length = (
#         0.5 * element_length**2 / KINEMATIC_VISCOSITY
#     )
#     if TIME_STEP_LENGTH > STABILITY_SAFETY_FACTOR * maximum_possible_time_step_length:
#         raise RuntimeError("Stability is not guaranteed. Reduce TIME_STEP_LENGTH.")

#     # Main loop
#     for _ in tqdm(range(N_ITERATIONS)):
#         d_u_prev_d_x = central_difference_x(u_prev)
#         d_u_prev_d_y = central_difference_y(u_prev)
#         d_v_prev_d_x = central_difference_x(v_prev)
#         d_v_prev_d_y = central_difference_y(v_prev)
#         laplace_u_prev = laplace(u_prev)
#         laplace_v_prev = laplace(v_prev)

#         # Tentative velocities
#         u_tent = (
#             u_prev
#             + TIME_STEP_LENGTH * (
#                 - (u_prev * d_u_prev_d_x + v_prev * d_u_prev_d_y)
#                 + KINEMATIC_VISCOSITY * laplace_u_prev
#             )
#         )
#         v_tent = (
#             v_prev
#             + TIME_STEP_LENGTH * (
#                 - (u_prev * d_v_prev_d_x + v_prev * d_v_prev_d_y)
#                 + KINEMATIC_VISCOSITY * laplace_v_prev
#             )
#         )

#         # Apply boundary conditions
#         u_tent[0, :] = 0.0
#         u_tent[:, 0] = 0.0
#         u_tent[:, -1] = 0.0
#         u_tent[-1, :] = HORIZONTAL_VELOCITY_TOP

#         v_tent[0, :] = 0.0
#         v_tent[:, 0] = 0.0
#         v_tent[:, -1] = 0.0
#         v_tent[-1, :] = 0.0

#         d_u_tent_d_x = central_difference_x(u_tent)
#         d_v_tent_d_y = central_difference_y(v_tent)

#         # Pressure correction
#         rhs = (
#             DENSITY / TIME_STEP_LENGTH
#             * (d_u_tent_d_x + d_v_tent_d_y)
#         )
#         p_next = np.zeros_like(p_prev)
#         for _ in range(N_PRESSURE_POISSON_ITERATIONS):
#             p_next[1:-1, 1:-1] = 0.25 * (
#                 p_prev[1:-1, :-2]
#                 + p_prev[:-2, 1:-1]
#                 + p_prev[1:-1, 2:]
#                 + p_prev[2:, 1:-1]
#                 - element_length**2 * rhs[1:-1, 1:-1]
#             )
#             # Pressure boundary conditions
#             p_next[:, -1] = p_next[:, -2]
#             p_next[0, :] = p_next[1, :]
#             p_next[:, 0] = p_next[:, 1]
#             p_next[-1, :] = 0.0

#             p_prev = p_next

#         d_p_next_d_x = central_difference_x(p_next)
#         d_p_next_d_y = central_difference_y(p_next)

#         # Final velocities
#         u_next = (
#             u_tent - TIME_STEP_LENGTH / DENSITY * d_p_next_d_x
#         )
#         v_next = (
#             v_tent - TIME_STEP_LENGTH / DENSITY * d_p_next_d_y
#         )

#         # Apply boundary conditions
#         u_next[0, :] = 0.0
#         u_next[:, 0] = 0.0
#         u_next[:, -1] = 0.0
#         u_next[-1, :] = HORIZONTAL_VELOCITY_TOP

#         v_next[0, :] = 0.0
#         v_next[:, 0] = 0.0
#         v_next[:, -1] = 0.0
#         v_next[-1, :] = 0.0

#         # Update variables
#         u_prev = np.clip(u_next, -1e6, 1e6)
#         v_prev = np.clip(v_next, -1e6, 1e6)
#         p_prev = np.clip(p_next, -1e6, 1e6)

#     # Vorticity Calculation
#     vorticity = central_difference_x(v_next) - central_difference_y(u_next)

#     # Plot 1: Streamlines
#     plt.figure()
#     plt.streamplot(X, Y, u_next, v_next, color="black")
#     plt.title("Streamlines")
#     plt.xlabel("X")
#     plt.ylabel("Y")
#     plt.show()

#     # Plot 2: Vorticity Contours
#     plt.figure()
#     plt.contourf(X, Y, vorticity, cmap="coolwarm")
#     plt.colorbar(label="Vorticity")
#     plt.title("Vorticity Contours")
#     plt.xlabel("X")
#     plt.ylabel("Y")
#     plt.show()

#     # Plot 3: Pressure Distribution
#     plt.figure()
#     plt.contourf(X, Y, p_next, cmap="viridis")
#     plt.colorbar(label="Pressure")
#     plt.title("Pressure Distribution")
#     plt.xlabel("X")
#     plt.ylabel("Y")
#     plt.show()

#     # Plot 4: U-Velocity vs Y at X=0.5
#     x_index = int(0.5 / element_length)
#     plt.figure()
#     plt.plot(u_next[:, x_index], y, label="U-velocity")
#     plt.title("U-Velocity vs Y at X=0.5")
#     plt.xlabel("U-Velocity")
#     plt.ylabel("Y")
#     plt.legend()
#     plt.grid()
#     plt.show()

#     # Plot 5: V-Velocity vs X at Y=0.5
#     y_index = int(0.5 / element_length)
#     plt.figure()
#     plt.plot(x, v_next[y_index, :], label="V-velocity")
#     plt.title("V-Velocity vs X at Y=0.5")
#     plt.xlabel("X")
#     plt.ylabel("V-Velocity")
#     plt.legend()
#     plt.grid()
#     plt.show()

# if __name__ == "_main_":
#     main()

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# Constants
N_POINTS = 51
DOMAIN_SIZE = 1.0
N_ITERATIONS = 1000
TIME_STEP_LENGTH = 0.001
KINEMATIC_VISCOSITY = 0.01
DENSITY = 1.0
HORIZONTAL_VELOCITY_TOP = 1.0

N_PRESSURE_POISSON_ITERATIONS = 50
STABILITY_SAFETY_FACTOR = 0.5
# Discretization
element_length = DOMAIN_SIZE / (N_POINTS - 1)
x = np.linspace(0.0, DOMAIN_SIZE, N_POINTS)
y = np.linspace(0.0, DOMAIN_SIZE, N_POINTS)
X, Y = np.meshgrid(x, y)

    # Initialize arrays
u_prev = np.zeros_like(X)
v_prev = np.zeros_like(X)
p_prev = np.zeros_like(X)

print(np.shape(u_prev))
def central_difference_x(f):
  diff = np.zeros_like(f)
  diff[1:-1, 1:-1] = (
    f[1:-1, 2:] - f[1:-1, :-2]
    ) / (2 * element_length)
  return diff

def central_difference_y(f):
  diff = np.zeros_like(f)
  diff[1:-1, 1:-1] = (
      f[2:, 1:-1] - f[:-2, 1:-1]
  ) / (2 * element_length)
  return diff

def laplace(f):
    diff = np.zeros_like(f)
    diff[1:-1, 1:-1] = (
        f[1:-1, :-2]
        + f[:-2, 1:-1]
        + f[1:-1, 2:]
        + f[2:, 1:-1]
        - 4 * f[1:-1, 1:-1]
    ) / (element_length**2)
    return diff

    # Check time step stability
maximum_possible_time_step_length = (
        0.5 * element_length**2 / KINEMATIC_VISCOSITY
    )
if TIME_STEP_LENGTH > STABILITY_SAFETY_FACTOR * maximum_possible_time_step_length:
        raise RuntimeError("Stability is not guaranteed. Reduce TIME_STEP_LENGTH.")

    # Main loop
for _ in tqdm(range(N_ITERATIONS)):
        d_u_prev_d_x = central_difference_x(u_prev)
        d_u_prev_d_y = central_difference_y(u_prev)
        d_v_prev_d_x = central_difference_x(v_prev)
        d_v_prev_d_y = central_difference_y(v_prev)
        laplace_u_prev = laplace(u_prev)
        laplace_v_prev = laplace(v_prev)

        # Tentative velocities
        u_tent = (
            u_prev
            + TIME_STEP_LENGTH * (
                - (u_prev * d_u_prev_d_x + v_prev * d_u_prev_d_y)
                + KINEMATIC_VISCOSITY * laplace_u_prev
            )
        )
        v_tent = (
            v_prev
            + TIME_STEP_LENGTH * (
                - (u_prev * d_v_prev_d_x + v_prev * d_v_prev_d_y)
                + KINEMATIC_VISCOSITY * laplace_v_prev
            )
        )

        # Apply boundary conditions
        u_tent[0, :] = 0.0
        u_tent[:, 0] = 0.0
        u_tent[:, -1] = 0.0
        u_tent[-1, :] = HORIZONTAL_VELOCITY_TOP

        v_tent[0, :] = 0.0
        v_tent[:, 0] = 0.0
        v_tent[:, -1] = 0.0
        v_tent[-1, :] = 0.0

        d_u_tent_d_x = central_difference_x(u_tent)
        d_v_tent_d_y = central_difference_y(v_tent)

        # Pressure correction
        rhs = (
            DENSITY / TIME_STEP_LENGTH
            * (d_u_tent_d_x + d_v_tent_d_y)
        )
        p_next = np.zeros_like(p_prev)
        for _ in range(N_PRESSURE_POISSON_ITERATIONS):
            p_next[1:-1, 1:-1] = 0.25 * (
                p_prev[1:-1, :-2]
                + p_prev[:-2, 1:-1]
                + p_prev[1:-1, 2:]
                + p_prev[2:, 1:-1]
                - element_length**2 * rhs[1:-1, 1:-1]
            )
            # Pressure boundary conditions
            p_next[:, -1] = p_next[:, -2]
            p_next[0, :] = p_next[1, :]
            p_next[:, 0] = p_next[:, 1]
            p_next[-1, :] = 0.0

            p_prev = p_next

        d_p_next_d_x = central_difference_x(p_next)
        d_p_next_d_y = central_difference_y(p_next)

        # Final velocities
        u_next = (
            u_tent - TIME_STEP_LENGTH / DENSITY * d_p_next_d_x
        )
        v_next = (
            v_tent - TIME_STEP_LENGTH / DENSITY * d_p_next_d_y
        )

        # Apply boundary conditions
        u_next[0, :] = 0.0
        u_next[:, 0] = 0.0
        u_next[:, -1] = 0.0
        u_next[-1, :] = HORIZONTAL_VELOCITY_TOP

        v_next[0, :] = 0.0
        v_next[:, 0] = 0.0
        v_next[:, -1] = 0.0
        v_next[-1, :] = 0.0

        # Update variables
        u_prev = np.clip(u_next, -1e6, 1e6)
        v_prev = np.clip(v_next, -1e6, 1e6)
        p_prev = np.clip(p_next, -1e6, 1e6)

    # Vorticity Calculation
vorticity = central_difference_x(v_next) - central_difference_y(u_next)

    # Plot 1: Streamlines
plt.figure()
plt.streamplot(X, Y, u_next, v_next, color="black")
plt.title("Streamlines")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

    # Plot 2: Vorticity Contours
plt.figure()
plt.contourf(X, Y, vorticity, cmap="coolwarm")
plt.colorbar(label="Vorticity")
plt.title("Vorticity Contours")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

    # Plot 3: Pressure Distribution
plt.figure()
plt.contourf(X, Y, p_next, cmap="viridis")
plt.colorbar(label="Pressure")
plt.title("Pressure Distribution")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

    # Plot 4: U-Velocity vs Y at X=0.5
x_index = int(0.5 / element_length)
plt.figure()
plt.plot(u_next[:, x_index], y, label="U-velocity")
plt.title("U-Velocity vs Y at X=0.5")
plt.xlabel("U-Velocity")
plt.ylabel("Y")
plt.legend()
plt.grid()
plt.show()

    # Plot 5: V-Velocity vs X at Y=0.5
y_index = int(0.5 / element_length)
plt.figure()
plt.plot(x, v_next[y_index, :], label="V-velocity")
plt.title("V-Velocity vs X at Y=0.5")
plt.xlabel("X")
plt.ylabel("V-Velocity")
plt.legend()
plt.grid()
plt.show()