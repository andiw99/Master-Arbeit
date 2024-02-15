import numpy as np


def solve_linear_system(equation_list, energies):
    unknowns = set()
    for equation in equation_list:
        unknowns.update(equation.keys())
    unknowns = sorted(unknowns)
    print(unknowns)
    num_unknowns = len(unknowns)
    A = np.zeros((num_unknowns, num_unknowns))
    b = np.zeros(num_unknowns)

    #for i, unknown in enumerate(unknowns):
    #    for j, equation in enumerate(equation_list):
    #        coef = equation.get(unknown, 0)
    #        b[i] -= coef * equation.get(inhomogeneous_key, 0)
    #        A[i, j] = coef

    for j, equation in enumerate(equation_list):
        # so i guess j is now just the row
        for i, unknown in enumerate(unknowns):
            A[j, i] = equation[unknown]
        b[j] = energies[j]

    print("A:")
    print(A)
    print("b:")
    print(b)

    try:
        x = np.linalg.solve(A, b)
        solution = {unknowns[i]: x[i] for i in range(num_unknowns)}
        return solution
    except np.linalg.LinAlgError:
        print("Error: Singular matrix - system is not solvable.")
        return None

def main():
    # PBEsol
    theta4x2 = 1.2247
    theta2x2 = 1.2249
    theta2x1 = 1.2366
    theta4x1 = 1.2374
    theta_SDF_5L = 1.2343
    theta_SDF_3L = 1.2343
    theta_SDF_4L = 1.2587

    # PBE
    # theta4x2 = 1.236
    # theta2x2 = 1.238
    # theta2x1 = 1.243
    # theta4x1 = 1.247
    # theta_SDF_3L = 1.244
    # theta_SDF_4L = 1.260
    # theta_SDF_5L = 1.243


    # wait no, we calculate p just from the thet
    p = np.pi / theta2x1

    # the slap seems to have periodic boundary conditions so this prefactor is wrong? It is just two
    #d iagonal_prefactor = 2      # accounts for the fact that there are some diagonal interactions missing in the slap

    print("p = ", p )

    p2x2 = {
        "J_parallel": np.cos(4 * theta2x2) - np.cos(4 * theta4x2),
        "J_perp": 1 - np.cos(4 * theta4x2),
        "J_x": 2 * (np.cos(4 * theta2x2) - 1),
        "h": np.cos(p * theta2x2) - np.cos(p * theta4x2),
    }

    p2x1 = {
        "J_parallel": 1 - np.cos(4 * theta4x2),
        "J_perp": 1 - np.cos(4 * theta4x2),
        "J_x": 0,
        "h": np.cos(p * theta2x1) - np.cos(p * theta4x2),
    }

    p4x1 = {
        "J_parallel": 1 - np.cos(4 * theta4x2),
        "J_perp": np.cos(4 * theta4x1) - np.cos(4 * theta4x2),
        "J_x": 2 * (np.cos(4 * theta4x1) - 1),
        "h": np.cos(p * theta4x1) - np.cos(p * theta4x2),
    }

    sdf = {
        "J_parallel": np.cos(2 * (theta_SDF_5L - theta_SDF_4L)) + np.cos(2 * (theta_SDF_3L - theta_SDF_4L))  +
                      np.cos(2 * (theta_SDF_3L - theta4x2)) + np.cos(2 * (theta_SDF_5L - theta4x2)) - 4 * np.cos(4 * theta4x2),
        "J_perp": 2 * (np.cos(2 * (theta4x2 - theta_SDF_4L)) + np.cos(2 * (theta4x2 - theta_SDF_5L)) + np.cos(2 * (theta4x2 - theta_SDF_3L))) - 6 * np.cos(4 * theta4x2),
        "J_x": 2 * (np.cos(2 * (theta4x2 - theta_SDF_4L)) + np.cos(2 * (theta4x2 - theta_SDF_5L)) + np.cos(2 * (theta4x2 - theta_SDF_3L)) - 3),
        "h": np.cos(p * theta_SDF_4L) + np.cos(p * theta_SDF_5L) + np.cos(p * theta_SDF_3L) - 3 * np.cos(p * theta4x2),
    }

    c4x2_equilibrium_angle = {
        "J_parallel": 4 * np.sin(4 * theta4x2),
        "J_perp": 4 * np.sin(4 * theta4x2),
        "J_x": 0,
        "h": p * np.sin(p * theta4x2),
    }

    equations = [p2x2, p2x1, p4x1, c4x2_equilibrium_angle]
    # PBEsol
    energies = [1.36, 74.21, 101.68, 0]           # Okay you devided two by twelve, one by 24 and one by two?
    # PBE
    # energies = [1.03, 78.5, 107.8, 131.4]
    config_names = ["p(2x2)", "p(2x1)", "p(4x1)", "angle"]
    for E, (equation, name) in enumerate(zip(equations, config_names)):
        eq_str = f"{name}:\t"
        for i, parameter in enumerate(equation):
            if i != 0:
                eq_str += " + "
            eq_str += f"{equation[parameter]:.7f} {parameter}"
        eq_str += f" = {energies[E]:.2f} meV"
        print(eq_str)

    equation_list = [p2x2, p2x1, p4x1, sdf]

    print(solve_linear_system(equation_list, energies))


if __name__ == "__main__":
    main()