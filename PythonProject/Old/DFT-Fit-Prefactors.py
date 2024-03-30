import numpy as np
from scipy.optimize import fsolve
import functools


def solve_linear_system(equation_list, energies, p):
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
            if unknown == "h":
                A[j, i] = equation[unknown](p)
            else:
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


def h_p_functions(hp, theta4x2, E_B, J_para, J_perp):
    h, p = hp
    return (h * (1 + np.cos(p * np.pi / 2)) - E_B,
            h * p * np.sin(theta4x2 * p) + 4 * np.sin(4 * theta4x2) * (J_para * J_perp))

def h_p_new(hp, E, J_para, J_perp, J_x, theta):
    h, p = hp
    return (h * (np.cos(np.pi / 2 * p) -np.cos(theta * p)) +
                 J_para * (2 * np.cos(2 * (theta  - np.pi / 2)) - 1 - np.cos(4 * theta)) +
            2 * J_x * (2 * np.cos(np.pi + 2 * theta) - 1 - np.cos(4 * theta)) - E,
                h * p * np.sin(theta * p) + 4 * np.sin(4 * theta) * (J_para * J_perp))

def h_final(h, p, E, J_para, J_x, theta):
    return (h * (np.cos(np.pi / 2 * p) -np.cos(theta * p)) +
                 J_para * (2 * np.cos(2 * (theta  - np.pi / 2)) - 1 - np.cos(4 * theta)) +
            2 * J_x * (2 * np.cos(np.pi + 2 * theta) - 1 - np.cos(4 * theta)) - E)


def eq_system(paras, equations, energies):
    """
    This is somehow supposed to represent this huge nonlinear system
    :param paras:
    :param equations:
    :param energies:
    :return:
    """
    J_para, J_perp, J_x, h, p = paras
    eq_system = ()
    # I just noticed that equaiton["h"] is a function of p
    for equation, energy in zip(equations, energies):
        eq = (equation["J_parallel"] * J_para + equation["J_perp"] * J_perp +
              equation["J_x"] * J_x + equation["h"](p) * h) - energy
        eq_system += (eq,)
    return eq_system

def h_p_flip(p, theta):
    return np.cos(np.pi / 2 * p) - np.cos(theta * p)

def h_p_sdf(p, theta_SDF_4L, theta_SDF_5L, theta_SDF_3L, theta4x2):
    return np.cos(p * theta_SDF_4L) + np.cos(p * theta_SDF_5L) + np.cos(p * theta_SDF_3L) - 3 * np.cos(p * theta4x2)

def h_p_p4x1(p, theta4x1, theta4x2):
    return np.cos(p * theta4x1) - np.cos(p * theta4x2)

def h_p_p2x1(p, theta2x1, theta4x2):
    return np.cos(p * theta2x1) - np.cos(p * theta4x2)

def h_p_p2x2(p, theta2x2, theta4x2):
    return np.cos(p * theta2x2) - np.cos(p * theta4x2)


def main():
    # PBEsol
    # theta4x2 = 1.2247
    # theta2x2 = 1.2249
    # theta2x1 = 1.2366
    # theta4x1 = 1.2374
    # theta_SDF_5L = 1.2343
    # theta_SDF_3L = 1.2343
    # theta_SDF_4L = 1.2587

    # PBE
    # theta4x2 = 1.236
    # theta2x2 = 1.238
    # theta2x1 = 1.243
    # theta4x1 = 1.247
    # theta_SDF_3L = 1.244
    # theta_SDF_4L = 1.260
    # theta_SDF_5L = 1.243

    # Same angles \vartheta = 71
    theta4x2 = 1.238
    theta2x2 = 1.238
    theta2x1 = 1.238
    theta4x1 = 1.238
    theta_SDF_3L = 1.244
    theta_SDF_4L = 1.260
    theta_SDF_5L = 1.243

    # Same angles \vartheta = 72
    theta4x2 = 1.2566
    theta2x2 = 1.2566
    theta2x1 = 1.2566
    theta4x1 = 1.2566
    theta_SDF_3L = 1.244
    theta_SDF_4L = 1.260
    theta_SDF_5L = 1.243

    # wait no, we calculate p just from the thet
    p = np.pi / theta4x2
    p = 2.5
    # the slap seems to have periodic boundary conditions so this prefactor is wrong? It is just two
    #d iagonal_prefactor = 2      # accounts for the fact that there are some diagonal interactions missing in the slap

    print("p = ", p )

    p2x2 = {
        "J_parallel": np.cos(4 * theta2x2) - np.cos(4 * theta4x2),
        "J_perp": 1 - np.cos(4 * theta4x2),
        "J_x": 2 * (np.cos(4 * theta2x2) - 1),
        #"h": functools.partial(h_p_p2x2, theta2x2=theta2x2, theta4x2=theta4x2),
    }

    p2x1 = {
        "J_parallel": 1 - np.cos(4 * theta4x2),
        "J_perp": 1 - np.cos(4 * theta4x2),
        "J_x": 0,
        #"h": functools.partial(h_p_p2x1, theta2x1=theta2x1, theta4x2=theta4x2),
    }

    p4x1 = {
        "J_parallel": 1 - np.cos(4 * theta4x2),
        "J_perp": np.cos(4 * theta4x1) - np.cos(4 * theta4x2),
        "J_x": 2 * (np.cos(4 * theta4x1) - 1),
        #"h": functools.partial(h_p_p4x1, theta4x1=theta4x1, theta4x2=theta4x2),
    }

    sdf = {
        "J_parallel": np.cos(2 * (theta_SDF_5L - theta_SDF_4L)) + np.cos(2 * (theta_SDF_3L - theta_SDF_4L))  +
                      np.cos(2 * (theta_SDF_3L - theta4x2)) + np.cos(2 * (theta_SDF_5L - theta4x2)) - 4 * np.cos(4 * theta4x2),
        "J_perp": 2 * (np.cos(2 * (theta4x2 - theta_SDF_4L)) + np.cos(2 * (theta4x2 - theta_SDF_5L)) + np.cos(2 * (theta4x2 - theta_SDF_3L))) - 6 * np.cos(4 * theta4x2),
        "J_x": 2 * (np.cos(2 * (theta4x2 - theta_SDF_4L)) + np.cos(2 * (theta4x2 - theta_SDF_5L)) + np.cos(2 * (theta4x2 - theta_SDF_3L)) - 3),
        "h": functools.partial(h_p_sdf, theta_SDF_4L=theta_SDF_4L, theta_SDF_5L=theta_SDF_5L, theta_SDF_3L=theta_SDF_3L,
                               theta4x2=theta4x2),
    }

    flip = {
        "J_parallel": 2 * np.cos(2 * (theta4x2 - np.pi / 2)) - 1 - np.cos(4 * theta4x2),
        "J_perp": 0,
        "J_x": 2 * (2 * np.cos(np.pi + 2 * theta4x2) - 1 - np.cos(4 * theta4x2)),
        "h": functools.partial(h_p_flip, theta=theta4x2)
    }

    c4x2_equilibrium_angle = {
        "J_parallel": 4 * np.sin(4 * theta4x2),
        "J_perp": 4 * np.sin(4 * theta4x2),
        "J_x": 0,
        "h": p * np.sin(p * theta4x2),
    }
    # angle equations
    # equations = [p2x2, p2x1, p4x1, c4x2_equilibrium_angle]
    # PBEsol
    equations = [p2x2, p2x1, p4x1, sdf, flip]
    energies = [1.36, 74.21, 101.68, 111.3, 74]           # Okay you devided two by twelve, one by 24 and one by two?
    # PBE
    # energies = [1.03, 78.5, 107.8, 131.4]
    config_names = ["p(2x2)", "p(2x1)", "p(4x1)", "sdf"]
    for E, (equation, name) in enumerate(zip(equations, config_names)):
        eq_str = f"{name}:\t"
        for i, parameter in enumerate(equation):
            if i != 0:
                eq_str += " + "
            if parameter != "h":
                eq_str += f"{equation[parameter]:.7f} {parameter}"
            else:
                eq_str += f"{equation[parameter](p):.7f} {parameter}"
        eq_str += f" = {energies[E]:.2f} meV"
        print(eq_str)

    equation_list = [p2x2, p2x1, p4x1]#, sdf]

    solution = solve_linear_system(equation_list, energies, p=p)
    print(solution)
    J_para = solution["J_parallel"]
    J_perp = solution["J_perp"] - 2 * solution['J_x']
    J_x = solution["J_x"]
    print("J_para = ", J_para)
    print("J_perp = ", J_perp)
    
    # okay we want to try to estimate h and p with the energy barrier equation aswell as equilibrium position equation
    theta2x1 = 1.243
    theta4x2 = 1.236
    E_B = 170

    x, y = fsolve(h_p_functions, (400, 2.5), args=(theta4x2, E_B, J_para, J_perp))
    print(x, y)
    E_B = 74
    x, y = fsolve(h_p_new, (400, 2.5), args=(E_B, J_para, J_perp, J_x, theta4x2))
    print(x, y)
    h = fsolve(h_final, x0=np.array([1000]), args=(p, E_B, J_para, J_x, theta4x2))
    print(f"h = {h} with fixed p = {p}")
    equations = [p2x2, p2x1, p4x1, sdf, flip]
    energies = [1.36, 74.21, 101.68, 111.3, 74]
    J_para, J_perp, J_x, h, p = fsolve(eq_system, (100, -20, -10, 400, 2.5), args=(equations, energies))
    print(J_para, J_perp, J_x, h, p)


if __name__ == "__main__":
    main()