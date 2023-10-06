import numpy as np

def dipol_int(p1, p2, r):
    """
    calculates the ineraction energy of two dipoles
    :param p1: first dipole moment as 3-vector (numpy array)
    :param p2: second dipole moment as 3-vector (numpy array)
    :param r: distance vector between the vectors
    :return: interaction energy as float
    """
    r_abs = np.linalg.norm(r)
    E = (np.dot(p1, p2) - 3 * (np.dot(p1, r / r_abs) * np.dot(p2, r / r_abs))) / (r_abs ** 3)
    return E

def calc_energy(configuration, inter):
    """
    calulates the total energy of a configuration
    :param configuration: map with positions as tuples as keys and whatever is important for the interaction as value
    :param inter: interaction, can be replaced
    :return: total energy
    """
    E_total = 0
    for pos1 in configuration.keys():
        for pos2 in configuration.keys():
            # only add if the positions are not the same
            if pos1 != pos2:
                r = np.array(pos2) - np.array(pos1)
                # TODO currently counts double but for our considerations it does not matter
                E_total += inter(configuration[pos1], configuration[pos2], r)
    return E_total


def main():
    equil_angle = 47 / 450 * np.pi
    interaction = dipol_int

    right_up = np.array([np.cos(equil_angle), np.sin(equil_angle), 0])
    left_up = np.array([-np.cos(equil_angle), np.sin(equil_angle), 0])
    right_down = np.array([np.sin(equil_angle), -np.cos(equil_angle), 0])
    left_down = np.array(-[np.sin(equil_angle), -np.cos(equil_angle), 0])

    lowest_config = {
        (0, 0, 0): left_up, (2, 0, 0): right_up, (4, 0, 0): left_up,
        (0, 0, 1): right_up, (2, 0, 1): left_up, (4, 0, 1): right_up,
        (0, 0, 2): left_up, (2, 0, 2): right_down, (4, 0, 2): left_up,
        (0, 0, 3): right_up, (2, 0, 3): left_up, (4, 0, 3): right_up,
        (0, 0, 4): left_up, (2, 0, 4): right_up, (4, 0, 4): left_up
    }

    low_config = {
        (0, 0, 0): left_up, (2, 0, 0): right_up, (4, 0, 0): left_up,
        (0, 0, 1): right_up, (2, 0, 1): left_up, (4, 0, 1): right_up,
        (0, 0, 2): left_up, (2, 0, 2): right_up, (4, 0, 2): left_up,
        (0, 0, 3): right_up, (2, 0, 3): left_up, (4, 0, 3): right_up,
        (0, 0, 4): left_up, (2, 0, 4): right_up, (4, 0, 4): left_up
    }

    high_config = {
        (0, 0, 0): left_up, (2, 0, 0): left_up, (4, 0, 0): left_up,
        (0, 0, 1): right_up, (2, 0, 1): right_up, (4, 0, 1): right_up,
        (0, 0, 2): left_up, (2, 0, 2): left_up, (4, 0, 2): left_up,
        (0, 0, 3): right_up, (2, 0, 3): right_up, (4, 0, 3): right_up,
        (0, 0, 4): left_up, (2, 0, 4): left_up, (4, 0, 4): left_up
    }

    very_high_config = {
        (0, 0, 0): right_up, (2, 0, 0): right_up, (4, 0, 0): right_up,
        (0, 0, 1): right_up, (2, 0, 1): right_up, (4, 0, 1): right_up,
        (0, 0, 2): right_up, (2, 0, 2): right_up, (4, 0, 2): right_up,
        (0, 0, 3): right_up, (2, 0, 3): right_up, (4, 0, 3): right_up,
        (0, 0, 4): right_up, (2, 0, 4): right_up, (4, 0, 4): right_up
    }

    left_very_high_config = {
        (0, 0, 0): left_up, (2, 0, 0): left_up, (4, 0, 0): left_up,
        (0, 0, 1): left_up, (2, 0, 1): left_up, (4, 0, 1): left_up,
        (0, 0, 2): left_up, (2, 0, 2): left_up, (4, 0, 2): left_up,
        (0, 0, 3): left_up, (2, 0, 3): left_up, (4, 0, 3): left_up,
        (0, 0, 4): left_up, (2, 0, 4): left_up, (4, 0, 4): left_up
    }

    low_energy = calc_energy(low_config, interaction)
    high_energy = calc_energy(high_config, interaction)
    very_high_energy = calc_energy(very_high_config, interaction)
    left_very_high_energy = calc_energy(left_very_high_config, interaction)

    print(low_energy)
    print(high_energy)
    print(very_high_energy)
    print(left_very_high_energy)

if __name__ == "__main__":
    main()