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
                # todo we only look at NNN interactions now, the NNN distance is sqrt(5)
                if not np.linalg.norm(r) > np.sqrt(5):
                    E_total += inter(configuration[pos1], configuration[pos2], r)
    return E_total


def construct_4x2(rows,  cols, equil_angle=47 / 450 * np.pi):
    right_up = np.array([np.cos(equil_angle), np.sin(equil_angle), 0])
    left_up = np.array([-np.cos(equil_angle), np.sin(equil_angle), 0])

    config = {}
    for i in range(rows):
        for j in range(cols):
            pos = (2 * j, 0, i)
            # I think if row + col is even, we have one orientation and if odd, we have the other
            if (i + j)  % 2 == 0:
                config[pos] = right_up
            else:
                config[pos] = left_up

    return config

def construct_2x2(rows, cols, equil_angle=47 / 450 * np.pi):
    right_up = np.array([np.cos(equil_angle), np.sin(equil_angle), 0])
    left_up = np.array([-np.cos(equil_angle), np.sin(equil_angle), 0])

    config = {}
    for i in range(rows):
        for j in range(cols):
            pos = (2 * j, 0, i)
            # now the orientation only depends on the row

            if i % 2 == 0:
                config[pos] = right_up
            else:
                config[pos] = left_up

    return config


def main():
    equil_angle = 47 / 450 * np.pi
    interaction = dipol_int

    right_up = np.array([np.cos(equil_angle), np.sin(equil_angle), 0])
    left_up = np.array([-np.cos(equil_angle), np.sin(equil_angle), 0])
    right_down = np.array([np.sin(equil_angle), -np.cos(equil_angle), 0])
    left_down = np.array([-np.sin(equil_angle), -np.cos(equil_angle), 0])

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

    not_as_high_config = {
        (0, 0, 0): right_up, (2, 0, 0): left_up, (4, 0, 0): right_up,
        (0, 0, 1): right_up, (2, 0, 1): left_up, (4, 0, 1): right_up,
        (0, 0, 2): right_up, (2, 0, 2): left_up, (4, 0, 2): right_up,
        (0, 0, 3): right_up, (2, 0, 3): left_up, (4, 0, 3): right_up,
        (0, 0, 4): right_up, (2, 0, 4): left_up, (4, 0, 4): right_up
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

    four_low_config = {
        (0, 0, 0): left_up,                     (2, 0, 0): right_up,
        (0, 0, 1): right_up, (2, 0, 1): left_up,
                             (2, 0, 2): right_up,
    }

    four_high_config = {
        (0, 0, 0): left_up,                                (2, 0, 0): left_up,
        (0, 0, 1): right_up,    (2, 0, 1): right_up,
                                (2, 0, 2): left_up,
    }

    low_energy = calc_energy(low_config, interaction)
    high_energy = calc_energy(high_config, interaction)
    very_high_energy = calc_energy(very_high_config, interaction)
    not_as_high_energy = calc_energy(not_as_high_config, interaction)
    left_very_high_energy = calc_energy(left_very_high_config, interaction)
    four_low_config_energy = calc_energy(four_low_config, interaction)
    four_high_config_energy = calc_energy(four_high_config, interaction)

    print(low_energy)
    print(high_energy)
    print(very_high_energy)
    print(left_very_high_energy)
    print("4 - config low: ", four_low_config_energy)
    print("4 - config high: ", four_high_config_energy)
    print(not_as_high_energy)

    pos_0 = np.array((0, 0, 0))
    mom_0 = right_up

    pos_10 = np.array((2, 0, -1))
    pos_11 = np.array((2, 0, 0))
    pos_12 = np.array((2, 0, 1))
    mom_10 = right_up
    mom_11 = left_up
    mom_12 = right_up

    E_outer = dipol_int(mom_0, mom_10, pos_10 - pos_0) + dipol_int(mom_0, mom_12, pos_12 - pos_0)
    E_outer_unaligned = dipol_int(mom_0, left_up, pos_10 - pos_0) + dipol_int(mom_0, left_up, pos_12 - pos_0)
    E_middle_unaligned = dipol_int(mom_0, mom_11, pos_11 - pos_0)
    E_middle_aligned = dipol_int(mom_0, right_up, pos_11 - pos_0)


    print(E_outer)
    print(E_outer_unaligned)
    print(E_middle_unaligned)
    print(E_middle_aligned)

    print(f"config antisymmetric: E = {E_outer + E_middle_unaligned}")
    print(f"config symmetric: E = {E_outer_unaligned + E_middle_aligned}")

    large_4x2 = construct_4x2(1000, 5)
    large_2x2 = construct_2x2(1000, 5)

    E_large_4x2 = calc_energy(large_4x2, interaction)
    E_large_2x2 = calc_energy(large_2x2, interaction)

    print(f"large config antisymmetric: E = {E_large_4x2}")
    print(f"large config symmetric: E = {E_large_2x2}")

if __name__ == "__main__":
    main()