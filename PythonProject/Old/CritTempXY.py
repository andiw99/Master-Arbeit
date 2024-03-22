from FunctionsAndClasses import *

J_para = 10
J_perp = 0.1

print(T_c_est(J_para, J_perp), T_c_est(J_para, J_perp)[0] / 1300)

print(T_c_est_Ising_eff_ratio(31))
print(T_c_est_Ising_eff(T_c_est_Ising_eff_ratio(31)))

print(derivative_F_BJ_Ising(T_c_est_Ising_eff_ratio(31)))

print(calculate_angle_from_slopes(-10.25, 31))
print(calculate_angle_from_slopes(-12.34, 31))
