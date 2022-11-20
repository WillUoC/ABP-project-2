import numpy as np
import matplotlib.pyplot as plt

# Initial values
P_t1 = 101325                   # Total ambient pressure (Pa) = N/m^2 
T_t1 = 288.15                   # Total ambient temperature (K)
a_a = 340.3                     # Ambient speed of sound (m/s)
rho_t1 = 1.225                  # Ambient density (kg/m^3) 
R = 287                         # Universal gas constant (J/(kg*K))
T_t4 = 2150                     # Turbine inlet temperature (K)
n_pt = 0.93                     # Turbine polytropic efficiency
pi_f = 1.6                      # Fan pressure ratio
n_pf = 0.91                     # Fan polytropic efficency
n_pc = 0.9                      # Compressor polytropic efficiency
pi_c = 18                       # Compressor pressure ratio
BPR = 7                         # Bypass ratio
md_a = 520                      # Total mass flow (kg/s)
md_h = md_a/(BPR+1)             # Hot/core mass flow (kg/s)
md_c = (md_a*BPR)/(BPR+1)       # Cold/bypass mass flow (kg/s)
gamma_c = 1.4                   # Cold gamma
gamma_h = 1.333                 # Hot gamma
cpa = 1.005 * 1000 #J/(kg*K)    # Cold ratio of specific heats
cpg = 1.148 * 1000 #J/(kg*K)    # Hot ratio of specific heats



# RESOLUTION = 1000



# Iterate
#Varying Variables


# n_pc_nominal = 0.92 # 0.9 - 0.99
# Ca_nominal = 150 #m/s 150 - 200
# t_h_rat_nominal = 0.3 # 0 - 1
# M_nominal = 1.1 # 0 - 1.2

# n_pc_a = np.linspace(0.9, 0.99, RESOLUTION)
# Ca_a = np.linspace(150, 200, RESOLUTION)
# t_h_rat_a = np.linspace(0, 1, RESOLUTION)
# M_a = np.linspace(0, 1.2, RESOLUTION)

# iter_array = np.dstack((
#     np.vstack((n_pc_a, Ca_nominal * np.ones(RESOLUTION), t_h_rat_nominal * np.ones(RESOLUTION), M_nominal * np.ones(RESOLUTION))).T,
#     np.vstack((n_pc_nominal * np.ones(RESOLUTION), Ca_a, t_h_rat_nominal * np.ones(RESOLUTION), M_nominal * np.ones(RESOLUTION))).T,
#     np.vstack((n_pc_nominal * np.ones(RESOLUTION), Ca_nominal * np.ones(RESOLUTION), t_h_rat_a, M_nominal * np.ones(RESOLUTION))).T,
#     np.vstack((n_pc_nominal * np.ones(RESOLUTION), Ca_nominal * np.ones(RESOLUTION), t_h_rat_nominal * np.ones(RESOLUTION), M_a)).T
# ))

# Varying Variables
C_a = 150 # m/s 150 - 200   Axial Velocity
t_h_rat = 0.3 # 0 - 1       Tip-Hub Ratio
M = 1.1 # 0 - 1.2           Tip Mach Number
lam = 1     #               Degree of reaction
deHaller = 0.764    #       deHaller Criteria

# Thermo analysis
T_1 = T_t1 - ((C_a**2)/(2*cpa))
P_1 = P_t1*((T_1/T_t1)**(gamma_c/(gamma_c-1)))
rho_1 = P_1/(R*T_1)
a_1 = np.sqrt(gamma_c*R*T_1)

# Fan
P_t2 = P_t1*(pi_f)

T_t2 = T_t1 + T_t1*((pi_f**((gamma_c-1)/(gamma_c*n_pf)))-1)

T_2 = T_t2 - ((C_a**2)/(2*cpa))
P_2 = P_t2*((T_2/T_t2)**(gamma_c/(gamma_c-1)))
rho_2 = P_2/(R*T_2)
a_2 = np.sqrt(gamma_c*R*T_2)

w_f = cpa*T_t1*((pi_f**((gamma_c-1)/(gamma_c*n_pf)))-1) # Work for driving the fan

# Compressor 
P_t3 = (pi_c)*P_t2

T_t3 = T_t2*((pi_c**((gamma_c-1)/(gamma_c*n_pc)))-1) + T_t2

T_3 = T_t3 - ((C_a**2)/(2*cpa))
P_3 = P_t3*((T_3/T_t3)**(gamma_c/(gamma_c-1)))
rho_3 = P_3/(R*T_3)

A_ec = md_h/(rho_3*C_a)                                 # Compressor exit area

w_c = cpa*(T_t3 - T_t2)                                 # Work for diving the compressor

#### Design Process ####

# Dependent Variables
r_ft = np.sqrt(md_a/(np.pi*rho_1*C_a*(1-t_h_rat**2)))   # fan tip radius
A_t = np.pi*r_ft**2                                     # total area
A_c = A_t/(BPR+1)                                       # compressor area
r_fh = r_ft * t_h_rat                                   # fan hub radius
r_ch = r_fh                                             # compressor hub radius equals fan hub radius
r_ct = np.sqrt(A_c/(2*np.pi) + r_ch**2)                 # compressor tip radius
ht_rat_c = r_ch/r_ct                                    # compressor hub-to-tip ratio
r_m = (r_ct + r_ch)/2                                   # compressor blade mean radius


## Compressor
print("### Compressor ###")

U_ct = M * a_2                                          # compressor tip speed
N = U_ct/(2*np.pi*r_ct)                                 # revolutions per second

# Compressor Mean Line
U_m = 2*np.pi*N*r_m                                     # blade mean radius tangential velocity
beta_1_m = np.arctan(U_m/C_a)                           # angle of relative velocity
V_2 = C_a/np.cos(beta_1_m)                              # relative velocity
V_3 = deHaller * V_2                                    # velocity with maximum turining
beta_2_m = np.arccos(C_a/V_3)                           # angle after blade with max turning

# Compute stages based on constant mean line
delt_T_stage = lam*U_m*C_a*(np.tan(beta_1_m)-np.tan(beta_2_m))/cpa  # stagnation temperature change over stage
num_stages = (T_t3 - T_t2)/delt_T_stage                             # calculate stages by total temp change divided by stage temp change

print("\nTemperature change per stage:",round(delt_T_stage, 2),"K")
print("Total change in temp across compressor:", round(T_t3 - T_t2, 2))
print("Number of stages:", int(num_stages+1))

# Compressor Tip
U_t = U_ct
beta_1_t = np.arctan(U_t/C_a)
V_2_t = C_a/np.cos(beta_1_t)
V_3_t = deHaller * V_2_t
beta_2_t = np.arccos(C_a/V_3)

# Compressor Root
U_r = 2*np.pi*N*r_ch
beta_1_r = np.arctan(U_r/C_a)
V_2_r = C_a/np.cos(beta_1_r)
V_3_r = deHaller * V_2_r
beta_2_r = np.arccos(C_a/V_3_r)

print("\n\nCompressor Blade Angles:")

print("\n  Tip Line")
print("   Beta 1:", round(beta_1_t*180/np.pi, 2), "degrees")
print("   Beta 2:", round(beta_2_t*180/np.pi, 2), "degrees")

print("\n  Mean Line")
print("   Beta 1:", round(beta_1_m*180/np.pi, 2), "degrees")
print("   Beta 2:", round(beta_2_m*180/np.pi, 2), "degrees")

print("\n  Root Line")
print("   Beta 1:", round(beta_1_r*180/np.pi, 2), "degrees")
print("   Beta 2:", round(beta_2_r*180/np.pi, 2), "degrees")

h = A_ec/(2*np.pi*r_m) # Blade height at exit

r_te = r_m + h/2 # Tip radius at exit
r_re = r_m - h/2 # Root radius at exit



## Turbine
print("\n\n### Turbine ###")

# Varying Variables
psi = 2.8                                       # Blade loading coefficient, cannot exceed 3
phi = 0.77                                      # Flow coefficient, cannot be less than 0.75
C_a = 300                                       # Axial flow in turbine (m/s)


delta_T0s = w_c/cpg                             # Temperature change over turbine section
print("\nTurbine Temperature Drop:", round(delta_T0s, 2), "K")


# Degree of reaction
lam_m = 0.5                                     # 50% degree of reaction on the mean line

########################

# Need to figure out pressure/temperature/density/area changes for more accurate math
# Max diameter of turbine outlet cannot be larger than first compressor stage diameter
# Need to specify tip/hub ratio for turbine blades
# Find how area changes with density to conserve mass flow and axial velocity

########################

# Turbine Outlet Area
T_e = T_t4 - delta_T0s - C_a**2/(2*cpg)                         # Static temperature at turbine exit (K)
P_t4 = P_t3                                                     # No pressure loss in combustor (Pa)

inv_P_rat = (1 - delta_T0s/T_t4)**(gamma_h/(n_pt*(gamma_h-1)))  # Pressure ratio P_te/P_t4
P_te = P_t4 * inv_P_rat                                         # Stagnation exit pressuer (Pa)
P_e = P_te - C_a**2/(2*cpg)                                     # Static exit pressure (Pa)

rho_e = P_e/(R * T_e)                                           # Static exit density (kg/m^3)

A_et = md_h/(rho_e * C_a)                                       # Turbine exit area (m^2)
print("Turbine Exit Area:       ",round(A_et,6),"m^2")

#   Make exit turbine stage blade tip equal to compressor tip
r_turb_t = r_ct                                                 # Largest turbine tip radius (m)
r_turb_r = np.sqrt(r_turb_t**2 - A_et/np.pi)                    # Turbine exit root radius (m)
r_turb_m = (r_turb_t + r_turb_r)/2                              # Turbine exit mean radius (m)
print("Turbine Exit Tip Radius: ", round(r_turb_t, 4),"m")
print("Turbine Exit Mean Radius:", round(r_turb_m, 4),"m")
print("Turbine Exit Root Radius:", round(r_turb_r, 4),"m")

#   Use mean radius to find blade speed to calculate temperature change per stage
U_tm = 2*np.pi*N*r_turb_m                                       # Turbine mean line velocity (m/s)
D_T0s = psi*U_tm**2/(2*cpg)                                     # Change in stagnation temperature per stage

#   Calculate number of stages by dividing the total temperature drop by the drop per stage
t_stages = delta_T0s/D_T0s
print("Turbine Stages:", int(t_stages+1))


# Blade angles
# Assume flow is axial coming into nozzle blades
#   Calculate phi values for tip, mean, root
phi_t = C_a/(2*np.pi*N*r_turb_t)
phi_m = C_a/(2*np.pi*N*r_turb_m)
phi_r = C_a/(2*np.pi*N*r_turb_r)

print("Tip Phi: ", round(phi_t,3))
print("Mean Phi:", round(phi_m,3))
print("Root Phi:", round(phi_r,3))

# Mean Line
beta_2_m = np.arctan(1/(2*phi_m)*(psi/2-2*lam_m))
alpha_2_m = np.arctan(np.tan(beta_2_m)+1/phi_m)
beta_3_m = np.arctan(1/(2*phi_m)*(psi/2+2*lam_m))
alpha_3_m = np.arctan(np.tan(beta_3_m)-1/phi_m)

# Tip Line
beta_2_t = np.arctan(1/(2*phi_t)*(psi/2-2*lam_m))
alpha_2_t = np.arctan(np.tan(beta_2_m)+1/phi_t)
beta_3_t = np.arctan(1/(2*phi_t)*(psi/2+2*lam_m))
alpha_3_t = np.arctan(np.tan(beta_3_m)-1/phi_t)

# Root Line
beta_2_r = np.arctan(1/(2*phi_r)*(psi/2-2*lam_m))
alpha_2_r = np.arctan(np.tan(beta_2_m)+1/phi_r)
beta_3_r = np.arctan(1/(2*phi_r)*(psi/2+2*lam_m))
alpha_3_r = np.arctan(np.tan(beta_3_m)-1/phi_r)

#Absolute Velocity 
C_2_m = C_a / np.cos(alpha_2_m)
C_2_t = C_a / np.cos(alpha_2_t)
C_2_r = C_a / np.cos(alpha_2_r)

#Temperature Calculation
T_2_m = T_t4 - (C_2_m**2)/(2*cpg)
T_2_t = T_t4 - (C_2_t**2)/(2*cpg)
T_2_r = T_t4 - (C_2_r**2)/(2*cpg)   

#Pressure Ratio Calculation
Pr_2_m = (T_t4/T_2_m)**(gamma_h/(gamma_h-1))
Pr_2_t = (T_t4/T_2_t)**(gamma_h/(gamma_h-1))     #These values are below the critical pressure ratio and so unchoked
Pr_2_r = (T_t4/T_2_r)**(gamma_h/(gamma_h-1))

#Pressure Calculation
P_2_m = P_t3/Pr_2_m
P_2_t = P_t3/Pr_2_t
P_2_r = P_t3/Pr_2_r

# Display angles

print("\n\nTurbine Blade Angles:")

print("\n  Tip Line")
print("   Beta 2:", round(beta_2_t*180/np.pi, 2), "degrees")
print("   Beta 3:", round(beta_3_t*180/np.pi, 2), "degrees")

print("\n  Mean Line")
print("   Beta 2:", round(beta_2_m*180/np.pi, 2), "degrees")
print("   Beta 3:", round(beta_3_m*180/np.pi, 2), "degrees")

print("\n  Root Line")
print("   Beta 2:", round(beta_2_r*180/np.pi, 2), "degrees")
print("   Beta 3:", round(beta_3_r*180/np.pi, 2), "degrees")



















