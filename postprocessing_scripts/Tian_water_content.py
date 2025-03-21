import numpy as np
import matplotlib.pyplot as plt

def water_concentration_func(pressure, temperature, cutoff_pressure, lithology):
    """ 
    This function takes polynomials used to define the phase diagram for hydration of a bulk rock composition from Tian et al. 2019, as well as
    a temperature and pressure, and outputs the partition coefficient of the rock at the specified P-T conditions.
    
    lithology: A string that specifies which rock composition from Tian et al., 2019 is being used. Can either be: "sediment", "MORB", "gabbro", or "peridotite".
    cutoff_pressure: The pressure, in GPa, that defines the maximum allowed pressure to be used in the Tian parameterization. This is necessary because the polynomials break down at high pressures leading to infinite maximum bound water contents.
    pressure: The pressure, in GPa, of the given lithology
    temperature: The temperature, in K, of the given lithology
    """
    
    # The polynomials defined in Tian et al., 2019 break down above certain pressures. If the input pressure is above this pressure,
    # set it to the highest pressure where the polynomials do NOT breakdown
    if lithology == "sediment":
        H2O_LR_poly = np.array([-2.03283, 10.8186, -21.2119, 18.3351, -6.48711, 8.32459])
        H2O_c_sat_poly = np.array([-0.150662, 0.301807, 1.01867])
        H2O_Td_poly = np.array([2.83277, -24.7593, 85.9090, 524.898])
        
    if lithology == "gabbro":
        H2O_LR_poly = np.array([-1.81745, 7.67198, -10.8507, 5.09329, 8.14519])
        H2O_c_sat_poly = np.array([-0.0176673, 0.0893044, 1.52732])
        H2O_Td_poly = np.array([-1.72277, 20.5898, 637.517])
        if pressure > cutoff_pressure:
            pressure = cutoff_pressure
            
    if lithology == "MORB":
        H2O_LR_poly = np.array([-1.78177, 7.50871, -10.4840, 5.19725, 7.96365])
        H2O_c_sat_poly = np.array([0.0102725, -0.115390, 0.324452, 1.41588])
        H2O_Td_poly = np.array([-3.81280, 22.7809, 638.049])
        if pressure > cutoff_pressure:
            pressure = cutoff_pressure
            
    if lithology == "peridotite":
        H2O_LR_poly = np.array([-19.0609, 168.983, -630.032, 1281.84, -1543.14, 1111.88, -459.142, 95.4143, 1.97246])
        H2O_c_sat_poly = np.array([0.00115628, 2.42179])
        H2O_Td_poly = np.array([-15.4627, 94.9716, 636.603])
        if pressure > cutoff_pressure:
            pressure = cutoff_pressure

    # Initialize the variables
    inv_pressure = 1/pressure
    ln_LR_val = 0
    ln_c_sat_val = 0
    Td_val = 0

    # Calculate the values for LR, c_sat, and Td
    for i in range(len(H2O_LR_poly)):
        ln_LR_val += H2O_LR_poly[i] * (inv_pressure**(len(H2O_LR_poly) - 1 - i))

    for j in range(len(H2O_c_sat_poly)):
        if lithology == "sediment":
            ln_c_sat_val += H2O_c_sat_poly[j] * (np.log10(pressure)**(len(H2O_c_sat_poly) - 1 - j))
            
        else:
            ln_c_sat_val += H2O_c_sat_poly[j] * (pressure**(len(H2O_c_sat_poly) - 1 - j))

    for k in range(len(H2O_Td_poly)):
        Td_val += H2O_Td_poly[k] * (pressure**(len(H2O_Td_poly) - 1 - k))

    LR_val = np.exp(ln_LR_val)
    c_sat_val = np.exp(ln_c_sat_val)

    # Calculate the partition coefficient (amount of water that can be stored in the rock)
    partition_coeff = c_sat_val * np.exp(LR_val * (1/temperature - 1/Td_val))
    return partition_coeff


################## Create the PT diagram for a given lithology ##################

T_range_vals = np.linspace(300, 1000, 500) + 273 # K
P_range_vals = np.linspace(0.5, 10, 500) # GPa

max_water_content = 10.5 # wt%
cutoff_pressure = 26
lithology_string = "peridotite"

T_mesh, P_mesh = np.meshgrid(T_range_vals, P_range_vals)
H2O_mesh = T_mesh * 0
for i in range(len(H2O_mesh)):
    for j in range(len(H2O_mesh[i])):
        water_at_point = water_concentration_func(P_mesh[i][j], T_mesh[i][j], cutoff_pressure, lithology=lithology_string)
        # The polynomials can allow for massive amounts of bound water (close to 50%) which is not physical. Restrict the water content to some
        # reasonable upper bound.
        if water_at_point > max_water_content:
            H2O_mesh[i][j] = max_water_content
        else:
            H2O_mesh[i][j] = water_at_point
        
plt.figure(dpi=100)
plt.title('Peridotite', fontsize=25)
plt.contourf(T_mesh - 273, P_mesh, H2O_mesh, vmin=0, levels=32)
cb = plt.colorbar()
cb.set_label(label='wt% Water', size=20)
cb.ax.tick_params(labelsize=15) 


# cb.ax.set_title('Color Scale', fontsize=14)
plt.xlabel('Temperature - Celsius', fontsize=20)
plt.ylabel('Pressure - GPa', fontsize=20)
plt.tick_params(labelsize=15)
plt.savefig(lithology_string + "_phase_diagram.png", bbox_inches='tight')
# plt.show()
