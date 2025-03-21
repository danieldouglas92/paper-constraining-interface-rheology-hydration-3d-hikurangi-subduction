import Tian_water_content
import numpy as np
import matplotlib.pyplot as plt

T_range_vals = np.linspace(300, 1000, 500) + 273 # K
P_range_vals = np.linspace(0.5, 6, 500) # GPa

max_water_content = 5.5 # wt%

T_mesh, P_mesh = np.meshgrid(T_range_vals, P_range_vals)
H2O_mesh = T_mesh * 0
for i in range(len(H2O_mesh)):
    for j in range(len(H2O_mesh[i])):
        water_at_point = Tian_water_content.water_concentration_func(LR_poly_gabbro, c_sat_poly_gabbro, Td_poly_gabbro, P_mesh[i][j], T_mesh[i][j], lithology="gabbro")
        # The polynomials can allow for massive amounts of bound water (close to 50%) which is not physical. Restrict the water content to some
        # reasonable upper bound.
        if water_at_point > max_water_content:
            H2O_mesh[i][j] = max_water_content
        else:
            H2O_mesh[i][j] = water_at_point
        
plt.figure(dpi=100)
plt.title('Gabbro PT Diagram with Tian Approximation')
plt.contourf(T_mesh - 273, P_mesh, H2O_mesh, vmin=0, levels=16)
plt.colorbar(label='wt% water - Gabbro')
plt.xlabel('Temperature - Celsius')
plt.ylabel('Pressure - GPa')
plt.show()