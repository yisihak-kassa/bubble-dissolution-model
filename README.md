# Single Bubble Dissolution Model

This repository contains a model and a data analysis workflow I developed for my master's thesis research. The work was carried out as part of an experiment conducted at the Helmholtz Centre for Environmental Research in Germany which focused on quantifying the transport of gases by rising bubbles in a water column. The code for simulating gas exchange is contained in `bubble_tracker.py`.

## Important Notes

- **Data**: Due to ownership and privacy restrictions, the actual experimental data used in the thesis is not included in this repository. However, all scripts for data manipulation, analysis, and plotting are available.

The model simulates gas exchange for a moving bubble by integrating the following mass-transfer differential equation using Euler's method. For a rising bubble containing a mix of gases, the mass transfer of a gas $i$ is described by

```math
\frac{dM_i}{dz} = -K_{Li} \left( H_i P_i - C_i \right) \frac{4\pi r^2}{v_b}
```

where  
- $M_i$ is the mass of gas $i$ in the bubble (mol)  
- $K_{Li}$ is the mass transfer coefficient (m/s)  
- $H_i$ is Henry’s law constant (mol/L bar)  
- $P_i$ is the partial pressure of gas $i$ (bar)  
- $C_i$ is the concentration of gas $i$ in the liquid phase (mol/L)  
- $r$ is the bubble radius (m)  
- $v_b$ is the bubble rise velocity (m/s).  

The term $4\pi r^2$ represents the bubble surface area. Since gas exchange occurs at the interface between the gas and liquid phases, this surface area directly influences the rate of mass transfer. The term $H_i P_i - C_i$ represents the concentration gradient that drives gas transfer. If the dissolved gas profile remains constant throughout the depth of the water column, then $C_i$ remains unchanged. In this case, the concentration gradient is determined solely by the partial pressure of gas $i$ inside the bubble.  

At greater depths, the total pressure in the bubble is higher, which increases the partial pressure of the gas compared to lower depths. Consequently, $H_i P_i$ becomes larger, resulting in a greater concentration gradient at these depths. However, as the bubbles ascend, the concentration gradient gradually decreases and approaches zero. This leads to a reduction in the rate of gas transfer with increasing travel distance.

