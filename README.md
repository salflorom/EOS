# EOS
## Author: Santiago A. Flores Roman

## Description
It computes the fluid's properties according to a given EOS.
Available EOS:
- Peng-Robinson
- Soave (Soave-Redlich-Kwong)
- Ideal (Ideal gas)
- Johnson (Lennard-Jones)
- VdW (Van der Waals)

Available properties:
- phi (fugacity coefficient)
- fugacity
- mu (chemical potential)
- idealMu (ideal part of mu)
- excessMu (excess part of mu)
- phase (liquid or vapour)
- mass (fluid's mass)
- molarMass (fluid's molar mass)
- Pc (critical pressure)
- Tc (critical temperature)
- omega (acentric factor)
- P (system's pressure)
- T (system's temperature)
- Dmolar (system's molar density)
- P_sat (saturation pressure)
- bulkModulus (bulk modulus or reciprocal of isothermal compressibility)
- zFactor (compressibility factor)
- roots (solutions of EOS [fluid's densities])

Note 1: The script uses the following units: kg, kg/mol, m, J/mol, K (only for T, and Tc), Pa.
Note 2: As the EOSs implemented don't predict exactly the saturation pressure, the user can input this value when calling ThermodynamicState. See Example 1.

## Installation
1. Clone the git repo.
```bash
   cd EOS
   pip install .
```

## Instructions
1. Import the script.
   ```python
   from EOS import EOS
   ```
3. Create the fluid.
-  Example 1:
   ```python
   Tcrit = critical_temperature  # K
   Pcrit = critical_pressure  # Pa
   w = acentric_factor
   mM = molar_mass  # kg/mol
   benzene = EOS.EOS(Tc=Tcrit, Pc=Pcrit, omega=w, molarMass=mM)
   ```
-  Example 2:
    ```python
    eps = epsilon  # K
    sig = sigma  # m
    mM = molar_mass  # kg/mol
    benzene = EOS.EOS(epsilon=eps, sigma=sig, molarMass=mM)
    ```
3.  Call Thermodynamic State
    ```python
    benzene.ThermodynamicState(eos='Peng-Robinson', P=0.0136e6, T=298)
    benzene.ThermodynamicState(eos='Ideal Gas', T=300, Dmolar=...)
    benzene.ThermodynamicState(eos='Soave', P=300, Dmolar=...)
    benzene.ThermodynamicState(eos='VdW', P=1e5, T=..., P_sat=101235)
    benzene.ThermodynamicState(eos='Johnson', P=1e5, T=...)
    ```
4. Call any fluid properties.
    ```python
    print(benzene.phi)
    print(benzene.fugacity)  # Pa
    print(benzene.mu)  # J/mol
    print(benzene.idealMu)  # J/mol
    print(benzene.excessMu)  # J/mol
    print(benzene.P_sat)  # Pa
    print(benzene.roots)  # mol/m^3
    print(benzene.Dmolar)  # mol/m^3
    print(benzene.P)  # Pa
    print(benzene.T)  # K
    print(benzene.bulkModulus)  # Pa
    print(benzene.zFactor)
    print(benzene.phase)
    ```

## Build
EOS python package made by Geordy Jomon
```bash
python3 -m pip install --upgrade build
python3 -m build
pip install .
```
