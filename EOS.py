# Author: Santiago A. Flores Roman
#
# Description: It computes the fluid's properties according to a given EOS.
#   Available EOS: 
#   - Peng-Robinson
#   - Soave (Soave-Redlich-Kwong)
#   - Ideal
#
#   Available properties: phi (fugacity coefficient), fugacity, mu (chemical potential), 
#       idealMu (ideal part of mu), excessMu (excess part of mu), phase (liquid or vapour),
#       mass (fluid's mass), molarMass (fluid's molar mass), Pc (critical pressure),
#       Tc (critical temperature), omega (acentric factor), P (system's pressure),
#       T (system's temperature), Dmolar (system's molar density).
#       
#   Note: The script uses the following units: kg, kg/mol, m, J/mol, K (only for T, and Tc), Pa.
#
# Instructions: 
#   1.- Import the script.
#       Example: import EOS
#   2.- Create the fluid.
#       Example: benzene = EOS.EOS(Tc,Pc,omega,molarMass)
#   3.- Call ThermodynamicState.
#       Example: benzene.ThermodynamicState(eos='Peng-Robinson',P=0.0136e6,T=298)
#                benzene.ThermodynamicState(eos='Peng-Robinson',T=300,Dmolar=...)
#   4.- Call any fluid's property.
#       Example: print(benzene.phi)
#                print(benzene.fugacity)
#                print(benzene.mu)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import scipy.constants as const
from scipy.optimize import fsolve
import numpy as np
from numpy.polynomial import Polynomial

Rg = const.R
kb = const.Boltzmann
planck = const.h
avogadro = const.Avogadro

class EOS():
    def __init__(self,Tc,Pc,omega,molarMass):
        self.molarMass = molarMass #kg/mol
        self.mass = molarMass/avogadro #kg
        self.Pc = Pc #Pa
        self.Tc = Tc #K
        self.omega = omega
        self.bLow, self.aLowC, self.alpha, self.kappa = 0, 0, 0, 0
        self.alpha, self.aUpp, self.bUpp, self.zFactor, self.phi, self.fugacity = 0, 0, 0, 0, 0, 0
        self.mu, self.excessMu, self.idealMu = 0, 0, 0 #J/mol
        self.phase, self.eos = '', ''
        self.state = {}

    def ZFactors(self):
        coeffs = self.coeffs
        zFactor = coeffs.roots()
        zFactor = zFactor[np.imag(zFactor) == 0]
        zFactor = np.sort(np.real(zFactor))
        self.zFactor = zFactor

    def Phi(self):
        eos = self.eos
        zFactor = self.zFactor
        molarMass = self.molarMass #kg/mol
        bUpp = self.bUpp
        aUpp = self.aUpp
        sqrt2 = np.sqrt(2)
        phi = -1
        if eos == 'Peng-Robinson':
            np.seterr(invalid='ignore')
            phi = np.exp((zFactor-1)-np.log(zFactor-bUpp)-aUpp/(2*sqrt2*bUpp)*np.log((zFactor+(1+sqrt2)*bUpp)/(zFactor+(1-sqrt2)*bUpp)))
        elif eos == 'Soave':
            np.seterr(invalid='ignore')
            phi = np.exp((zFactor-1)-np.log(zFactor-bUpp)-aUpp*np.log(1+bUpp/zFactor)/bUpp)
        elif eos == 'VdW': phi = -1 # Check!
        elif eos == 'Ideal': phi = zFactor
        self.phi = phi

    def Phase(self):
        zFactor = self.zFactor
        phi = self.phi
        Tc = self.Tc
        Pc = self.Pc
        phase = self.phase
        state = self.state
        T, P = state['T'], state['P'] #K, Pa
        if (len(zFactor) == 1) or (zFactor[0] <= 0): 
            phi = phi[-1]
            zFactor = zFactor[-1]
            if (T > Tc) and (P > Pc): phase = 'super_critical_fluid'
            elif (T < Tc) and (P < Pc): phase = 'Vapour'
            elif (T < Tc) and (P > Pc): phase =  'Liquid'
        else:
            # Vapour (stable), liquid (metastable)
            if phi[2] < phi[0]: 
                phi = phi[2]
                zFactor = zFactor[2]
                phase = 'Vapour (stable) - Liquid (metastable)'
            # Vapour (metastable), liquid (stable)
            elif phi[2] > phi[0]: 
                phi = phi[0]
                zFactor = zFactor[0]
                phase = 'Vapour (metastable) - Liquid (stable)'
            # Vapour (stable), liquid (stable)
            else: 
                phi = phi[2]
                zFactor = zFactor[2]
                phase = 'Vapour (stable) -Liquid (stable)'
        if not 'Dmolar' in state.keys(): 
            rhoMolar = (P/(zFactor*Rg*T)) #mol/m^3
            # rhoMolar = PropsSI('Dmolar','T',T,'P',P,'Benzene')
            state['Dmolar'] = rhoMolar
        rhoN = P/(zFactor*kb*T) #m^-3
        # rhoN = P/(kb*T) #m^-3 Check!
        state['Dn'] = rhoN
        self.phi = phi
        self.zFactor = zFactor
        self.phase = phase
        self.state = state

    def Fugacity(self): 
        P = self.state['P']
        phi = self.phi
        fugacity = phi*P #Pa
        self.fugacity = fugacity #Pa

    def ChemicalPotential(self):
        state = self.state
        fugacity = self.fugacity #Pa
        phi = self.phi
        mass = self.mass #kg
        T, P, rhoN = state['T'], state['P'], state['Dn']
        thermalWL = planck/np.sqrt(2*np.pi*mass*kb*T) #m
        # Ideal Mu uses the density calculated from an EOS.
        idealMu = Rg*T*np.log(rhoN*thermalWL**3) #J/mol
        excessMu = Rg*T*np.log(phi) #J/mol
        mu = idealMu+excessMu #J/mol
        self.idealMu = idealMu
        self.excessMu = excessMu
        self.mu = mu

    def PengRobinsonEOS(self):
        def EOSTV(T,V,P0,ac,kappa,b): 
            Tc = self.Tc
            alpha = (1+kappa*(1-np.sqrt(T/Tc)))**2
            a = ac*alpha
            P = Rg*T/(V-b)-a/(V*(V+b)+b*(V-b))-P0
            return P
        state = self.state
        Tc = self.Tc
        Pc = self.Pc
        omega = self.omega
        bLow = 0.0777960739*Rg*Tc/Pc
        aLowC = 0.4572355289*(Rg*Tc)**2/Pc
        kappa = 0.37464+1.54226*omega-0.26992*omega**2
        T, rho, P = 0, 0, 0
        if not 'P' in state.keys():
            T = state['T']
            rho = state['Dmolar']
            V = 1/rho
            P = EOSTV(T,V,0,aLowC,kappa,bLow)
            state['P'] = P
        elif not 'T' in state.keys():
            P = state['P']
            rho = state['Dmolar']
            V = 1/rho
            T = fsolve(EOSTV,x0=1,args=(V,P,aLowC,kappa,bLow))
            state['T'] = T
        else: T, P = state['T'], state['P']
        alpha = (1+kappa*(1-np.sqrt(T/Tc)))**2
        aUpp = alpha*aLowC*P/(Rg*T)**2
        bUpp = bLow*P/(Rg*T)
        coeffs = Polynomial([-(aUpp*bUpp-bUpp**2-bUpp**3), (aUpp - 2*bUpp - 3*bUpp**2), bUpp-1, 1])
        self.bLow = bLow
        self.aLowC = aLowC
        self.alpha = alpha
        self.kappa = kappa
        self.aUpp = aUpp
        self.bUpp = bUpp
        self.coeffs = coeffs
        self.state = state

    def VdWEOS(self):
        def EOSTV(T,V,P0,a,b): 
            P = Rg*T/(V-b)-a/V**2-P0
            return P
        Tc = self.Tc
        Pc = self.Pc
        aLow = 27*(Rg*Tc)**2/(64*Pc)
        bLow = Rg*Tc/(8*Pc)
        T, rho, P = 0, 0, 0
        if not 'P' in state.keys():
            T = state['T']
            rho = state['Dmolar']
            V = 1/rho
            P = EOSTV(T,V,0,aLow,bLow)
            state['P'] = P
        elif not 'T' in state.keys():
            P = state['P']
            rho = state['Dmolar']
            V = 1/rho
            T = fsolve(EOSTV,x0=1,args=(V,P,aLow,bLow))
            state['T'] = T
        else: T, P = state['T'], state['P']
        coeffs = [] # Check!
        self.aLow = aLow
        self.bLow = bLow
        self.coeffs = coeffs
        self.state = state

    def SoaveEOS(self):
        def EOSTV(T,V,P0,ac,kappa,b): 
            Tc = self.Tc
            alpha = (1.0+kappa*(1.0-np.sqrt(T/Tc)))**2
            a = ac*alpha
            P = Rg*T/(V-b)-a/(V*(V+b))-P0
            return P
        state = self.state
        Tc = self.Tc
        Pc = self.Pc
        omega = self.omega
        bLow = 0.08664*Rg*Tc/Pc
        aLowC = 0.42747*(Rg*Tc)**2/Pc
        kappa = 0.480+1.574*omega-0.176*omega**2
        T, rho, P = 0, 0, 0
        if not 'P' in state.keys():
            T = state['T']
            rho = state['Dmolar']
            V = 1/rho
            P = EOSTV(T,V,0,aLowC,kappa,bLow)
            state['P'] = P
        elif not 'T' in state.keys():
            P = state['P']
            rho = state['Dmolar']
            V = 1/rho
            T = fsolve(EOSTV,x0=1,args=(V,P,aLowC,kappa,bLow))
            state['T'] = T
        else: T, P = state['T'], state['P']
        alpha = (1+kappa*(1-np.sqrt(T/Tc)))**2
        aUpp = aLowC*alpha*P/(Rg*T)**2
        bUpp = bLow*P/(Rg*T)
        coeffs = [1, -1, aUpp - bUpp - bUpp**2, -aUpp*bUpp]
        self.bLow = bLow
        self.aLowC = aLowC
        self.alpha = alpha
        self.kappa = kappa
        self.aUpp = aUpp
        self.bUpp = bUpp
        self.coeffs = coeffs
        self.state = state

    def IdealEOS(self):
        def EOSTV(T,V,P0): 
            P = Rg*T/V-P0
            return P
        state = self.state
        T, rho, P = 0, 0, 0
        if not 'P' in state.keys():
            T = state['T']
            rho = state['Dmolar']
            V = 1/rho
            P = EOSTV(T,V,0)
            state['P'] = P
        elif not 'T' in state.keys():
            P = state['P']
            rho = state['Dmolar']
            V = 1/rho
            T = fsolve(EOSTV,x0=1,args=(V,P))
            state['T'] = T
        else: T, P = state['T'], state['P']
        coeffs = [1, -1]
        self.coeffs = coeffs
        self.state = state

    def Errors(self):
        state = self.state
        eos = self.eos
        if len(state.keys()) <= 1:
            raise KeyError ('Neither pressure (P), molar density (Dmolar) nor temperature (T) was given. '
                            'Give at least two.')
        if eos == '':
            raise KeyError ('Equation of state was not defined. Equations available: Peng-Robinson, Soave, Ideal.')

    def ThermodynamicState(self,eos,**state):
        self.state = state
        self.eos = eos
        self.Errors()
        if eos == 'Peng-Robinson': self.PengRobinsonEOS()
        elif eos == 'Soave': self.SoaveEOS()
        elif eos == 'Ideal': self.IdealEOS()
        self.ZFactors()
        self.Phi()
        self.Phase()        
        self.Fugacity()
        self.ChemicalPotential()

def Help():
    print('Description: It computes the fluid\'s properties according to a given EOS.\n'
          '\tAvailable EOS:\n'
          '\t- Peng-Robinson\n'
          '\t- Soave (Soave-Redlich-Kwong)\n'
          '\t- Ideal\n'
          '\n'
          '\tAvailable properties: phi (fugacity coefficient), fugacity, mu (chemical potential),\n'
          '\t\tidealMu (ideal part of mu), excessMu (excess part of mu), phase (liquid or vapour),\n'
          '\t\tmass (fluid\'s mass), molarMass (fluid\'s molar mass), Pc (critical pressure),\n'
          '\t\tTc (critical temperature), omega (acentric factor), P (system\'s pressure),\n'
          '\t\tT (system\'s temperature), Dmolar (system\'s molar density).\n'
          '\n'          
          '\tNote: The script only works with SI units (kg, m, J, K, Pa, mol, ...).\n'
          '\n\n'
          'Instructions:\n'
          '\tFrom a script:\n'
          '\t\t1.- Import the script:\n'
          '\t\t\tExample: import EOS\n'
          '\t\t2.- Create the fluid.\n'
          '\t\t\tExample: benzene = EOS.EOS(Tc,Pc,omega,molarMass)\n'
          '\t\t3.- Call ThermodynamicState.\n'
          '\t\t\tExample: benzene.ThermodynamicState(eos=\'Peng-Robinson\',P=0.0136e6,T=298)\n'
          '\t\t\t         benzene.ThermodynamicState(eos=\'Peng-Robinson\',T=298,Dmolar=...)\n'
          '\t\t4.- Call any fluid\'s property.\n'
          '\t\t\tExample: print(benzene.phi)\n'
          '\t\t\t         print(benzene.fugacity)\n'
          '\t\t\t         print(benzene.mu)\n'
          '\n'
          '\tFrom command line:\n'
          '\t\tExecute this script.\n'
          '\t\t\tpython3 EOS.py [-flag] [value(s)]'
          '\t\tFlags: Tc (critical temperature), Pc (critical pressure), w (acentric factor), mm'
          '\t\t\t(molar mass), P (pressure), T (temperature), D (molar density), EOS (eq. od state),'
          '\t\t\tprops (fluid\'s properties to extract).\n'
          '\t\tNote:\n'
          '\t\t\t- Available properties: mu (chemical potential), idmu (ideal chemical potential), \n'
          '\t\t\t\texmu (excess chemical potential), f (fugacity), phi (fugacity coefficient), phase'
          '\t\t\t\t(fluid\'s phase).\n')

def ReadInputValues(argv):
    Pc, Tc, omega, molarMass, P, T, Dmolar, eos, properties = 0, 0, 0, 0, 0, 0, 0, '', []
    for i in range(len(argv)):
        argv[i] = argv[i].lower()
        if argv[i] == '-pc': Pc = float(argv[i+1])
        elif argv[i] == '-tc': Tc = float(argv[i+1])
        elif argv[i] == '-w': omega = float(argv[i+1])
        elif argv[i] == '-mm': molarMass = float(argv[i+1])
        elif argv[i] == '-p': P = float(argv[i+1])
        elif argv[i] == '-t': T = float(argv[i+1])
        elif argv[i] == '-d': Dmolar = float(argv[i+1])
        elif argv[i] == '-eos': eos = argv[i+1]
        elif argv[i] == '-props':
            for j in range(i+1,len(argv)): 
                if '-' in argv[j]: break
                properties.append(argv[j])
    return Pc, Tc, omega, molarMass, P, T, Dmolar, eos, properties
        
def PrintOuput(molecule,props):
    for i in range(len(props)): props[i] = props[i].lower()
    if 'mu' in props: print(molecule.mu)
    if 'idmu' in props: print(molecule.idealMu)
    if 'exmu' in props: print(molecule.excessMu)
    if 'phi' in props: print(molecule.phi)
    if 'f' in props: print(molecule.fugacity)
    if 'phase' in props: print(molecule.phase)

################################################################################
if __name__ == '__main__':
    from sys import argv

    if len(argv) > 1:
        Pc, Tc, omega, molarMass, P, T, Dmolar, eos, props = ReadInputValues(argv)
        molecule = EOS(Pc,Tc,omega,molarMass)
        if not Dmolar: molecule.ThermodynamicState(eos=eos,P=P,T=T)
        elif not P: molecule.ThermodynamicState(eos=eos,Dmolar=Dmolar,T=T)
        elif not T: molecule.ThermodynamicState(eos=eos,Dmolar=Dmolar,P=P)
        PrintOuput(molecule,props)
    else:
        # Benzene
        T_c = 562.05
        P_c = 4.894E6
        w = 0.2092
        Mm =  78.11E-3 #kg/mol

        T_0 = 298 #K
        P_0 = 1e7 #Pa 

        system = EOS(Tc=T_c, Pc=P_c, omega=w, molarMass=Mm)
        system.ThermodynamicState(eos='Peng-Robinson',T=T_0,P=P_0)
        print('Peng-Robinson')
        print('P[Pa]:',P_0)
        print('phi:',system.phi)
        print('Mu[kJ/mol]:',system.mu*1e-3)
        print('phase:',system.phase)
        # system.ThermodynamicState(eos='Soave',T=T_0,P=P_0)
        # print('phi =',system.phi)
        # print('f =',system.fugacity)
        # print('phase:',system.phase)
        # system.ThermodynamicState(eos='Ideal',T=T_0,P=P_0)
        # print('Ideal')
        # print('P[Pa]:',P_0)
        # print('Mu:',system.mu/(Rg*epsilon))

