# Author: Santiago A. Flores Roman
#
# Description: It computes the fluid's properties according to a given EOS.
#   Available EOS:
#   - Peng-Robinson
#   - Soave (Soave-Redlich-Kwong)
#   - Ideal (Ideal gas)
#   - Johnson (Lennard-Jones)
#   - VdW (Van der Waals)
#
#   Available properties: phi (fugacity coefficient), fugacity, mu (chemical potential),
#       idealMu (ideal part of mu), excessMu (excess part of mu), phase (liquid or vapour),
#       mass (fluid's mass), molarMass (fluid's molar mass), Pc (critical pressure),
#       Tc (critical temperature), omega (acentric factor), P (system's pressure),
#       T (system's temperature), Dmolar (system's molar density), P_sat (saturation pressure),
#       bulkModulus (bulk modulus or reciprocal of isothermal compressibility),
#       zFactor (compressibility factor), roots (solutions of EOS [fluid's densities]).
#
#   Note 1: The script uses the following units: kg, kg/mol, m, J/mol, K (only for T, and Tc), Pa.
#   Note 2: As the EOSs implemented don't predict exactly the saturation pressure, the user can
#       input this value when calling ThermodynamicState. See Example 1.
#
# Instructions:
#   1.- Import the script.
#       Example: import EOS
#   2.- Create the fluid.
#       Example 1:
#           Tcrit = critical_temperature #K
#           Pcrit = critical_pressure #Pa
#           w = acentric_factor
#           mM = molar_mass #kg/mol
#           benzene = EOS.EOS(Tc=Tcrit, Pc=Pcrit, omega=w, molarMass=mM)
#       Example 2:
#           eps = epsilon #K
#           sig = sigma #m
#           mM = molar_mass #kg/mol
#       benzene = EOS.EOS(epsilon=eps, sigma=sig, molarMass=mM)
#   3.- Call ThermodynamicState.
#       Example 1:
#           benzene.ThermodynamicState(eos='Peng-Robinson', P=0.0136e6, T=298)
#           benzene.ThermodynamicState(eos='Ideal Gas', T=300, Dmolar=...)
#           benzene.ThermodynamicState(eos='Soave', P=300, Dmolar=...)
#           benzene.ThermodynamicState(eos='VdW', P=1e5, T=..., P_sat=101235)
#       Example 2:
#           benzene.ThermodynamicState(eos='Johnson',P=1e5,T=...)
#   4.- Call any fluid's properties.
#       Example: print(benzene.phi)
#                print(benzene.fugacity) #Pa
#                print(benzene.mu) #J/mol
#                print(benzene.idealMu) #J/mol
#                print(benzene.excessMu) #J/mol
#                print(benzene.P_sat) #Pa
#                print(benzene.roots) #mol/m^3
#                print(benzene.Dmolar) #mol/m^3
#                print(benzene.P) #Pa
#                print(benzene.T) #K
#                print(benzene.bulkModulus) #Pa
#                print(benzene.zFactor)
#                print(benzene.phase)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import scipy.constants as const
from scipy.optimize import root
import numpy as np
from numpy.polynomial import Polynomial

Rg = const.R
kb = const.Boltzmann
planck = const.h
avogadro = const.Avogadro

class EOS():
    def __init__(self,Tc=0,Pc=0,omega=0,molarMass=0,sigma=0,epsilon=0):
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
        self.P, self.P_sat, self.bulkModulus, self.T, self.Dmolar = 0, 0, 0, 0, 0 #Pa, Pa, Pa, K, mol/m^3
        self.roots = 0 #mol/m^3
        #For Johnson EOS---------------------------------------------------------------------------
        self.epsilon = epsilon #K
        self.sigma = sigma #m
        self.gamma, self.rcut = 3, 5
        self.aCoeff, self.bCoeff, self.cCoeff = np.zeros(9), np.zeros(7), np.zeros(9)
        self.dCoeff, self.gCoeff = np.zeros(7), np.zeros(7)
        self.xCoeff = [0.0, 0.8623085097507421, 2.976218765822098, -8.402230115796038,
                       0.1054136629203555, -0.8564583828174598, 1.582759470107601,
                       0.7639421948305453, 1.753173414312048, 2.798291772190376e3,
                       -4.8394220260857657e-2, 0.9963265197721935, -3.698000291272493e1,
                       2.084012299434647e1, 8.305402124717285e1, -9.574799715203068e2,
                       -1.477746229234994e2, 6.398607852471505e1, 1.603993673294834e1,
                       6.805916615864377e1, -2.791293578795945e3, -6.245128304568454,
                       -8.116836104958410e3, 1.488735559561229e1, -1.059346754655084e4,
                       -1.131607632802822e2, -8.867771540418822e3, -3.986982844450543e1,
                       -4.689270299917261e3, 2.593535277438717e2, -2.694523589434903e3,
                       -7.218487631550215e2, 1.721802063863269e2]
        self.helmholtz, self.rho = [], []

    def ZFactors(self):
        eos = self.eos
        if eos == 'Johnson':
            roots = self.rho
            state = self.state
            epsilon, sigma = self.epsilon, self.sigma #K, m^3
            P, T = state['P']*sigma**3/(kb*epsilon), state['T']/epsilon
            zFactor = P/(roots*T)
        else:
            coeffs = self.coeffs
            zFactor = coeffs.roots()
        zFactor = zFactor[np.imag(zFactor) == 0]
        zFactor = np.sort(np.real(zFactor))
        self.zFactor = zFactor

    def Phi(self):
        eos = self.eos
        zFactor = self.zFactor
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
        elif eos == 'VdW':
            np.seterr(invalid='ignore')
            phi = np.exp(bUpp/(zFactor-bUpp)-2*aUpp/zFactor-np.log(1-aUpp*(zFactor-bUpp)/zFactor**2)) #Check!
        elif eos == 'Ideal': phi = zFactor
        elif eos == 'Johnson':
            helmholtz = self.helmholtz
            state = self.state
            epsilon = self.epsilon #m, K
            zFactor = self.zFactor
            T = state['T']/epsilon
            muResidual = helmholtz+zFactor*T-T
            phi = np.exp(muResidual/T)
        phi = np.nan_to_num(phi,nan=1e3)
        self.phi = phi

    def FindSaturationPressure(self):
        pSat = self.P_sat #Initially, pSat=1e-5. If pSat=0, then pSat is not computed.
        eos = self.eos
        if pSat == 0: return
        if eos == 'Peng-Robinson': self.PengRobinsonEOS()
        elif eos == 'VdW': self.VdWEOS()
        elif eos == 'Soave': self.SoaveEOS()
        elif eos == 'Ideal': self.IdealEOS()
        elif eos == 'Johnson': self.JohnsonEOS()
        self.ZFactors()
        self.Phi()
        if eos == 'Johnson':
            # This piece of code was developed by Maximov, Max.
            def P_LJ(Rho):
                epsilon = self.epsilon #K
                aCoeff, bCoeff = self.aCoeff, self.bCoeff
                rcut, gamma = self.rcut, self.gamma
                T = self.state['T']/epsilon
                fExp = np.exp(-gamma*Rho**2)
                sum8, sum6 = 0, 0
                for i in range(1,9): sum8 += aCoeff[i]*Rho**(i+1)
                for i in range(1,7): sum6 += bCoeff[i]*Rho**(2*i+1)
                pFull = Rho*T+sum8+fExp*sum6
                pTail = 32/9 * np.pi*Rho**2 * ((1/rcut)**9-3/2*(1/rcut)**3)
                P = pFull-pTail
                return P
            def Helmholtz(Rho):
                aCoeff, bCoeff, gCoeff = self.aCoeff, self.bCoeff, self.gCoeff
                rcut, gamma = self.rcut, self.gamma
                fExp = np.exp(-gamma*Rho**2)
                gCoeff = [0] * 7
                gCoeff[1] = (1-fExp)/(2*gamma)
                gCoeff[2] = -(fExp*Rho**2 - 2*gCoeff[1])/(2*gamma)
                gCoeff[3] = -(fExp*Rho**4 - 4*gCoeff[2])/(2*gamma)
                gCoeff[4] = -(fExp*Rho**6 - 6*gCoeff[3])/(2*gamma)
                gCoeff[5] = -(fExp*Rho**8 - 8*gCoeff[4])/(2*gamma)
                gCoeff[6] = -(fExp*Rho**10 - 10*gCoeff[5])/(2*gamma)
                sum8, sum6 = 0, 0
                for i in range(1,9): sum8 += aCoeff[i]*Rho**i/i
                for i in range(1,7): sum6 += bCoeff[i]*gCoeff[i]
                helmholtzFull = sum8+sum6
                helmholtzTail = 32/9 * np.pi * Rho * ((1/rcut)**9 - 3/2*(1/rcut)**3)
                helmholtz = helmholtzFull-helmholtzTail
                return helmholtz
            def Mu_LJ(Rho):
                epsilon = self.epsilon #K
                sigma = self.sigma #m
                mass = self.mass #kg
                helmholtz = Helmholtz(Rho)
                T = self.state['T']/epsilon
                P = P_LJ(Rho)
                thermalWL = planck/np.sqrt(2*np.pi*mass*kb*T) #m
                thermalWL /= sigma
                idMu = T*np.log(Rho*thermalWL**3)
                exMu = helmholtz+P/Rho-T
                return idMu+exMu
            def logspace(start, stop, num):
                log10 = np.log(10.0)
                delta_dens = (np.log(stop) - np.log(start)) / (num - 1) / log10
                return np.array([pow(10.0, np.log(start) / log10 + delta_dens * (i - 1)) for i in range(num)])
            def gas_spinodal_point(mu, rho):
                inflection_idx = np.argwhere(np.diff(mu(rho)) < 0)
                if len(inflection_idx) == 0:
                    raise ValueError("No spinodal points. Check your input parameters "\
                                      "or give a saturation pressure among your input parameters.")
                    pSat = 0
                    self.pSat = pSat
                return rho[inflection_idx[0]][0]
            def liq_spinodal_point(mu, rho):
                rho_inv = rho[::-1]
                inflection_idx = np.argwhere(np.diff(mu(rho_inv)) > 0)
                if len(inflection_idx) == 0:
                    raise ValueError("No spinodal points. Check your input parameters.")
                    pSat = 0
                    self.pSat = pSat
                return rho_inv[inflection_idx[0]][0]
            epsilon = self.epsilon #K
            sigma = self.sigma #m
            rho = np.linspace(1e-4,0.9,16000)
            gas_spinodal_rho = gas_spinodal_point(Mu_LJ, rho)
            liq_spinodal_rho = liq_spinodal_point(Mu_LJ, rho)
            # min_rho = rho[0]
            max_rho = rho[-1]
            for i_pass in range(1, 3):
                n_points = n_gas_points = n_liq_points = n_interp = len(rho)
                # we divide the gas branch logarithmically and the liquid branch linearly
                if i_pass == 1: gas_rho = logspace(rho[0], gas_spinodal_rho, n_gas_points)
                else: gas_rho = np.linspace(rho[0], gas_spinodal_rho, n_gas_points)
                gas_mu = Mu_LJ(gas_rho)
                gas_p = P_LJ(gas_rho)
                # calculate chem.potential on the liquid branch
                liq_rho = np.linspace(liq_spinodal_rho, max_rho, n_liq_points)
                liq_mu = Mu_LJ(liq_rho)
                liq_p = P_LJ(liq_rho)
                # Determination of mu0, p0, satur_liq_dens and satur_gas_dens
                # here we consider ONLY 2-phase region, between mu_liq[0] and mu_gas[-1]
                interp_mu = np.linspace(liq_mu[0], gas_mu[-1], n_interp)
                gas_interp_p = np.interp(interp_mu, gas_mu, gas_p)
                liq_interp_p = np.interp(interp_mu, liq_mu, liq_p)
                gas_interp_rho = np.interp(interp_mu, gas_mu, gas_rho)
                liq_interp_rho = np.interp(interp_mu, liq_mu, liq_rho)
                # Solving p_gas_spl[mu] - p_liq_spl[mu] = 0
                i = np.argwhere(gas_interp_p - liq_interp_p < 0)[0]
                p_diff = gas_interp_p[i] - liq_interp_p[i]
                prev_dens = gas_interp_p[i - 1] - liq_interp_p[i - 1]
                if p_diff < 0 and prev_dens > 0:
                    i_eq = i if (abs(p_diff) < abs(prev_dens)) else i - 1
                    i_eq = i_eq[0]
                    # muSat = interp_mu[i_eq]
                    # dens = np.linspace(min_rho, max_rho, n_points)
                    # mu0_intersection_i = np.argwhere(np.diff(np.sign(Mu_LJ(dens) - muSat)) != 0)

                    # rho0 = dens[mu0_intersection_i[1]]
                    # densities = [gas_interp_rho[i_eq], rho0, liq_interp_rho[i_eq]],
                    pSat = liq_interp_p[i_eq]

                    # min_rho = gas_interp_rho[i_eq - n_points // 300]
                    max_rho = liq_interp_rho[i_eq + n_points // 300]
                    gas_spinodal_rho = gas_interp_rho[i_eq + n_points // 300]
                    liq_spinodal_rho = liq_interp_rho[i_eq - n_points // 300]
                else:
                    raise ValueError("Unable to find vapor-liquid equilibrium on the %d pass." % i_pass)
                    pSat = 0
            pSat *= kb*epsilon/sigma**3 #Pa
        else:
            targetP, pSat = self.state['P'], self.P_sat #Pa, Pa
            tolerance = 1e-2
            for trialP in np.linspace(1e-5,1e9,int(1e7)):
                self.state['P'] = trialP #Pa
                if eos == 'Peng-Robinson': self.PengRobinsonEOS()
                elif eos == 'VdW': self.VdWEOS()
                elif eos == 'Soave': self.SoaveEOS()
                elif eos == 'Ideal': self.IdealEOS()
                elif eos == 'Johnson': self.JohnsonEOS()
                self.ZFactors()
                self.Phi()
                phi = self.phi
                if (len(phi) == 1):
                    print('Warning: Unable to find saturation pressure.\n'\
                          '\tTry using a lower temperature.')
                    pSat = 0
                    break
                if (abs(phi[0]-phi[-1]) < tolerance):
                    pSat = trialP; break #Pa
            self.state['P'] = targetP #Pa
        self.P_sat = pSat #Pa

    def Phase(self):
        zFactor = self.zFactor
        phi = self.phi
        state = self.state
        Tc = self.Tc
        Pc = self.Pc
        T, P, pSat = state['T'], state['P'], self.P_sat #K, Pa, Pa
        roots = P/(zFactor*Rg*T) #mol/m^3
        if pSat != 0:
            if (len(zFactor) == 1) or (zFactor[0] <= 0):
                phi = phi[-1]
                zFactor = zFactor[-1]
                if (T > Tc) and (P > Pc): phase = 'super_critical_fluid'
                elif (T < Tc) and (P < Pc): phase = 'Vapour'
                elif (T < Tc) and (P > Pc): phase =  'Liquid'
            else:
                if P < pSat: #if phi[0] > phi[-1]:
                    phi, zFactor = phi[-1], zFactor[-1]
                    phase = 'Vapour (stable) - Liquid (metastable)'
                # Vapour (metastable), liquid (stable)
                elif P > pSat: #elif phi[0] < phi[-1]:
                    phi, zFactor = phi[0], zFactor[0]
                    phase = 'Vapour (metastable) - Liquid (stable)'
                # Vapour (stable), liquid (stable)
                else:
                    phi, zFactor = phi[0], zFactor[0]
                    phase = 'Vapour (stable) - Liquid (stable)'
        else: phase = 'Unpredicted'
        if not 'Dmolar' in state.keys():
            rhoMolar = P/(zFactor*Rg*T) #mol/m^3
            state['Dmolar'] = rhoMolar
        rhoN = P/(zFactor*kb*T) #m^-3
        state['Dn'] = rhoN
        self.phi = phi
        self.zFactor = zFactor
        self.phase = phase
        self.state = state
        self.P, self.T, self.Dmolar = state['P'], state['T'], state['Dmolar']
        self.roots = roots

    def ChemicalPotential(self):
        state = self.state
        phi = self.phi
        mass = self.mass #kg
        T, P, rhoN = state['T'], state['P'], state['Dn']
        fugacity = phi*P #Pa
        thermalWL = planck/np.sqrt(2*np.pi*mass*kb*T) #m
        # Ideal Mu uses the density calculated from an EOS.
        idealMu = Rg*T*np.log(rhoN*thermalWL**3) #J/mol
        excessMu = Rg*T*np.log(phi) #J/mol
        mu = idealMu+excessMu #J/mol
        self.idealMu = idealMu
        self.excessMu = excessMu
        self.mu = mu
        self.fugacity = fugacity #Pa

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
            T = root(EOSTV,x0=1,args=(V,P,aLowC,kappa,bLow)).x[0]
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
        state = self.state
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
            T = root(EOSTV,x0=1,args=(V,P,aLow,bLow)).x[0]
            state['T'] = T
        else: T, P = state['T'], state['P']
        aUpp = aLow*P/(Rg*T)**2
        bUpp = bLow*P/(Rg*T)
        coeffs = Polynomial([-aUpp*bUpp, aUpp, -(1+bUpp), 1]) # Check!
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
            T = root(EOSTV,x0=1,args=(V,P,aLowC,kappa,bLow)).x[0]
            state['T'] = T
        else: T, P = state['T'], state['P']
        alpha = (1+kappa*(1-np.sqrt(T/Tc)))**2
        aUpp = aLowC*alpha*P/(Rg*T)**2
        bUpp = bLow*P/(Rg*T)
        coeffs = Polynomial([-aUpp*bUpp, aUpp - bUpp - bUpp**2, -1, 1])
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
            T = root(EOSTV,x0=1,args=(V,P)).x[0]
            state['T'] = T
        else: T, P = state['T'], state['P']
        coeffs = Polynomial([-1, 1])
        self.coeffs = coeffs
        self.state = state

    # Information about the EOS:
    # Johnson, J. K., Zollweg, J. A., & Gubbins, K. E. (1993).
    # The Lennard-Jones equation of state revisited. Molecular Physics, 78(3), 591-618.
    #
    # All the expressions in this EOS must be manipulated in reduced units.
    def JohnsonEOS(self):
        def EOSTRho(T,Rho,P0=0):
            aCoeff, bCoeff, cCoeff, dCoeff = self.aCoeff, self.bCoeff, self.cCoeff, self.dCoeff
            xCoeff, rcut, gamma = self.xCoeff, self.rcut, self.gamma
            # Table 5
            aCoeff[1] = xCoeff[1]*T + xCoeff[2]*np.sqrt(T) + xCoeff[3] + xCoeff[4]/T + xCoeff[5]/T**2
            aCoeff[2] = xCoeff[6]*T + xCoeff[7] + xCoeff[8]/T + xCoeff[9]/T**2
            aCoeff[3] = xCoeff[10]*T + xCoeff[11] + xCoeff[12]/T
            aCoeff[4] = xCoeff[13]
            aCoeff[5] = xCoeff[14]/T + xCoeff[15]/T**2
            aCoeff[6] = xCoeff[16]/T
            aCoeff[7] = xCoeff[17]/T + xCoeff[18]/T**2
            aCoeff[8] = xCoeff[19]/T**2
            # Table 6
            bCoeff[1] = xCoeff[20]/T**2 + xCoeff[21]/T**3
            bCoeff[2] = xCoeff[22]/T**2 + xCoeff[23]/T**4
            bCoeff[3] = xCoeff[24]/T**2 + xCoeff[25]/T**3
            bCoeff[4] = xCoeff[26]/T**2 + xCoeff[27]/T**4
            bCoeff[5] = xCoeff[28]/T**2 + xCoeff[29]/T**3
            bCoeff[6] = xCoeff[30]/T**2 + xCoeff[31]/T**3 + xCoeff[32]/T**4
            fExp = np.exp(-gamma*Rho**2)
            sum8, sum6 = 0, 0
            for i in range(1,9): sum8 += aCoeff[i]*Rho**(i+1)
            for i in range(1,7): sum6 += bCoeff[i]*Rho**(2*i+1)
            pFull = Rho*T+sum8+fExp*sum6
            pTail = 32/9 * np.pi*Rho**2 * ((1/rcut)**9-3/2*(1/rcut)**3)
            P = pFull-pTail-P0
            self.fExp = fExp
            self.aCoeff, self.bCoeff, self.cCoeff, self.dCoeff = aCoeff, bCoeff, cCoeff, dCoeff
            return P
        def Helmholtz(Rho):
            aCoeff, bCoeff, gCoeff = self.aCoeff, self.bCoeff, self.gCoeff
            rcut, gamma = self.rcut, self.gamma
            fExp = np.exp(-gamma*Rho**2)
            gCoeff[1] = (1-fExp)/(2*gamma)
            gCoeff[2] = -(fExp*Rho**2 - 2*gCoeff[1])/(2*gamma)
            gCoeff[3] = -(fExp*Rho**4 - 4*gCoeff[2])/(2*gamma)
            gCoeff[4] = -(fExp*Rho**6 - 6*gCoeff[3])/(2*gamma)
            gCoeff[5] = -(fExp*Rho**8 - 8*gCoeff[4])/(2*gamma)
            gCoeff[6] = -(fExp*Rho**10 - 10*gCoeff[5])/(2*gamma)
            sum8, sum6 = 0, 0
            for i in range(1,9): sum8 += aCoeff[i]*Rho**i/i
            for i in range(1,7): sum6 += bCoeff[i]*gCoeff[i]
            helmholtzFull = sum8+sum6
            helmholtzTail = 32/9 * np.pi * Rho * ((1/rcut)**9 - 3/2*(1/rcut)**3)
            helmholtz = helmholtzFull-helmholtzTail
            self.gCoeff = gCoeff
            return helmholtz
        state = self.state
        state = self.state
        epsilon = self.epsilon #K
        sigma = self.sigma #m
        Tc, Pc = 1.313, 0.13 #K, Pa
        T, rho, P = 0, 0, 0
        roots, helmholtz = [], []
        if not 'P' in state.keys():
            T = state['T']/epsilon
            rho = state['Dmolar']*avogadro*sigma**3
            P = EOSTRho(T,rho)
            state['P'] = P*kb*epsilon/sigma**3 #Pa
            roots.append(rho)
        elif not 'T' in state.keys():
            P = state['P']*sigma**3/(kb*epsilon)
            rho = state['Dmolar']*avogadro*sigma**3
            T = root(EOSTRho,x0=1,args=(rho,P)).x[0]
            state['T'] = T*epsilon #K
            roots.append(rho)
        else:
            def EOSRhoT(Rho,T,P0=0):
                aCoeff, bCoeff, cCoeff, dCoeff = self.aCoeff, self.bCoeff, self.cCoeff, self.dCoeff
                xCoeff, rcut, gamma = self.xCoeff, self.rcut, self.gamma
                # Table 5
                aCoeff[1] = xCoeff[1]*T + xCoeff[2]*np.sqrt(T) + xCoeff[3] + xCoeff[4]/T + xCoeff[5]/T**2
                aCoeff[2] = xCoeff[6]*T + xCoeff[7] + xCoeff[8]/T + xCoeff[9]/T**2
                aCoeff[3] = xCoeff[10]*T + xCoeff[11] + xCoeff[12]/T
                aCoeff[4] = xCoeff[13]
                aCoeff[5] = xCoeff[14]/T + xCoeff[15]/T**2
                aCoeff[6] = xCoeff[16]/T
                aCoeff[7] = xCoeff[17]/T + xCoeff[18]/T**2
                aCoeff[8] = xCoeff[19]/T**2
                # Table 6
                bCoeff[1] = xCoeff[20]/T**2 + xCoeff[21]/T**3
                bCoeff[2] = xCoeff[22]/T**2 + xCoeff[23]/T**4
                bCoeff[3] = xCoeff[24]/T**2 + xCoeff[25]/T**3
                bCoeff[4] = xCoeff[26]/T**2 + xCoeff[27]/T**4
                bCoeff[5] = xCoeff[28]/T**2 + xCoeff[29]/T**3
                bCoeff[6] = xCoeff[30]/T**2 + xCoeff[31]/T**3 + xCoeff[32]/T**4
                fExp = np.exp(-gamma*Rho**2)
                sum8, sum6 = 0, 0
                for i in range(1,9): sum8 += aCoeff[i]*Rho**(i+1)
                for i in range(1,7): sum6 += bCoeff[i]*Rho**(2*i+1)
                pFull = Rho*T+sum8+fExp*sum6
                pTail = 32/9 * np.pi*Rho**2 * ((1/rcut)**9-3/2*(1/rcut)**3)
                P = pFull-pTail-P0
                self.aCoeff, self.bCoeff, self.cCoeff, self.dCoeff = aCoeff, bCoeff, cCoeff, dCoeff
                return P
            state = self.state
            epsilon, sigma = self.epsilon, self.sigma #K, m
            P, T = state['P']*sigma**3/(kb*epsilon), state['T']/epsilon
            linspace = np.linspace(0,1,10,endpoint=False)
            for i in range(len(linspace)): #Check!
                findRoots = root(EOSRhoT,x0=linspace[i],args=(T,P))
                if findRoots.success:
                    roots.append(round(findRoots.x[0],9))
            roots = np.unique(np.array(roots))
        for i in roots: helmholtz.append(Helmholtz(i))
        self.helmholtz = np.sort(helmholtz)
        self.rho = np.array(roots)
        Tc *= epsilon #K
        Pc *= kb*epsilon/sigma**3 #Pa
        self.Tc, self.Pc = Tc, Pc
        self.state = state

    def BulkModulus(self):
        eos = self.eos
        state = self.state
        bulkModulus = 0
        if eos == 'Ideal':
            P = state['P'] #Pa
            bulkModulus = P
        elif eos == 'VdW': bulkModulus = -1 #Check!
        elif eos == 'Soave': bulkModulus = -1 #Check!
        elif eos == 'Peng-Robinson':
            aLowC, alpha, bLow = self.aLowC, self.alpha, self.bLow
            aLow = aLowC*alpha
            T, rho = state['T'], state['Dmolar'] #K, mol/m^3
            Vm = 1/rho #m^3/mol
            bulkModulus = Vm*(Rg*T/(Vm-bLow)**2 - 2*aLow*(Vm+bLow) / (Vm*(Vm+bLow)+bLow*(Vm-bLow))**2) #Pa
        elif eos == 'Johnson':
            epsilon, sigma = self.epsilon, self.sigma
            aCoeff, bCoeff, gamma = self.aCoeff, self.bCoeff, self.gamma
            P = state['P']*sigma**3/(kb*epsilon)
            T = state['T']/epsilon
            rho = state['Dmolar']*avogadro*sigma**3
            fExp = np.exp(-gamma*rho**2)
            sum8, sum6 = 0, 0
            for i in range(1,9): sum8 += i*aCoeff[i]*rho**(i+1)
            for i in range(1,7): sum6 += bCoeff[i]*(i-gamma*rho**2)*rho**(2*i+1)
            bulkModulus = P + sum8 + 2*fExp*sum6
            bulkModulus *= kb*epsilon/sigma**3 #Pa
        self.bulkModulus = bulkModulus

    def Errors(self):
        state = self.state
        eos = self.eos
        if len(state.keys()) <= 1:
            raise KeyError ('Neither pressure (P), molar density (Dmolar) nor temperature (T) was given. '
                            'Give at least two.')
        if eos == '':
            raise KeyError ('Equation of state was not defined. Equations available: Peng-Robinson, Soave, Ideal.')

    def ThermodynamicState(self, eos, P_sat=1e-5, **state):
        self.state = state
        self.eos = eos
        self.P_sat = P_sat #Pa
        self.Errors()
        if P_sat == 1e-5: self.FindSaturationPressure()
        if eos == 'Peng-Robinson': self.PengRobinsonEOS()
        elif eos == 'VdW': self.VdWEOS() #Check!
        elif eos == 'Soave': self.SoaveEOS() #Check!
        elif eos == 'Ideal': self.IdealEOS()
        elif eos == 'Johnson': self.JohnsonEOS()
        self.ZFactors()
        self.Phi()
        self.Phase()
        self.ChemicalPotential()
        self.BulkModulus()

def Help():
    print('Description: It computes the fluid\'s properties according to a given EOS.\n'
          '\tAvailable EOS:\n'
          '\t- Peng-Robinson\n'
          '\t- Soave (Soave-Redlich-Kwong)\n'
          '\t- VdW (Van der Waals)\n'
          '\t- Johnson (for Lennard-Jones fluids)\n'
          '\t- Ideal (Ideal Gas)\n'
          '\n'
          '\tAvailable properties: phi (fugacity coefficient), fugacity, mu (chemical potential),\n'
          '\t\tidealMu (ideal part of mu), excessMu (excess part of mu), phase (liquid or vapour),\n'
          '\t\tmass (fluid\'s mass), molarMass (fluid\'s molar mass), Pc (critical pressure),\n'
          '\t\tTc (critical temperature), omega (acentric factor), P (system\'s pressure),\n'
          '\t\tT (system\'s temperature), Dmolar (system\'s molar density),\n'
          '\t\tbulkModulus (fluid\'s bulk modulus).\n'
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
          '\t\t\t         print(benzene.mu)\n')

################################################################################
if __name__ == '__main__':
    from CoolProp.CoolProp import PropsSI, PhaseSI
    molName = 'Argon'
    T_0 = 87.3 #K
    P_0 = 2*PropsSI('P','T',T_0,'Q',1,molName) #Pa. At T_0 = 119.6 K. Saturation pressure.
    Dm_0 = PropsSI('Dmolar','T',T_0,'P',P_0,molName) #kg/mol. At T_0 = 119.6 K.
    Psat = PropsSI('P','T',T_0,'Q',1,molName)
    Pc = PropsSI('Pcrit',molName) #Pa
    Tc = PropsSI('Tcrit',molName) #K
    omega = PropsSI('acentric',molName)
    Mm = 39.948e-3 #kg/mol
    epsilon = 119.6 #K
    sigma = 0.34e-9 #m
    print('P[Pa]:',P_0)
    print('\nCoolProp') #----------------------------------------------------------
    print('Psat[Pa]:',Psat)
    print('Z:',PropsSI('Z','T',T_0,'P',P_0,molName))
    print('Dmolar[mol/m^3]:',PropsSI('Dmolar','T',T_0,'P',P_0,molName))
    print('phase:',PhaseSI('T',T_0,'P',P_0,molName))
    print('K_T[GPa]:',1e-9/PropsSI('ISOTHERMAL_COMPRESSIBILITY','T',T_0,'P',P_0,molName))
    print('P:',round(PropsSI('P','T',T_0,'Dmolar',Dm_0,molName)*sigma**3/(kb*epsilon),4))
    system = EOS(Pc=Pc, Tc=Tc, omega=omega, epsilon=epsilon, sigma=sigma, molarMass=Mm)
    print('\nSoave EOS') #-------------------------------------------------------
    system.ThermodynamicState(eos='Soave', T=T_0, P=P_0)
    print('Psat[Pa]:',system.P_sat)
    print('phi:',system.phi)
    print('Mu[kJ/mol]:',system.mu*1e-3)
    print('Z:',system.zFactor)
    print('Dmolar[mol/m^3]:',system.Dmolar)
    print('phase:',system.phase)
    print('Roots[mol/m^3]:',system.roots)
    print('\nPeng-Robinson EOS') #-----------------------------------------------
    system.ThermodynamicState(eos='Peng-Robinson', T=T_0, P=P_0)
    print('Psat[Pa]:',system.P_sat)
    print('phi:',system.phi)
    print('Mu[kJ/mol]:',system.mu*1e-3)
    print('Z:',system.zFactor)
    print('Dmolar[mol/m^3]:',system.Dmolar)
    print('phase:',system.phase)
    print('K_T[GPa]:',system.bulkModulus*1e-9)
    print('Roots[mol/m^3]:',system.roots)
    print('\nVan der Waals EOS') #-----------------------------------------------
    system.ThermodynamicState(eos='VdW', P=P_0, T=T_0, P_sat=Psat)
    print('Psat[Pa]:',system.P_sat)
    print('phi:',system.phi)
    print('Mu[kJ/mol]:',system.mu*1e-3)
    print('Z:',system.zFactor)
    print('Dmolar[mol/m^3]:',system.Dmolar)
    print('phase:',system.phase)
    print('Roots[mol/m^3]:',system.roots)
    print('\nJohnson EOS') #-----------------------------------------------------
    system.ThermodynamicState(eos='Johnson', P=P_0, T=T_0)
    print('Psat[Pa]:',system.P_sat)
    print('phi:',system.phi)
    print('Mu[kJ/mol]:',system.mu*1e-3)
    print('Z:',system.zFactor)
    print('Dmolar[mol/m^3]:',system.Dmolar)
    print('phase:',system.phase)
    print('K_T[GPa]:',system.bulkModulus*1e-9)
    print('Roots[mol/m^3]:',system.roots)
    print('\nIdeal Gas EOS') #---------------------------------------------------
    system.ThermodynamicState(eos='Ideal', T=T_0, P=P_0, P_sat=Psat)
    print('Psat[Pa]:',system.P_sat)
    print('phi:',system.phi)
    print('Mu[kJ/mol]:',system.mu*1e-3)
    print('Z:',system.zFactor)
    print('Dmolar[mol/m^3]:',system.Dmolar)
    print('phase:',system.phase)
    print('K_T[GPa]:',system.bulkModulus*1e-9)
    print('Roots[mol/m^3]:',system.roots)
