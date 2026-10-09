from dataclasses import dataclass

@dataclass
class IrrepNamingPlot:
    texName: str
    
@dataclass
class QuantumNumber:
    theParticleType: str
    theIsospin: str
    theIsospinStr: str
    theStrangeness: str
    theBaryonNumber: str
    
@dataclass
class HadronParams:
    longName: str
    mass: str
    texName : str
    spinParity: str
    

the_list_irreps = {
    "G1g": IrrepNamingPlot(r'$G_{1g}$'),
    "G1u": IrrepNamingPlot(r'$G_{1u}$'),
    "G2g": IrrepNamingPlot(r'$G_{2g}$'),
    "G2u": IrrepNamingPlot(r'$G_{2u}$'),
    "Hg": IrrepNamingPlot(r'$H_{g}$'),
    "Hu": IrrepNamingPlot(r'$H_{u}$'),
    "G1": IrrepNamingPlot(r'$G_{1}$'),
    "G2": IrrepNamingPlot(r'$G_{2}$'),
    "G": IrrepNamingPlot(r'$G$'),
    "F1": IrrepNamingPlot(r'$F_{1}$'),
    "F2": IrrepNamingPlot(r'$F_{2}$'),
    "A1um": IrrepNamingPlot(r'$A_{1u}^{-}$'),
    "A1u": IrrepNamingPlot(r'$A_{1u}$'),
    "A1g": IrrepNamingPlot(r'$A_{1g}$'),
    "A2m": IrrepNamingPlot(r'$A_{2m}$'),
    "A2g": IrrepNamingPlot(r'$A_{2g}$'),
    "A1u": IrrepNamingPlot(r'$A_{2u}$'),
    "Eg": IrrepNamingPlot(r'$E_{g}$'),
    "Eu": IrrepNamingPlot(r'$E_{u}$'),
    "T1g": IrrepNamingPlot(r'$T_{1g}$'),
    "T1u": IrrepNamingPlot(r'$T_{1u}$'),
    "T2g": IrrepNamingPlot(r'$T_{2g}$'),
    "T2u": IrrepNamingPlot(r'$T_{2u}$'),
    "A1": IrrepNamingPlot(r'$A_{1}$'),
    "A2": IrrepNamingPlot(r'$A_{2}$'),
    "A2u": IrrepNamingPlot(r'$A_{2u}$'),
    "B1": IrrepNamingPlot(r'$B_{1}$'),
    "E": IrrepNamingPlot(r'$E$'),
    "B2": IrrepNamingPlot(r'$B_{2}$')}


    
the_quantum_numbers=[
    QuantumNumber('fermionic', '2I1', 'I=1/2', 'S0', 'B1'),
    QuantumNumber('fermionic', '2I1', 'I=1/2', 'Sm2', 'B1'),
    QuantumNumber('fermionic', '2I3', 'I=3/2', 'S0', 'B1'),
    QuantumNumber('fermionic', '2I3', 'I=3/2', 'Sm2', 'B1'),
    QuantumNumber('fermionic', '2I5', 'I=5/2', 'S0', 'B1'),
    QuantumNumber('fermionic', 'I0', 'I=0', 'Sm1', 'B1'),
    QuantumNumber('fermionic', 'I0', 'I=0', 'Sm3', 'B1'),
    QuantumNumber('fermionic', 'I1', 'I=1', 'Sm1', 'B1'),
    QuantumNumber('fermionic', 'I1', 'I=1', 'Sm3', 'B1'),
    QuantumNumber('fermionic', 'I2', 'I=2', 'Sm1', 'B1'),
    QuantumNumber('bosonic', 'I0', 'I=0', 'S0', 'B0'),
    QuantumNumber('bosonic', 'I0', 'I=0', 'S2', 'B0'),
    QuantumNumber('bosonic', '2I1', 'I=1/2', 'S1', 'B0'),
    QuantumNumber('bosonic', '2I3', 'I=3/2', 'S1', 'B0'),
    QuantumNumber('bosonic', '2I5', 'I=5/2', 'S1', 'B0'),
    QuantumNumber('bosonic', 'I1', 'I=1', 'S0', 'B0'),
    QuantumNumber('bosonic', 'I1', 'I=1', 'S2', 'B0'),
    QuantumNumber('bosonic', 'I2', 'I=2', 'S0', 'B0'),
    QuantumNumber('bosonic', 'I2', 'I=2', 'S2', 'B0'),
    QuantumNumber('bosonic', 'I3', 'I=3', 'S0', 'B0')]
    

    
the_list_hadrons = {
    'pi': HadronParams(longName='Pion', mass='pion_mass', texName=r'$\pi$', spinParity='$0^{-}$'),
    'K': HadronParams(longName='Kaon', mass='kaon_mass', texName=r'$K$', spinParity='$0^{-}$'),
    'KB': HadronParams(longName='K-bar', mass='kaon_mass', texName=r'$\bar{K}$', spinParity='$0^{-}$'),
    'eta': HadronParams(longName='eta', mass='eta_mass', texName=r'$\eta$', spinParity='$0^{-}$'),
    'N': HadronParams(longName='Nucleon', mass='nucleon_mass', texName=r'$N$', spinParity='$1/2^{+}$'),
    'etaprime-958': HadronParams(longName='eta-prime', mass='etaprime_mass', texName=r"$\eta'$", spinParity='$0^{-}$'),
    'Lambda': HadronParams(longName='Lambda', mass='lambda_mass', texName=r'$\Lambda$', spinParity='$1/2^{+}$'), 
    'Sigma': HadronParams(longName='Sigma', mass='sigma_mass', texName=r'$\Sigma$', spinParity='$1/2^{+}$'),
    'Delta-1232': HadronParams(longName='Delta', mass='delta_mass', texName=r'$\Delta$', spinParity='$3/2^{+}$'),
    'Xi': HadronParams(longName='Xi', mass='xi_mass', texName=r'$\Xi$', spinParity='$1/2^{+}$'),
    'Sigma-1385': HadronParams(longName='Sigma-star', mass='sigmastar_mass', texName=r'$\Sigma^{*+}$', spinParity='$3/2^{+}$'),
    'Xi-1530': HadronParams(longName='Xi-star', mass='xistar_mass', texName=r'$\Xi^{*0}$', spinParity='$3/2^{+}$'),
    'Omega': HadronParams(longName='Omega', mass='omega_mass', texName=r'$\Omega^{-}$', spinParity='$3/2^{+}$' ),
    'D': HadronParams(longName='Dmeson', mass='Dmeson_mass', texName=r'$D$', spinParity='$0^{+}$')}
