unit_dict={
    "eV":1,
    "Ry":13.60570397635,
    "Ha":27.2114079527,
    "ps":1,
    "A":1,
    "nm":10,
    "hbar":6.582119569e-28, 
}

eV=0.03674929286048857
Ry=0.5 
Ha=1
ps=1 
Angstrom=1.8897261246257702
nm=10*Angstrom
bohr=1


def unit_change(elist,I0,natoms,n=1,target_unit="eV**(-1)*natoms**(-1)*ps**(-1)"):
    """
    This program changes the unit in hartree into unit of eV^{-1}natoms^{-1}ps^{-1}
    elist must be in units of hartree. 
    I0: the electric field intensity, in SI unit
    n: refractive index 
    target_units: not implemented, in the future release, researchers can write the targeted unit as python expression, and the program do the calculation automatically. 
    """
    HaeV=27.211386245988
    eps_0=8.85418782e-12
    e=1.602176634e-19
    hbar=1.054571817e-34
    HaSI=HaeV*e#Hartree in Si units 
    r0=5.29177210903e-11
    c=299792458#ms^-1
    #a0=7.291*r0/1e-9



    E=(I0*2/c/n/eps_0)**0.5#in unit of V/m, we want to change it into VH/r0, where 1V=1/27.2VH, r0=5.29177210903e-11
    EvsEh=HaeV/r0
    th=hbar/HaSI
    E_in_h=E/EvsEh
    elist*=E_in_h**2
    elist/=th
    elist*=1e-12
    elist/=natoms 

    #Now let's look at the units:
    

    eV=1 
    Ha=27.211386245988
    ps=1e-12
    ns=1e-9  



    return elist

