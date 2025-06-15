"""
Maxwell Bolztmann function document
"""
import numpy as np
import matplotlib.pyplot as plt
def maxwell(speed, kbt, mass=1):
    """
    Maxwell Boltz probability of given speed

    Args:
        speed: speed with which to find probability
        kbt: k*temperture called kbt at which distibution of defined
        mass: mass of particles in system

    Returns:
        p: probability of the particle having initialized speed given boltzmann
    """
    #end
    p=(mass/kbt)*speed*np.exp((-1*0.5*mass*speed*speed)/kbt)

    return p