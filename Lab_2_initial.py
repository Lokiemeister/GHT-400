# Lab 2 workspace
# we will convert our existing Excel model to Python, and begin adding modules

# importing required libraries
import numpy as np
import matplotlib.pyplot as plt

# input variables
# Controller variables
ghArea = 1  # greenhouse area, m^2
timeStep = 1000  # time per loop, seconds
nSteps = 20  # number of steps we are testing

# initial temperatures
T_glazing_init = 13
T_air_init = 10
T_ground_init = 15

# outside air
T_outside = -7  # outside temperature, degrees Celsius
RAD_Solar = 300  # solar radiation, W
h_outside = 10  # outside air convection coefficient, W/(m^2*K)
e_sky = 0.75  # sky emissivity, no units

# glazing values
x_glass = 10  # glass thickness, mm
e_glazing = 0.2  # glass emissivity, no units
tau_glazing = 0.7  # glass transmissivity, no units
rho_glazing = 2500  # glass density, kg/m^3
cp_glazing = 840  # glass specific heat capacity, J/(kg*K)

# air values
x_air = 20  # air thickness, feet
rho_air = 1.2  # air density, kg/m^3
cp_air = 1006  # air specific heat capacity, J/(kg*K)
h_insideAir = 5  # inside air convection coefficient, W/(m^2*K)

# ground values
e_ground = 1  # ground emissivity, no units
x_ground = 1  # ground thickness, m
rho_ground = 2400  # ground density, kg/m^3
cp_ground = 760  # ground specific heat capacity, J/(kg*K)

# calculated glazing values
chi_glazing = 1 - e_glazing - tau_glazing
m_glazing = x_glass / 1000 * rho_glazing

# calculated air values
x_air = x_air * 12 * 25.4 / 1000  # convert air thickness to m from feet
m_air = x_air * rho_air

# calculated ground values
m_ground = x_ground * rho_ground

# constants
Stephan_Boltzmann = 5.67E-8  # Stephan-Boltzmann constant, W/(m^2*K^4)

# initial temperatures
T_air = np.zeros(nSteps)
T_ground = np.zeros(nSteps)
T_glazing = np.zeros(nSteps)
T_air[0] = T_air_init  # greenhouse air temperature, degrees Celsius
T_ground[0] = T_ground_init  # ground temperature, degrees Celsius
T_glazing[0] = T_glazing_init  # glazing temperature, degrees Celsius


# setup of functions
def Qrad(e, A, Ts, Tsurr):
    Qrad = e * A * 5.67E-8 * ((Ts + 273.15) ** 4 - (Tsurr+273.15) ** 4)

    return Qrad


def Qconv(h, A, Ts, Tinfin):
    Qconv = h * A * (Ts - Tinfin)

    return Qconv


def T2(T1, m, cp, Q, timeStep):
    T2 = Q / (m * cp) * timeStep + T1

    return T2


# create a loop to calculate over time
# first, create loop variables. These are the same calculations performed in excel.
QRAD_sun_glazing = np.zeros(nSteps)
QRAD_sun_ground = np.zeros(nSteps)
QRAD_sun_air = np.zeros(nSteps)

QRAD_ground_glazing = np.zeros(nSteps)
QRAD_ground_sky = np.zeros(nSteps)
QRAD_ground_air = np.zeros(nSteps)

QRAD_glazing_sky = np.zeros(nSteps)

QCONV_ground_air = np.zeros(nSteps)

QCONV_glazing_air = np.zeros(nSteps)

QCONV_glazing_sky = np.zeros(nSteps)

Qsum_glazing = np.zeros(nSteps)
Qsum_air = np.zeros(nSteps)
Qsum_ground = np.zeros(nSteps)

Qvent = np.zeros(nSteps)

i = 1
while i < nSteps:
    # Radiation calculations
    QRAD_sun_glazing[i] = RAD_Solar * e_glazing
    QRAD_sun_ground[i] = RAD_Solar * tau_glazing
    QRAD_sun_air[i] = RAD_Solar - QRAD_sun_glazing[i] - QRAD_sun_ground[i]

    QRAD_ground_glazing[i] = Qrad(e_ground, ghArea, T_ground[i - 1], T_glazing[i - 1])
    QRAD_ground_sky[i] = Qrad(e_ground, ghArea, T_ground[i - 1], T_outside)
    QRAD_ground_air[i] = Qrad(e_ground, ghArea, T_ground[i - 1], T_air[i - 1])

    QRAD_glazing_sky[i] = Qrad(e_ground, ghArea, T_glazing[i - 1], T_outside)

    # create a summation of heat inputs for all temperature calculations
    Qsum_glazing[i] = QRAD_sun_glazing[i - 1] + QRAD_ground_glazing[i - 1] - QRAD_glazing_sky[i - 1] - \
                      QCONV_glazing_air[i - 1] - QCONV_glazing_sky[i - 1]
    Qsum_air[i] = QRAD_sun_air[i - 1] + QRAD_ground_air[i - 1] + QCONV_ground_air[i - 1] + QCONV_glazing_air[i - 1]
    Qsum_ground[i] = QRAD_sun_ground[i - 1] - QRAD_ground_glazing[i - 1] - QRAD_ground_sky[i - 1] - QRAD_ground_air[
        i - 1] - QCONV_ground_air[i - 1]

    T_glazing[i] = T2(T_glazing[i - 1], m_glazing, cp_glazing, Qsum_glazing[i], timeStep)
    T_air[i] = T2(T_air[i - 1], m_air, cp_air, Qsum_air[i], timeStep)
    T_ground[i] = T2(T_ground[i - 1], m_ground, cp_ground, Qsum_ground[i], timeStep)

    # Convection calculations
    QCONV_ground_air[i] = Qconv(h_insideAir, ghArea, T_ground[i], T_air[i])

    QCONV_glazing_air[i] = Qconv(h_insideAir, ghArea, T_glazing[i], T_air[i])

    QCONV_glazing_sky[i] = Qconv(h_outside, ghArea, T_glazing[i], T_outside)

    i += 1

plt.figure()
plt.plot(T_glazing)
plt.plot(T_air)
plt.plot(T_ground)
plt.legend(["Glazing", "Air", "Ground"])
plt.title('Greenhouse Temperatures')
plt.xlabel('Iteration')
plt.ylabel('Temperature (Celsius)')

plt.figure()
plt.plot(Qsum_glazing)
plt.plot(Qsum_air)
plt.plot(Qsum_ground)
plt.legend(["Glazing", "Air", "Ground"])
plt.title('Greenhouse Heat Transfer')
plt.xlabel('Iteration')
plt.ylabel('Heat Transfer (W)')

plt.show()

