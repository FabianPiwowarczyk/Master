import numpy as np

# constants
mdry = 28.9644e-3	#molar mass of dry air in kg
mh2o = 18.01528e-3	#molar mass of h2o in kg
mco2 = 44.01e-3	#molar mass of co2 in kg, TBC
mch4 = 16.0425e-3	#molar mass of ch4 in kg, TBC

avogadro = 6.0221367e23	#avogadro constant
boltzmann = 1.38064852e-23	#boltzmann constant
gasc_uni = avogadro * boltzmann		#universal gas constant
gasc_air = gasc_uni / mdry		#specific gas constant for dry air
gasc_h2o = gasc_uni / mh2o		#specific gas constant for water vapor
eam = 5.97219e24		#earth mass in k
ear = 6.371e6		#mean earth radius in m
gra = 6.67384e-11		#gravitational constant in m3 / kg / s2

fac = -2. * gra * eam / ear ** 3.	#free air correction constant (i.e. the reduction of g per m, linear approximation)
#   cp_air = np.double(1003.)                                            #constant for potential temperature = R/cp for dry air
#   xkappa = gasc_air / cp_air
#   p0 = np.double(1000.)                                                #reference pressure in hPa

