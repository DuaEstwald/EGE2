import numpy as np
import scipy.integrate as integrate


# =========================================
#          PUNTO 1. APARTADO 2
# =========================================


sigma_g = 2.63
sigma_c = 0.904


# Probabilidad de que se forme una galaxia en un punto escogido al azar.

P2 = lambda x2: np.exp(-(x2**2)/(2*sigma_g**2))/(np.sqrt(2*np.pi)*sigma_g)

# Integramos desde 2g hasta un numero muy grande, ya que no podemos integrar desde 2g hasta infinito numericamente, sin embargo el valor de la integral converge cuando es muy grande a un mismo valor.

P2result = integrate.quad(P2,2*sigma_g,1000)[0]


# Probabilidad de que se forme una galaxia en un punto a distancia R de donde se forme un cumulo. 

c20 = 3.23e-2
sigma_gp = sigma_g*np.sqrt(1.-c20**2)
x1 = 2*sigma_c

Pconj = lambda x2: np.exp(-((x2-(c20*sigma_g*sigma_c)*x1)**2)/(2.*sigma_gp**2))/(np.sqrt(2.*np.pi)*sigma_gp)

Pconjresult = integrate.quad(Pconj,2*sigma_g,1000)[0]

# Correlacion entre cumulos y galaxias a 20 Mpc en las condiciones iniciales.
coor = Pconjresult/P2result - 1



print('El cociente de probabilidades, apartado 1: ',coor)

# =======================================
#                PUNTO 3.
# =======================================


# a. Calcular la raiz cuadratica media del campo filtrado (top-hat) en una escala de 30hMpc 

# Si dividimos el Pk por sigma8**2 y calculamos el valor que debiera tener, para el caso de r = 8h-1Mpc



PkWTk = lambda k : ((k**2.)*(3.*(np.sin(r*k)-r*k*np.cos(r*k))/(r*k)**3.)**2.)*(19843.*(np.log(1+11.14*k)**2)*(1+18.5*k+5880*k**2+17580*k**3+1.04e6*k**4)**(-0.5)*k**(-1.))



r = 8.
dr0 = 0.83
sigma8 = dr0/np.sqrt((integrate.quad(PkWTk,0.,1000)[0])/(2.*np.pi**2.))


r = 30.

dr2 = sigma8**2.*(integrate.quad(PkWTk,0.,1000)[0])/(2.*np.pi**2.)

dr = np.sqrt(dr2)

print('sigma_8: ',sigma8)
print('Raiz cuadratica media para r = 30: ',dr)
# b. Calcular el promedio cuadratico del campo filtrado en escala r


PkWTkb = lambda k : ((3.*(np.sin(r*k)-r*k*np.cos(r*k))/(r*k)**3.)**2.)*(19843.*(np.log(1+11.14*k)**2)*(1+18.5*k+5880*k**2+17580*k**3+1.04e6*k**4)**(-0.5)*k**(-1.))



r = 1.
dr2b1 = sigma8**2.*(integrate.quad(PkWTkb,0.,1000)[0])*(((2./3.)*70.*0.3**(0.6))**2.)/(2.*np.pi**2.)

drb1 = np.sqrt(dr2b1)
print('Velocidad cuadratica media para r = 1: ',drb1)



r = 10.
dr2b2 = sigma8**2.*(integrate.quad(PkWTkb,0.,1000)[0])*(((2./3.)*70.*0.3**(0.6))**2.)/(2.*np.pi**2.)

drb2 = np.sqrt(dr2b2)
print('Velocidad cuadratica media para r = 10: ',drb2)
