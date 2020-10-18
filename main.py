import numpy as np
from scipy.integrate import odeint, quad

from matplotlib import pyplot as plt
import matplotlib as mpl

from excel import *

# initial data
N_0, S_0, I_0, R_0 = get_initial_data()
N_0 = S_0 + I_0 + R_0
y_0 = (N_0, S_0, I_0, R_0)

# parametros iniciales para regresion exponencial
beta_0 = 5e-3/N_0
mu_star_0 = 2.2e-4

# Calculamos beta y mu_star segun los datos del excel
nuevos_casos = get_datos_diarios('C')
nuevos_fallecimientos = get_datos_diarios('D')
nuevas_recuperaciones = get_datos_diarios('E')

S = S_0
I = I_0
R = R_0
N = N_0

data_length = len(nuevos_casos)

beta_obs = []
mu_star_obs = []
gamma_obs = []

for i in range(data_length):
    beta_obs.append(nuevos_casos[i]/(I*N_0))
    mu_star_obs.append(nuevos_fallecimientos[i]/I)
    gamma_obs.append(nuevas_recuperaciones[i]/I)

    # datos acumulados
    S += - nuevos_casos[i]
    I += nuevos_casos[i] - nuevos_fallecimientos[i] - nuevas_recuperaciones[i]
    R += nuevas_recuperaciones[i]
    N += - nuevos_fallecimientos[i]

    # checkiamos que vaya bien
    assert (N == S + I + R)


def reg_log(x, y):
    b = np.sum(np.multiply(x, np.log(y))) - np.log(y).mean() * np.sum(x)
    b /= np.sum(np.square(x)) - x.mean() * np.sum(x)

    a = np.exp(np.log(y).mean() - b * x.mean())

    return b, np.log(a)


# calculamos los parametros de la regresion exponencial
x = np.array(range(data_length))
a, b = reg_log(x, np.array(beta_obs))
c, d = reg_log(x, np.array(mu_star_obs))

# intervalo de tiempo en el que hacemos predicciones
linspace = range(241)


# estas son las curvitas
def beta(t):
    if t <= 150:
        return beta_0 + np.exp(a*t + b)
    elif t <= 180:
        return beta(150) + (10*beta_0 - beta(150))*(t - 150)/30
    else:
        return 10*beta_0


def mu_star(t):
    if t <= 150:
        return mu_star_0 + np.exp(c*t + d)
    elif t <= 180:
        return mu_star(150) + (2*mu_star_0 - mu_star(150))*(t - 150)/30
    else:
        return 2*mu_star_0


def mu(t):
    return 3.08e-5 + 5.4e-6*np.cos(1.67e-2*t + 0.735)


# Sacamos el valor de gamma, media de los últimos 20 días
gamma = np.sum(gamma_obs[-20:])/20

# Graficamos los parámetros
beta_list = []
mu_star_list = []
mu_list = []
for t in linspace:
    beta_list.append(beta(t))
    mu_star_list.append(mu_star(t))
    mu_list.append(mu(t))

beta_list = np.array(beta_list)
mu_star_list = np.array(mu_star_list)

plt.figure(1)

plt.subplot(2, 2, 1, title='β')
plt.plot(beta_obs, color='b')
plt.plot(linspace, beta_list, color='r')

plt.subplot(2, 2, 2, title='μ*')
plt.plot(mu_star_obs, color='b')
plt.plot(linspace, mu_star_list, color='r')

plt.subplot(2, 2, 3, title='μ')
plt.plot(linspace, mu_list, color='k')

plt.subplot(2, 2, 4, title='γ')
plt.plot(gamma_obs, color='b')
plt.plot(gamma*np.ones(data_length), color='r')

# La tasa de nacimiento es un parámetro fijo
lambd = 1.9e-5


# derivada en un momento t para los valores y = (N, S, I, R)
def f(y, t, beta, lambd, mu, mu_star, gamma):
    N, S, I, R = y

    dN = N*lambd - mu_star(t)*I - N*mu(t)
    dS = -beta(t)*S*I + N*lambd - S*mu(t)
    dI = beta(t)*S*I - I*mu_star(t) - gamma*I - I*mu(t)
    dR = I*gamma - R*mu(t)

    return dN, dS, dI, dR


# Sacamos la solucion de la EDO
sol = odeint(f, y_0, linspace, args=(beta, lambd, mu, mu_star, gamma))
N, S, I, R = np.split(sol, range(1, 4), axis=1)


def plot_SIR_curves(a, b):
    sub_linspace = linspace[a:b+1]

    # Graficamos las curvas N, S, I, R de nuestra simulación
    plt.figure(3)
    plt.title('Predicciones evolución epidemia')

    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['k', 'b'])
    plt.subplot(1, 2, 1, title='Susceptibles y totales')
    plt.plot(sub_linspace, sol[a:b+1, :2])
    plt.legend(['Totales', 'Susceptibles'], loc='lower right')

    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['r', 'g'])
    plt.subplot(1, 2, 2, title='Infectados y recuperados')
    plt.plot(sub_linspace, sol[a:b+1, 2:])
    plt.legend(['Infectados', 'Recuperados'], loc='lower right')


# Función para sacar los infectados en el momento t a partir de las curvas
def contagios(t):
    return beta(t)*S[t, 0]*I[t, 0]


# Función para sacar los fallecidos por covid en el momento t a partir de las curvas
def fallecidos_covid(t):
    return mu_star(t)*I[t, 0]


# Función para sacar los fallecidos totales en el momento t a partir de las curvas
def fallecidos_totales(t):
    return mu_star(t)*I[t, 0] + N[t, 0]*mu(t)


# Calculamos los parámetros diarios y los graficamos
def parametros_diarios(a, b, verbose=False):
    contagios_list = []
    fallecidos_covid_list = []
    fallecidos_total_list = []

    umbral_contagios = -1
    umbral_fallecimientos = -1

    sub_linspace = linspace[a:b+1]
    for t in sub_linspace:
        contagios_list.append(contagios(t))
        fallecidos_covid_list.append(fallecidos_covid(t))
        fallecidos_total_list.append(fallecidos_totales(t))

        if verbose:
            if umbral_contagios < 0 and contagios_list[-1] >= 5000:
                umbral_contagios = t
                print('A partir del día {} se superan los 5000 contagios diarios el día {}'.format(a, umbral_contagios))

            if umbral_fallecimientos < 0 and fallecidos_covid_list[-1] >= 50:
                umbral_fallecimientos = t
                print('A partir del día {} se superan los 50 fallecimientos diarios el día {}'.format(a, umbral_fallecimientos))

    # Los pasamos a numpy
    contagios_list = np.array(contagios_list)
    fallecidos_covid_list = np.array(fallecidos_covid_list)
    fallecidos_total_list = np.array(fallecidos_total_list)

    if verbose:
        plt.figure(4)

        plt.subplot(1, 2, 1, title='Contagios diarios')
        plt.plot(sub_linspace, contagios_list, color='r')

        plt.subplot(1, 2, 2, title='Fallecimientos diarios')
        plt.plot(sub_linspace, fallecidos_covid_list, color='m')
        plt.plot(sub_linspace, fallecidos_total_list, color='k')
        plt.legend(['Fallecidos por covid', 'Fallecidos totales'], loc='upper right')

        step = linspace[1] - linspace[0]

        print('Contagios en [{}, {}]: '.format(a, b), np.sum(contagios_list)*step)
        print('Fallecidos por covid en [{}, {}]: '.format(a, b), np.sum(fallecidos_covid_list)*step)
        print('Fallecidos totales en [{}, {}]: '.format(a, b), np.sum(fallecidos_total_list)*step)
        print('tasa de morición en [{}, {}]: '.format(a, b),
              np.sum(fallecidos_covid_list) / np.sum(contagios_list))
        print('Variación total de población en [{}, {}]: '.format(a, b),
              N[b, 0] - N[a, 0])

    return contagios_list, fallecidos_covid_list, fallecidos_total_list


# Función para sacar la incidencia acumulada de cualquier curva los últimos 14 días
def incidencia_acumulada(a, b):
    assert(a >= 13)
    contagios_list, _, _ = parametros_diarios(0, b)

    contagios_cum = []
    sub_linspace = linspace[a:b + 1]
    umbral = 'primero'
    for t in sub_linspace:
        contagios_cum.append(np.sum(contagios_list[t - 13:t + 1])/N[t, 0]*100000)

        if umbral == 'primero' and contagios_cum[-1] >= 50:
            print('Se supera la incidencia acumulada de 50 el día {}'.format(t))
            umbral = 'segundo'

        elif umbral == 'segundo' and contagios_cum[-1] >= 100:
            print('Se supera la incidencia acumulada de 100 el día {}'.format(t))
            umbral = 'acabado'

    plt.figure(5)

    plt.plot(sub_linspace, contagios_cum, color='r')
    plt.title('Incidencia acumulada de contagios')