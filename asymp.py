import numpy as np
from scipy.integrate import odeint

from matplotlib import pyplot as plt
import matplotlib as mpl


class ModeloAsintomaticos():
    def __init__(self, N_0, S_0, I_0, A_0, R_0, lambd, mu, mu_star, gamma,
                 gamma_star, beta_1, beta_2, beta_3, beta_4):
        assert(N_0 == S_0 + I_0 + A_0 + R_0)

        linspace = range(61)
        y_0 = (N_0, S_0, I_0, A_0, R_0)

        sol = odeint(self.f, y_0, linspace, args=(lambd, mu, mu_star, gamma,
                 gamma_star, beta_1, beta_2, beta_3, beta_4))

        N, S, I, A, R = np.split(sol, range(1, 5), axis=1)
        N = N[:, 0]
        S = S[:, 0]
        I = I[:, 0]
        A = A[:, 0]
        R = R[:, 0]


        plt.figure(1)
        plt.title('Predicciones evolución epidemia')

        plt.subplot(2, 2, 1, title='Totales')
        plt.plot(linspace, sol[:, 0], color='k')
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)

        plt.subplot(2, 2, 2, title='Susceptibles')
        plt.plot(linspace, sol[:, 0], color='b')
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)

        mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['r', 'y'])
        plt.subplot(2, 2, 3, title='Infectados y asintomáticos')
        plt.plot(linspace, sol[:, 2:4])
        plt.legend(['Infectados', 'Asintomáticos'], loc='upper right')

        plt.subplot(2, 2, 4, title='Recuperados')
        plt.plot(linspace, sol[:, 4], color='g')

        fallecidos = mu_star*I
        contagiados_por_infectado = beta_1*np.multiply(S, I)
        contagiados_por_asintomatico = beta_3*np.multiply(S, A)

        plt.figure(2)
        plt.title('Curvas diarias')

        plt.subplot(1, 2, 1, title='Fallecimientos diarios')
        plt.plot(fallecidos, color='k')

        plt.subplot(1, 2, 2, title='Contagios diarios')
        plt.plot(contagiados_por_infectado, color='r')
        plt.plot(contagiados_por_asintomatico, color='y')
        plt.plot(contagiados_por_infectado + contagiados_por_asintomatico, color='b')
        plt.legend(['Contagios por infectado', 'Contagios por asintomático', 'Contagios totales'])

        fallecidos_totales = np.sum(fallecidos)
        contagios_totales = np.sum(contagiados_por_infectado + contagiados_por_asintomatico)
        ratio_contagios_asintomaticos = np.sum(contagiados_por_asintomatico)/contagios_totales

        print('Fallecidos en los siguiente 60 días:', fallecidos_totales)
        print('Contagios en los siguiente 60 días:', contagios_totales)

        print('Porcentaje de contagios por asintomático:', ratio_contagios_asintomaticos*100)
        print('Fallecidos por asintomático:', ratio_contagios_asintomaticos*fallecidos_totales)

    def f(self, y, t, lambd, mu, mu_star, gamma,
                 gamma_star, beta_1, beta_2, beta_3, beta_4):
        N, S, I, A, R = y
        assert (abs(N - (S + I + A + R)) < 1e-3)

        dN = N*lambd - mu_star*I - N*mu
        dS = - (beta_1 + beta_2)*S*I - (beta_3 + beta_4)*S*A - S*mu + N*lambd
        dI = beta_1*S*I + beta_3*S*A - I*mu_star - gamma*I - mu*I
        dA = beta_2*S*I + beta_4*S*A - gamma_star*A - mu*A
        dR = I*gamma + A*gamma_star - R*mu
        assert (abs(dN - (dS + dI + dA + dR)) < 1e-3)

        return dN, dS, dI, dA, dR