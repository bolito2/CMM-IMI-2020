from asymp import *

# Valores iniciales
S_0 = 4.2e7
I_0 = 6e4
A_0 = 2.4e5
R_0 = 5e6
N_0 = S_0 + I_0 + A_0 + R_0

lambd = 2.1e-5
mu = 2.4e-5
mu_star = 5e-4
gamma = 1e-2
gamma_star = 7e-2
beta_1 = 2e-2/N_0
beta_2 = 8e-2/N_0
beta_3 = 1e-2/N_0
beta_4 = 4e-2/N_0

# II
print('ej 4 parte II\n')
model = ModeloAsintomaticos(N_0, S_0, I_0, A_0, R_0, lambd, mu, mu_star, gamma,
                 gamma_star, beta_1, beta_2, beta_3, beta_4)

plt.show()

# III
print('---------------------\n\n')
print('ej 4 parte III\n')

beta_1 /= 2
beta_2 /= 2
beta_3 /= 10
beta_4 /= 10

model2 = ModeloAsintomaticos(N_0, S_0, I_0, A_0, R_0, lambd, mu, mu_star, gamma,
                 gamma_star, beta_1, beta_2, beta_3, beta_4)

plt.show()