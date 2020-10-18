from main import *

graficar_parametros(0, 240)
# I
print('Ej 3 parte I\n')
parametros_diarios(150, 240, verbose=True)
plt.show()

# II
print('\n\nEj 3 parte II\n')
incidencia_acumulada(150, 240)
plt.show()

# III
print('\n\nEj 3 parte III')
parametros_diarios(180, 240, verbose=True)
plt.show()