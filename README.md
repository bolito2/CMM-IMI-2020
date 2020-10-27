# CMM-IMI-2020

Código usado para el [Concurso de Modelización Matemática (CMM-IMI). Edición 2020](http://blogs.mat.ucm.es/cmm/edicion-2020/)

[Nuestra participación](https://drive.google.com/file/d/14REuMg0-EDMIOYZ1lwpVq9YgJST6Yrfq/view?usp=sharing) ha sido creada por Óscar Álvarez Sánchez(grado en Matematicas UCM), Daniel Andrés Solís(grado en Ingeniería Matemática UCM) y Pablo Gómez Morales(grado en Matemáticas UCM) a partir de los resultados de este programa.

## Cómo usar
### Instalar dependencias
```bash
pip install matplotlib openpyxl scipy numpy
```
### Predicciones enunciado
Para ver los resultados de cada pregunta ejecutar en el terminal
```bash
python ej*.py
```
para el ejercicio deseado

### Predicciones a medida
Para obtener predicciones en otro periodo de tiempo se puede usar el modelo básico creando un script de la siguente forma 
```python
from main import *  # Calcular predicciones en el intervalo [0, 240]

# Graficar los parámetros intermedios(gamma, mu , etc) en el intervalo [x, y] ⊂ [0, 240]
graficar_parametros(x, y) 

# Graficar las curvas SIR en el intervalo [x, y] ⊂ [0, 240]
plot_SIR_curves(x, y) 

# Calcular los contagios/fallecimientos diarios en el intervalo [x, y] ⊂ [0, 240] 
# Si verbose=True se grafican las curvas
# Si contagios_diarios=True se representa el punto en el que los valores superan el umbral
parametros_diarios(x, y, verbose=True, contagios_diarios=True) 

# Graficar la incidencia acumulada en el intervalo [x, y] ⊂ [13, 240]
incidencia_acumulada(x, y)  

plt.show()  # Mostrar las gráficas obtenidas
```

Y para obtener predicciones con el modelo de asintomáticos
```python
from asymp import * # Importar modelo

# Introducir los valores iniciales
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

# Inicializar y ejecutar el modelo
model = ModeloAsintomaticos(N_0, S_0, I_0, A_0, R_0, lambd, mu, mu_star, gamma,
                 gamma_star, beta_1, beta_2, beta_3, beta_4)

plt.show()  # Mostrar los resultados
