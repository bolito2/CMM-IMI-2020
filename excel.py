import openpyxl

document = openpyxl.load_workbook('CMM_III_2020_Datos_Definitivo.xlsx', data_only=True)
sheet = document['Sheet1']


def get_initial_data():
    i_start = 14
    length = 4
    y_0 = []

    for i in range(i_start, i_start + length):
        y_0.append(sheet['H{}'.format(i)].value)

    return y_0


def get_datos_diarios(col):
    i_start = 2
    i_end = 44
    casos_diarios = []

    for i in range(i_start, i_end):
        casos_diarios.append(sheet['{}{}'.format(col, i)].value)

    return casos_diarios
