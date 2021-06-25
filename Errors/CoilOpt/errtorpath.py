"""

Rutina de código que discretiza el proceso de cálculo de error 
local y promedio.


Elaborado por: Bryam Nuñez
Ing. Física, TEC

"""


import numpy as np

def DatosCurvasOrdenados(datos):
    """"
    Función que recibe los datos y los ordena tal que  0<phi<2pi en orden creciente
    y 0<theta<2pi en orden creciente.
    """

    datos_vuelta = {0: []}
    # primer ciclo
    n = 0
    for i in range(1, len(datos)):
        if datos[i-1, 1] < datos[i, 1]:
            # se apendan los datos que cumplen phi<2pi
            datos_vuelta[n].append(datos[i-1])
        else:
            # si se pasa de vuelta se apenda phi<2pi y se crea una nueva entrada
            # en el diccionario para guardar datos correspondientes a otra vuelta
            datos_vuelta[n].append(datos[i-1])
            n += 1
            datos_vuelta[n] = []

    # segundo ciclo donde se ordena tal que 0<theta<2pi en orden creciente
    lista = []
    for i in range(len(datos_vuelta.keys())):
        # se apenda a la lista los valores de theta
        # y el número de llave del diccionario
        lista.append([datos_vuelta[i][0][2], i])
    # se ordenan de menor a mayor los valores de theta
    lista.sort(key=lambda x: x[0], reverse=False)
    # se genera una lista donde se guardan todos los datos ordenados
    datos_f = {}
    for i in range(len(datos_vuelta.keys())):
        datos_f[i] = datos_vuelta[lista[i][1]]

    return datos_f

def CurvasCartesianas(dnao):
    """
    Función que recibe los datos ordenados y los pasa a coordenadas cartesianas
    """
    R0 = 0.2477
    carte = {}
    for i in range(len(dnao.keys())):
        carte[i] = []
        for j in range(len(dnao[i])):
            Bmag = dnao[i][j][3]  # magnitud campo magnético
            Bx = dnao[i][j][4]  # componente en x del campo magnético
            By = dnao[i][j][5]  # componente en y del campo magnético
            Bz = dnao[i][j][6]  # componente en z del campo magnético
            R = R0+dnao[i][j][0]*np.cos(dnao[i][j][2])
            x = R*np.cos(dnao[i][j][1])  # componente en x del la partícula
            y = R*np.sin(dnao[i][j][1])  # componente en y del la partícula
            # componente en z del la partícula
            z = dnao[i][j][0]*np.sin(dnao[i][j][2])
            # componentes en cartesianas
            cartesianas = [x, y, z, Bmag, Bx, By, Bz]
            carte[i].append(cartesianas)
        carte[i] = np.array(carte[i])
    return carte

def DatosFinales(datos):
    """
    Función que recive el path de los datos y los ordena y transforma 
    a coordenadas cartesianas
    """
    dnao = DatosCurvasOrdenados(datos)
    datos = CurvasCartesianas(dnao)
    l = len(datos)
    # se agrega una última curva idéntica a la primera,para poder iterar sobre
    # la última curva utilizando la primera como referencia.
    datos[l] = datos[0]
    return datos

def ErrorLocal(dnao):

    l = len(dnao.keys())-1  # curvas a iterar para el error
    # se omite la última porque es idéntica a la primera. Esta
    # última curva sirve para para generar los vectores unitarios de la penúltima
    # curva que en realidad es la última curva que da el BS-Soltrac

    error = {}  # diccionario con error locales por curva

    # Ciclo principal donde se calculan los errores locales de todas las curvas.
    for i in range(0, l):

        error[i] = []  # se analiza la i-ésima curva
        # longitud de la i+1-ésima curva (con la que se compara)
        m = len(dnao[i+1])-1
        # longitud de la i-ésima curva(donde se calcula el error)
        n = len(dnao[i])-1

        # se itera sobre cada punto de la curva i-ésima.
        for j in range(n):
            # condición cuando se tienen una relación 1:1 de puntos entre curvas
            if j <= m:
                dx1 = dnao[i][j+1][0]-dnao[i][j][0]  # dx curva i-ésima
                dy1 = dnao[i][j+1][1]-dnao[i][j][1]  # dy curva i-ésima
                dz1 = dnao[i][j+1][2]-dnao[i][j][2]  # dz curva i-ésima
                r1_j = np.array([dx1, dy1, dz1])  # r curva i-ésima

                dx2 = dnao[i+1][j][0]-dnao[i][j][0]  # dx curva i+1-ésima
                dy2 = dnao[i+1][j][1]-dnao[i][j][1]  # dy curva i+1-ésima
                dz2 = dnao[i+1][j][2]-dnao[i][j][2]  # dz curva i+1-ésima
                r2_j = np.array([dx2, dy2, dz2])  # r curva i+1-ésima

                # producto cruz entre r_i-ésima y r_i+1-ésima
                N = np.cross(r1_j, r2_j)
                dA = np.sqrt(np.dot(N, N))  # área paralelogramo entre los r's
                n = N/dA  # vector unitario en dirección perpendicular a las dos curvas

                Bmag = dnao[i][j][3]  # magnitud campo magnético
                # componente en x del campo magnético en curva i
                Bx_j = dnao[i][j][4]
                # componente en y del campo magnético en curva i
                By_j = dnao[i][j][5]
                # componente en z del campo magnético en curva i
                Bz_j = dnao[i][j][6]
                # campo magnético en curva i-ésima
                B = np.array([Bx_j, By_j, Bz_j])

                # error local en el punto de anális
                e_j = abs(np.dot(B, n))/Bmag
                error[i].append([e_j, dA])

            # caso donde la curva i+1 tiene menos puntos que la curva i
            else:
                dx1 = dnao[i][j+1][0]-dnao[i][j][0]  # dx curva i-ésima
                dy1 = dnao[i][j+1][1]-dnao[i][j][1]  # dy curva i-ésima
                dz1 = dnao[i][j+1][2]-dnao[i][j][2]  # dz curva i-ésima
                r1_j = np.array([dx1, dy1, dz1])  # r curva i-ésima

                # se utiliza el último punto de la curva i+1
                dx2 = dnao[i+1][m][0]-dnao[i][j][0]  # dx curva i+1-ésima
                dy2 = dnao[i+1][m][1]-dnao[i][j][1]  # dy curva i+1-ésima
                dz2 = dnao[i+1][m][2]-dnao[i][j][2]  # dz curva i+1-ésima
                r2_j = np.array([dx2, dy2, dz2])  # r curva i+1-ésima

                # producto cruz entre r_i-ésima y r_i+1-ésima
                N = np.cross(r1_j, r2_j)
                dA = np.sqrt(np.dot(N, N))  # área paralelogramo entre los r's
                n = N/dA  # vector unitario en dirección perpendicular a las dos curvas

                Bmag = dnao[i][j][3]  # magnitud campo magnético
                # componente en x del campo magnético en curva i
                Bx_j = dnao[i][j][4]
                # componente en y del campo magnético en curva i
                By_j = dnao[i][j][5]
                # componente en j del campo magnético en curva i
                Bz_j = dnao[i][j][6]
                # campo magnético en curva i-ésima
                B = np.array([Bx_j, By_j, Bz_j])

                # error local en el punto de anális
                e_j = np.nan_to_num(abs(np.dot(B, n))/Bmag)
                error[i].append([e_j, dA])

    return error

def ErrorSup(dnao, area):
    """
    Función que da el error promedio en una superficie que cumpe JxB=nabla_p
    """
    # primero se calcula la cantidad de puntos en la superficie
    cantidad_puntos = 0
    for i in range(len(dnao.keys())):
        cantidad_puntos += len(dnao[i])
    A = area  # área envolvente
    n = cantidad_puntos -1   # número de puntos en la superficie
    dA = A/n  # diferencial de área
    error_dA = 0  # suma del error por diferencial de área
    for i in range(len(dnao.keys())):
        for j in range(len(dnao[i])):
            # se suman todas las contribuciones al error
            error_dA += dnao[i][j][0]*dA
    error_dA_sup = error_dA/A  # error promedio de la superficie
    return error_dA_sup