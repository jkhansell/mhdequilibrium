"""
Rutina de código para leer archivos de texto de filamentos de bobinas en el siguiente formato: 

x   y   z   current
####################

Ejemplo:

periods 2
begin filament
mirror NULL
3.412616e-01   9.144077e-02   0.000000e+00   4.350000e+03
3.415617e-01   9.024450e-02   2.010517e-03   4.350000e+03
.
.
.
2.331786e-01   2.654216e-01   0.000000e+00   0.000000e+00  1 modular
3.412616e-01   9.144077e-02   0.000000e+00   4.350000e+03
3.415617e-01   9.024450e-02   2.010517e-03   4.350000e+03
.
.
.
2.331786e-01   2.654216e-01   0.000000e+00   0.000000e+00  1 modular
end


Elaborado por: Johansell Villalobos
Ing. Física, TEC

"""
#Importación de librerías

import numpy as np

##################################

def read_coil(filename):
    """
    Función que devuelve un arreglo de tipo np.array con las posiciones y corrientes de un set de bobinas 
    dado un nombre o una ruta de archivo. 

    filename :: nombre o ruta de archivo .txt 

    """
          
    with open(filename, mode='r') as file:
        
        fileread = file.readlines()
        coils = []
        coil = []
        count = 0
        for iLine in range(len(fileread)):
            if iLine == 0: 
                periods = int(fileread[iLine][-2])
                continue
            elif iLine == 1 or iLine == 2:
                continue

            else:
                line = fileread[iLine].split()
                if len(line) == 6:
                    count += 1
                    coil.append([float(elem) for elem in line[0:4]])
                    coils.append(np.array(coil))
                    coil = []

                elif line[0] == 'end' or line[0] =='END':
                    continue
 
                else: 
                    coil.append([float(elem) for elem in line])

    

    return np.array(coils)
                

      
                
        


