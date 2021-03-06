
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lns
import csv
import subprocess
import matplotlib.cm as cmx
import sys

def place_angle(z,Rd):

	theta = np.arctan(z/Rd)
	perm = [[1, 1],[1, -1],[-1, 1],[-1, -1],[0, -1],[-1, 0],[0,1],[1,0]]
	top_sign = int(np.sign(z))
	bottom_sign = int(np.sign(Rd))
	perm_i = [top_sign,bottom_sign]

	if perm_i == perm[0]:
		theta = theta

	elif perm_i == perm[1]:
		theta = np.pi + theta

	elif perm_i == perm[2]:
		theta = 2*np.pi + theta

	elif perm_i == perm[3]:
		theta = np.pi

	elif perm_i == perm[4]:
		theta = np.pi

	elif perm_i == perm[5]:
		theta = 3*np.pi/2

	elif perm_i == perm[6]:
		theta = 0

	elif perm_i == perm[7]:
		theta = np.pi/2
	else:
		raise Exception("Math Error")

	return theta


def cart_to_tor(filename, R0, save_dir, main_dir):

	#os.chdir(main_dir)

	lineas = int(subprocess.check_output(['wc', '-l', filename]).split()[0])-2
	print('{}'.format(lineas))
	# part_path = np.loadtxt(filename, dtype=float,skiprows=1,delimiter='\t')

	part_path = np.loadtxt(filename, dtype=float, skiprows=1, delimiter=',', max_rows=lineas)

	#Cilindrical Toroidal Transform
	shp = part_path.shape
	#tor_coords = np.zeros(shp) #7 para los nuevos valores
	rotations = 0 
	os.chdir(save_dir)
	with open('tor'+filename,'w') as torfile:
		torfile.write("r,phi,theta,|B|, Br, Bphi, Bz\n")
		for i in range(shp[0]):
			xcart = part_path[i, 0]
			ycart = part_path[i, 1]
			zcart = part_path[i, 2]
			Bvalue = part_path[i, 3]
			Bx = part_path[i, 4]
			By = part_path[i, 5]
			Bz =part_path[i, 6]

			R = np.sqrt(xcart**2 + ycart**2)
			Rdiff = R - R0

			theta_tor = np.arctan2(zcart, Rdiff)+np.pi
			if theta_tor < 0: 
				theta_tor = theta_tor+2*np.pi
			
			phi_tor = np.arctan2(ycart, xcart)
			if phi_tor < 0: 
				phi_tor = phi_tor+2*np.pi
			
			
			if i != 0:
				xviejo = part_path[i-1,0]
				yviejo = part_path[i-1,1]
				posviejo = list(np.sign([xviejo, yviejo]))
				poscart = list(np.sign([xcart, ycart]))
				if posviejo == [1, -1] and poscart == [1, 1]: 
					rotations += 1 
					#print(rotations)
			
		
			#Sphi_tor = phi_tor + rotations*2*np.pi
			if zcart != 0: 
				r_tor = zcart/np.sin(theta_tor)
			else: 
				r_tor = abs(Rdiff)
				
			Br = Bx*np.cos(phi_tor)+By*np.sin(phi_tor)

			Bphi = -Bx*np.sin(phi_tor)+By*np.cos(phi_tor)

			#tor_coords[i] = r_tor, phi_tor, theta_tor, Bvalue, Br, Bphi, Bz
			torfile.write('{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n'.format(r_tor, phi_tor, theta_tor, Bvalue, Br, Bphi, Bz))
	
		torfile.close()
	os.chdir(main_dir)

	#return tor_coords


def tor_to_cart(filename,R0):


	part_path = np.loadtxt(filename, dtype=float,skiprows=1,delimiter='\t')

	#Transformar a coordenadas cartesianas
	shp = part_path.shape
	cart_coords = np.zeros(shp)

	with open('cart'+filename,'w') as cartfile:
		cartfile.write("x,y,z,|B|,Bx,By,Bz\n")
		for i in range(shp[0]):

			rtor = part_path[i, 0]
			phitor = part_path[i, 1]
			thetator = part_path[i, 2]
			Bvalue = part_path[i, 3]
			Bx = part_path[i, 4]
			By = part_path[i, 5]
			Bz =part_path[i, 6]

			x = (R0+rtor*np.cos(thetator))*np.sin(phitor)
			y = (R0+rtor*np.cos(thetator))*np.cos(phitor)
			z = rtor*np.sin(thetator)

			cart_coords[i] = x, y, z, Bvalue, Bx, By, Bz
			cartfile.write('{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n'.format(x, y, z, Bvalue, Bx, By, Bz))

		cartfile.close()

	#return cart_coords


def directorio_cart_a_tor(work_dir, save_dir):
	"""
	Entradas: save_dir -- str; directorio para guardar el archivo.
			, work_dir -- str; directorio donde se encuentran todos los archivos de las   
			trayectorias 
	Salidas: directorio con archivos en coordenadas toroidales. 
	Descrip: transfoma de coordenadas cartesianas a coordenads toroidales; se guardan en "save_dir".
	"""
	print('Transformando archivos en {}'.format(work_dir))
	lista_archivos = os.listdir(work_dir)
	print(lista_archivos)
	lista_archivos.sort()
	print(lista_archivos)
	R_0 = 0.2477
	main_dir = os.getcwd()
	for iArch in lista_archivos:
		os.chdir(work_dir)
		print('\t-> Cambiando coordenadas de archivo: {}'.format(iArch))
		cart_to_tor(iArch, R_0, save_dir, main_dir)

def errores_transf(file1, file2):
	
	part_path1 = np.loadtxt(file1, dtype=float,skiprows=2,delimiter=',')
	part_path2 = np.loadtxt(file2, dtype=float,skiprows=2,delimiter='\t')

	diferencia = (part_path1 - part_path2)/part_path1
	error = np.max(abs(diferencia))


	print(error, np.where(abs(diferencia) == error))

def calc_angle_grafico(y,z):

	R0 = 0.2477
	R = y - R0
	theta = place_angle(z,R)
	rad = np.sqrt(R**2+z**2)
	return theta, rad

def indv_plot_polvstor(filename,main_dir,save_dir):

	lineas = int(subprocess.check_output(['wc', '-l', filename]).split()[0])-2
	print('\t\tNumero de lineas de archivo '+filename+': {}'.format(lineas))
	data = np.loadtxt(filename, dtype=float, skiprows=1, delimiter='\t', max_rows=lineas)
	shp = data.shape
	_, phi, theta, B, _, _, _ = data.T

	for i in range(shp[0]):
		if theta[i] > np.pi : 
			theta[i] = theta[i]-2*np.pi


	theta = theta*180/np.pi
	phi = phi*180/np.pi
	os.chdir(save_dir)
	plt.scatter(phi[1400:3100],theta[1400:3100], marker=".", s=1, c="black",zorder=100)
	plt.scatter(phi, theta, marker=".", c=B, cmap=cmx.jet, zorder=0)	
	plt.title("Magnitud de B vs posicion angular "+filename)
	plt.xlabel(r"angulo toroidal $\phi$")
	plt.ylabel(r"angulo poloidal $\theta$")
	plt.xticks(np.linspace(0,360,7))
	plt.yticks(np.linspace(-180,180,7))
	cbar = plt.colorbar()
	cbar.set_label("Magnitud de campo magnetico B")
	plt.savefig("single"+str(filename)+".png")
	plt.close()
	os.chdir(main_dir)

def plot_polvstor(work_dir,save_dir):

	print('Graficando archivos en {}'.format(work_dir))
	lista_archivos = os.listdir(work_dir)	#Lista de archivos
	lista_archivos.sort()
	main_dir = os.getcwd()	#Directorio principal
	print("\tAccesando al directorio: "+work_dir+"\n")
	
	for iArch in lista_archivos:
		os.chdir(work_dir)
		print("\tAnalizando archivo "+iArch+"...")
		indv_plot_polvstor(iArch, main_dir, save_dir)

	os.chdir(main_dir)
	print("\tGraficos en"+work_dir+"listos.")
	
#------------------------------ECRH surface mapping ------------------------------#
from cycler import cycler

def Plasma_Poincare_map_o(file_dir, angle_cu):
	

	print('Graficando archivos en {}'.format(file_dir))

	main_dir = os.getcwd() 					
	#Se guarda el path del directorio inicial puesto que se necesitan accesar a los directorios deseados. 
	
	list_archiv = os.listdir(file_dir)
	list_archiv.sort()

	angle_cut = angle_cu*np.pi/180
	dphi = 0.01
	phi_up = angle_cut + dphi/2
	phi_down = angle_cut - dphi/2
	R0 = 0.2477 

		
	fig, ax = plt.subplots(figsize=(10,10), dpi=400)
	custom_cycler = cycler(color=['c', 'm', 'y', 'k'])
	ax.set_prop_cycle(custom_cycler)

	plotcoords = []
	with open('RZ'+str(angle_cu)+'.txt', 'w') as filea:
		
		filea.write('R\tZ\t@'+str(angle_cu))
		os.chdir(file_dir)
		for iArch in list_archiv: #Se itera sobre los archivos en el directorio
			print(iArch)
			data = np.loadtxt(iArch, dtype=float, skiprows=2, delimiter='\t')
			
			#Cilindrical Toroidal Transform
			iRange, _ = data.shape
		
			for i in range(iRange):

				if data[i,1] > phi_down and data[i,1] < phi_up:
					rtor = data[i,0]
					thetator = data[i,2]
					R = (R0+rtor*np.cos(thetator))
					z = rtor*np.sin(thetator)
					plotcoords.append([R,z])
			
			plotcoord = np.array(plotcoords)
			if plotcoord.shape[0]==0:
				continue
			else:
				radius, z = plotcoord.T
				for i in range(len(radius)):
					filea.write('{}\t{}\n'.format(radius[i], z[i]))
				filea.write(iArch+'\n')
				ax.scatter(radius, z, marker='.', s=1.3)

		filea.close()
	
	os.chdir(main_dir) 

	r = 0.18/2
	theta = np.linspace(0, 2*np.pi, 60)
	x = r*np.cos(theta) + R0
	y = r*np.sin(theta)
	#plt.scatter(R0,0, marker='.')
	ax.plot(x,y)
	ax.scatter(0.2477, 0, color='black', s=0.002)
	ax.set_xlabel('R [m]')
	ax.set_ylabel('Z [m]')
	fig.gca().set_aspect('equal', adjustable='box')
	ax.set_title('Diagrama de Poincare\n$\\theta = {:.2f}$'.format(angle_cut*180/np.pi))
	fig.savefig("poincare"+str(angle_cu)+".png")

	return fig, ax, phi_up, phi_down


def ECRHsurfmap(work_dir):

	R0 = 0.2477					#radio caracteristico SCR-1(magnetico)
	heatfreq = 2.45*10**9 		#frecuencia de calentamiento
	B2 = (2*np.pi*heatfreq*9.10938291*10**-31)/((1.60217*10**-19)*2)

	surfdata2 = list() 

	print('Analizando archivos en {}'.format(work_dir))
	lista_archivos = os.listdir(work_dir)
	lista_archivos.sort()
	main_dir = os.getcwd()
	print("\tAccesando al directorio: "+work_dir+"\n")

	os.chdir(work_dir)

	#Se compara para cada archivo si su valor de magnitud es igual a B2
	for iArch in lista_archivos:
		
		print("\tAnalizando archivo "+iArch+"...")
		
		#Analisis de los datos

		lineas = int(subprocess.check_output(['wc', '-l', iArch]).split()[0])-2
		print('\t Cantidad de lineas: {}'.format(lineas))
		data = np.loadtxt(iArch, dtype=float, skiprows=1, delimiter='\t', max_rows=lineas)
		shp = data.shape
		Bintwidth = 0.0002
		for i in range(shp[0]):
			if B2-Bintwidth/2 < data[i,3] and B2+Bintwidth/2 > data[i,3]:
				surfdata2.append(data[i])
		
	
	surf2 = np.array(surfdata2)
	os.chdir(main_dir)	 

	#Generar grafico poincare
	angle = 0#float(input("Ingresar angulo para la seccion: "))

	#Se grafica primero un mapa de poincare en el angulo dado para mostrar la superficie en el angulo dado
	fig, ax, phi_up, phi_down = Plasma_Poincare_map_o(work_dir, angle)

	iR2, _ = surf2.shape

	coords2 = list()

	for i in range(iR2):
		if phi_down < surf2[i,1] and surf2[i,1] < phi_up:
			R = R0+surf2[i,0]*np.cos(surf2[i,2])
			z = surf2[i,0]*np.sin(surf2[i,2])
			coords2.append([R,z])

	coords2 = np.array(coords2)
	
	R2, Z2 = coords2.T

	fit = np.polyfit(Z2, R2, 3)		

	#dependiendo de la orientacion es mejor intercambiar Z por R en np.polyfit, de esta manera np.polyfit(R2,Z2,3)

	points = np.linspace(min(Z2),max(Z2))
	fitpoly = np.polyval(fit ,points)

	#ax.scatter(R2, Z2, marker=".", color="blue", label="2do harmonico"), graficacion de puntos.  
	ax.plot(fitpoly, points, color="red")
	plt.savefig("poincareEH"+str(angle)+".png", dpi=300)
	plt.close()


def Calculo_iota_4(archivo_trayectoria, angulo):

	from math import isclose

	datos = np.loadtxt(archivo_trayectoria, skiprows = 2)
	R0 = 0.2477
	nPasos, _ = datos.shape
	punto_inicial = np.zeros(3)
	punto_final = np.zeros(3)
	angulo = angulo*np.pi/180

	nRev = 1

	# Variable a la que suman los delta_phi
	sumatoria_delta_phi = 0

	print('\n\nProcesando archivo {}'.format(archivo_trayectoria))

	rtor, _, _, _, _, _, _ = datos.T 
	R = np.max(rtor)+R0 


	for iPaso in range(nPasos):
		# si estamos en el primer paso se toma como punto inicial el primer punto
		if nRev < 200:
			if iPaso == 0:
				punto_inicial = datos[iPaso, :]


			elif datos[iPaso, 1] > angulo+2*np.pi*nRev:
				
				punto_final = datos[iPaso, :]

				# Se aumenta el contador
				nRev += 1

				delta_phi = abs(punto_final[2] - punto_inicial[2])

				sumatoria_delta_phi +=delta_phi

				punto_inicial = np.array(punto_final)
				
			else:
				continue


	iota = sumatoria_delta_phi/(nRev*2*np.pi)
	print('\tsumatoria: {:.7f}\tnRev: {}\tiota: {:.7f}\tlineas: {}\n\tRadio promedio: {}'.format(sumatoria_delta_phi, nRev, iota, nPasos, R))


	return iota, R

def grafico_iota(directorio, angulo): 

	print('Analizando archivos en {}'.format(directorio))
	lista_archivos = os.listdir(directorio)
	lista_archivos.sort()
	main_dir = os.getcwd()
	print("\tAccesando al directorio: "+directorio+"\n")

	iotas = list()
	radii = list()
	os.chdir(directorio)	
	for iArch in lista_archivos: 
		iota, rad= Calculo_iota_4(iArch, angulo)
		iotas.append(iota)
		radii.append(rad)
	
	vup = np.array([radii,iotas]).T
	print(vup)
	vup = vup[np.argsort(vup[:,0])]
	radii, iotas = vup.T
	os.chdir(main_dir)
	plt.plot(radii[0:len(radii)-1],iotas[0:len(iotas)-1], color="black") 
	plt.xlabel("R [m]")
	plt.ylabel("iota")
	plt.title(r"Transformada rotacional vs posicion radial para $\theta = {}$".format(angulo))
	plt.savefig("iots.png",dpi=300)
	plt.close()


	


