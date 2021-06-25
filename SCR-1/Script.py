import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lns
import csv
import subprocess
import matplotlib.cm as cmx
import sys
import plotly.graph_objects as go


from cycler import cycler

def Plasma_Poincare_map(file_dir, angle_cu):
	"""
	Entradas:file_dir (string) directorio del archivo con los datos de las superficies 				magneticas del plasma.

				angle_cut (float) ángulo en grados para el corte de la sección de poincaré. 

	Salidas: fig (matplotlib.pyplot figure) figura de mapa de poincare generado
				ax (matplotlib.pyplot axis) axis de la figura de mapa de poincare

				phi_up (float) limite superior de la pequeña caja que contiene las intersecciones de la trayectoria con el plano dado por phi = angle_cut

				phi_down (float) limite inferior de la pequeña caja que contiene las intersecciones de la trayectoria con el plano dado por phi = angle_cut

	Descrip: Función que toma las trayetorias generadas por el módulo BS-SOLCTRA y grafica un mapa de poincaré asociado al ángulo deseado. 

	NOTA: Si se utiliza en Linux no hay problema, con Anaconda presenta un WinError2, filenotfound en la linea []
	"""

	print('Graficando archivos en {}'.format(file_dir))

	main_dir = os.getcwd() 					
	#Se guarda el path del directorio inicial puesto que se necesitan accesar a los directorios deseados. 

	list_archiv = os.listdir(file_dir)
	list_archiv.sort()

	angle_cut = angle_cu*np.pi/180
	dphi = 0.03
	phi_up = angle_cut + dphi/2
	phi_down = angle_cut - dphi/2
	R0 = 0.2477 #Radio característico del SCR-1


	fig, ax = plt.subplots(figsize=(10,10), dpi=400)
	plotcoords = []
	with open('RZ'+str(angle_cu)+'.txt', 'w') as filea:
		
		filea.write('R\tZ\t@'+str(angle_cu)+'º\n')
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
					R = R0+rtor*np.cos(thetator)
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
				ax.scatter(radius, z, s=0.9)


		filea.close()

	os.chdir(main_dir) 
			
	r = 0.18/2
	theta = np.linspace(0, 2*np.pi, 60)
	x = r*np.cos(theta) + R0
	y = r*np.sin(theta)
	#plt.scatter(R0,0, marker='.'
	ax.legend(fontsize='small')
	ax.plot(x,y)
	ax.scatter(0.2477, 0, color='black', s=0.002)
	ax.set_xlabel('R [m]')
	ax.set_ylabel('Z [m]')
	fig.gca().set_aspect('equal', adjustable='box')
	ax.set_title('Diagrama de Poincare\n$\\theta = {:.2f}$'.format(angle_cu))
	fig.savefig("poincare"+str(angle_cu)+".png")

def indv_plot_polvstor(filename,main_dir,save_dir):
	"""
	Entradas: filename -- str; Nombre del archivo a analizar.
			, main_dir -- str; Path del directorio donde se está ejecutando el código.
			, save_dir -- str; Path del directorio donde se guardan los graficos.
	Salidas:  single_______.png imagen del gráfico hecho a partir de los datos
				encontrados en filename
	Descripción: Toma los archivos con datos del directorio "work_dir" y los analiza para
					generar graficas de B vs theta vs phi, y a su vez guardar estos archivos en "save_dir".

	Ej: indv_plot_polvstor(filename,"./results_6945","./graficos")
	"""


	data = np.loadtxt(filename, dtype=float, skiprows=1, delimiter='\t')
	shp = data.shape
	_, phi, theta, B, _, _, _ = data.T

	for i in range(shp[0]):
		if theta[i] > np.pi :
			theta[i] = theta[i]-2*np.pi

	theta = theta*180/np.pi
	phi = phi*180/np.pi
	os.chdir(save_dir)
	plt.scatter(phi, theta, marker=".", c=B, cmap=cmx.jet)
	plt.title("Magnitud de B vs posicion angular "+filename)
	plt.xlabel(r"Ángulo toroidal $\phi$")
	plt.ylabel(r"Ángulo poloidal $\theta$")
	plt.xticks(np.linspace(0,360,7))
	plt.yticks(np.linspace(-180,180,7))
	cbar = plt.colorbar()
	cbar.set_label("Magnitud de campo magnético B")
	plt.savefig("single"+str(filename)+".png")
	plt.close()
	os.chdir(main_dir)

def graphs(filename):
	data = np.loadtxt(filename, dtype=float, skiprows=2, delimiter='\t')
	x,y,z,B,bx,by,bz = data.T
	layout = go.Layout(scene=dict(aspectmode='data'))
	fig = go.Figure(data = go.Scatter3d(x=x,y=y,z=z,marker=dict(
		size=0.1,
		color=B,                # set color to an array/list of desired values
		colorscale='Viridis',   # choose a colorscale
		opacity=0.8),
		line=dict(
		color='darkblue',
		width=1)),layout=layout)
	fig.show()