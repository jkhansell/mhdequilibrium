import numpy as np
import plotly.graph_objects as go

class plot_switch_coils(object):

    def __init__(self, coildata):
        #coildata = read_coils("coils.txt") --- read read_coils.py
        self.coildata = coildata
        
    def plot_type(self, plottype): 
        """
        Dispatch method: 

        PLOT_COILS(coildata[,'simple'/'plane'/'other']) Plots coils
        PLOT_COILS plots the data read by READ_COILS from a VMEC XGRID file.

        Usage:
        coil_data=read_coils('coils.test');
        plot_coils(coil_data,'default');  % Plots each fillament
        plot_coils(coil_data,'plane');  %Plots intersection of coils (phi=0)
        plot_coils(coil_data,'simple'); %Plots volumetric rendering of groups
        plot_coils(coil_data,'field_period'); %Plots only for field period
        plot_coils(coil_data,'tube'); %Plots coils as volumetric tubes.
        plot_coils(coil_data,'tube','nedges',20,'tubewidth',0.2); %Plots coils 
            as volumetric tubes of width tubewidth and having nedges edges.
            Defaults are 4 edges (square) and tubewidth 0.15.

        """
        method_name = str(plottype)
        method = getattr(self, method_name, lambda: 'Invalid plot type')
        return method()
  
    def default(self):
        print("\tPlotting data...")
        
        layout = go.Layout(
            scene=dict(
                aspectmode='data'
            ))
        fig = go.Figure(layout=layout)

        for key in self.coildata.keys():
            
            if key != 'periods': 
                    x, y, z, _ = self.coildata[key]
                    fig.add_trace(go.Scatter3d(x=x, y=y, z=z))
                    
        print("\tDone!")
        fig.show()



def plot_coils(data, plottype): 
    plot = plot_switch_coils(data)
    plot.plot_type(plottype)