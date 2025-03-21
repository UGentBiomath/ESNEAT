# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 09:27:09 2021

@author: Cristian Camilo Gómez Cortés

### Description:
Create the user interface to interact with the automatic calibration process

### The input parameters are the following:
    - network: (str) - File path to the .inp file of the current network
    - format_measurement: (str) - File path to the .xlsx file which contains all the measurement related information.
    - format_infopipe: (str) - File path to the .xlsx file which contains all the age and material related information for each pipe.
### The output parameters are the following:
    - The user interface created with the implementation of tkinter
"""
# Load required packages
from tkinter import *
from tkinter import messagebox
from tkinter.ttk import *
from tkinter import filedialog
import wntr
import matplotlib.pyplot as plt
import pandas as pnd
from PIL import ImageTk, Image
import datetime
import threading
import numpy as np
import shutil

from sensitivity import sensitivity
from calibration import calibration
from validation import warning
from data_resampler import data_resample


network = ""
net_name = "Net"
format_measurement = ""
format_infopipe = ""
temp = ""
Elements = pnd.DataFrame(columns=["Type","ID"])
Real_data = []
time_step = 15
total_duration = 24

materials = ["DCI","GCI","PEH","PVC","STE","FCE","CON"]

measurements = False
pipeinfo = False

calib_roughness=False
calib_minorloss=False
calib_demand=False
calib_leak=False
calib_valves_coeff=False
calib_valves_stat=False
calib_demand_pattern=False

generations = 2
population = 10
loops= 1
hidden_layers = 0
disc_mass = 10
disc_energy = 4
percent_closure_valves = 20
rang_valve_coeff = [0.2,10.0]
rang_minorloss = [0.1,2.0]
rang_leak = [0,1]
rang_demand = [0.5,1.5]

class Aplication:
    def __init__(self):
        # Configuración de la raíz

        global roughness
        global minorloss
        global demand
        global leak
        global valves_coeff
        global valves_stat
        global demand_pattern
        
        global valves_coeff_min
        global minorloss_min
        global leak_coef_min
        global dem_min
        global valves_coeff_max
        global minorloss_max
        global leak_coef_max
        global dem_max
        global percent_closed_valves
        
        global generetions_total
        global population_total
        global flow_unit_factor
        global pressure_unit_factor
        global time_config

        self.root = Tk()

        roughness=IntVar()
        minorloss=IntVar()
        demand=IntVar()
        leak=IntVar()
        valves_coeff=IntVar()
        valves_stat = IntVar()
        demand_pattern = IntVar()
        
        valves_coeff_min=StringVar()
        minorloss_min=StringVar()
        leak_coef_min=StringVar()
        dem_min=StringVar()
        
        valves_coeff_max=StringVar()
        minorloss_max=StringVar()
        leak_coef_max=StringVar()
        dem_max=StringVar()
        
        percent_closed_valves = IntVar() 
        
        generetions_total = StringVar()
        population_total = StringVar()
        flow_unit_factor = DoubleVar()
        pressure_unit_factor = DoubleVar()
        time_config = StringVar()

        self.root.iconbitmap('Logos/Icono.ico')
        self.root.title("EXACTLY - Automatic Calibration Tool for Water Distribution Networks")
        self.root.geometry("800x450")
        # Main menu
        menubar = Menu(self.root)
        self.root.config(menu=menubar)
        # File menu
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="Load file", command=self.load_network)
        submenuform=Menu(menubar, tearoff=0)
        submenuform.add_command(label="Create Pipe information", command=self.format_infopipe)
        submenuform.add_command(label="Create Measurements", command= self.format_measurements)
        submenuform.add_command(label="Load Pipe Information", command = self.load_infopipe)
        submenuform.add_command(label="Load Measurements", command = self.load_measurements)
        filemenu.add_cascade(label="Formats",menu=submenuform)
        
        # Validation menu
        valmenu = Menu(menubar, tearoff=0)
        valmenu.add_command(label="Validation", command = self.validate)
        valmenu.add_command(label="Data resample", command = self.resample_data_rang)
        
        # Calibration menu
        simmenu = Menu(menubar, tearoff=0)
        simmenu.add_command(label="Sensitivity Analysis", command = self.bar_init)
        simmenu.add_command(label="Calibration Process", command = self.calibracion)
        
        # Help menu
        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="User Manual")
        helpmenu.add_command(label="About ...")
        
        # Headers menu
        menubar.add_cascade(label="Input File", menu=filemenu)
        menubar.add_cascade(label="Validation", menu=valmenu)
        menubar.add_cascade(label="Simulation", menu=simmenu)
        menubar.add_cascade(label="Help", menu=helpmenu)
        # Main loop of the app
        self.panel = Label(self.root)
        self.root.mainloop()
        
    #Load network inp file, create a copy of it and assign Hydraulics Save file
    def load_network(self):
        global network
        global net_name
        global temp
        network = filedialog.askopenfilename(initialdir="/",title = "Open",
                                         filetype = (("Network files(.inp)",".inp"),("All files","*.*")))
        net_name = network.split("/")[-1].split(".")[-2]
        u = network.split("/")
        u[-1] = "tempex.inp"
        temp = "/".join(u)
        shutil.copy(network, temp)
        busqueda = "[COORDINATES]"
        with open(temp, "r") as archivo_lectura:
            data = archivo_lectura.readlines()
        num = 0
        for linea in data:
            linea = linea.rstrip()
            if busqueda in linea:
                x = temp.split(".")
                x[-1] = "hyd"
                y = ".".join(x)
                data[num-1] = ' Hydraulics Save       "'+y+'"\n\n'
            num += 1
        with open(temp, 'w') as file:
            file.writelines( data )
        self.wn = wntr.network.WaterNetworkModel(temp)
        wntr.graphics.plot_network(self.wn)
        plt.savefig(net_name,dpi=300)
        net_plot = ImageTk.PhotoImage(Image.open(net_name+".png").resize((800,450)))
        self.panel.configure(image=net_plot)
        self.panel.photo_ref = net_plot
        self.panel.pack()
        
    # Create information of pipe format
    def format_infopipe(self):
        global network
        if network == "":
            messagebox.showerror("Unspecified file", "No file has been selected")
        else:
            headers = pnd.DataFrame(columns=["Pipes","Material","Age"])
            headers["Pipes"]= self.wn.pipe_name_list
            with pnd.ExcelWriter('PipeInfo_'+net_name+'.xlsx') as writer:  
                headers.to_excel(writer, index =False)
            messagebox.showinfo("Created file","Pipe information file has been successfully created")

    # Create measrument format
    def format_measurements(self):
        global network
        if network == "":
            messagebox.showerror("Unspecified file", "No file has been selected")
        else:
            global time_step
            global total_duration
            time_step = int(self.wn.options.time.hydraulic_timestep/60)
            total_duration = int(self.wn.options.time.duration/3600+time_step/60)
            s = datetime.datetime.strptime('12:00 AM', '%I:%M %p') 
            w = []
            w.append(s.strftime('%I:%M %p'))
            for i in range(time_step,60*total_duration+1,time_step):
                w.append((s+datetime.timedelta(minutes=i)).strftime('%I:%M %p'))
            headers = pnd.DataFrame(columns=["Type"])
            headers2 = pnd.DataFrame(columns=["ID"])
            headers2["ID"] = w
            with pnd.ExcelWriter('Measurements_'+net_name+'.xlsx') as writer:  
                headers.to_excel(writer, index =False)
                headers2.to_excel(writer,startrow=1, index =False)
            messagebox.showinfo("Created file","Measurements file has been successfully created")
            
    # Load pipe information and save it as tag for each pipe
    def load_infopipe(self):
        global format_infopipe
        global network
        global materials
        global pipeinfo
        if network == "":
            messagebox.showerror("Unspecified file", "No file has been selected")
        else:
            format_infopipe = filedialog.askopenfilename(initialdir="/",title = "Open",
                                         filetype = (("Format files(.xlsx)",".xlsx"),("All files","*.*")))
            sheet = pnd.ExcelFile(format_infopipe).sheet_names
            xl = pnd.read_excel(format_infopipe,sheet_name=None)
            for n in sheet:
                if xl[n].columns[0] != 'Pipes' or xl[n].columns[2] != 'Age':
                    messagebox.showerror("Invalid file information","The input file is invalid: Headers incorrect")
                    format_infopipe = ""
                    break
                if len(xl[n].iloc[:,0]) != len(self.wn.pipe_name_list):
                    messagebox.showerror("Invalid file information","The input file is invalid: Pipe names incorrect")
                    format_infopipe = ""
                    break
                if pnd.isnull(xl[n].iloc[:,1]).any() == True or np.isnan(xl[n].iloc[:,2]).any() == True:
                    messagebox.showerror("Invalid file information","The input file is invalid: Missing information in the age or material")
                    format_infopipe = ""
                    break
                t = 0
                for pipe_name, pipe in self.wn.pipes():
                    if xl[n].iloc[t,1] in materials and xl[n].iloc[t,2] <= 50:
                        pipe.tag = xl[n].iloc[t,1]+"_"+str(int(xl[n].iloc[t,2]))
                        t += 1
                    else:
                        messagebox.showerror("Invalid file information","The input file is invalid: Material or age incorrect in some pipes")
                        format_infopipe = ""
                        break
                if format_infopipe != "":
                    messagebox.showinfo("Created file","Pipe information file has been successfully added")
                    self.wn.write(self.wn.name)
                    pipeinfo = True
    # Cargar Formato de mediciones
    def load_measurements(self):
        global format_measurement
        global network
        global Real_data
        global Elements
        global total_duration
        global time_step
        global measurements
        if network == "":
            messagebox.showerror("Unspecified file", "No file has been selected")
        else:
            time_step = int(self.wn.options.time.hydraulic_timestep/60)
            total_duration = int(self.wn.options.time.duration/3600+time_step/60)
            format_measurement = filedialog.askopenfilename(initialdir="/",title = "Open",
                                         filetype = (("Format files(.xlsx)",".xlsx"),("All files","*.*")))
            sheet = pnd.ExcelFile(format_measurement).sheet_names
            xl = pnd.read_excel(format_measurement,sheet_name=None,header=None)
            x =[]
            y = []
            for n in sheet:
                if xl[sheet[0]].iloc[0,0] != 'Type' or xl[sheet[0]].iloc[1,0] != 'ID':
                    format_measurement = ""
                    messagebox.showerror("Invalid file information","The input file is invalid: Headers incorrect")
                    break
                if len(xl[sheet[0]].iloc[:,0]) != total_duration/(time_step/60)+3:
                    messagebox.showerror("Invalid file information","The input file is invalid: Time lenght incorrect")
                    format_measurement = ""
                    break

                if all(p == "T" or p == "N" for p in xl[sheet[0]].iloc[0,1:]):
                    x.extend(xl[sheet[0]].iloc[0,1:])
                else:
                    format_measurement = ""
                    messagebox.showerror("Invalid file information","The input file is invalid: Element type incorrect")
                    break

                if all(str(p) in self.wn.pipe_name_list or str(p) in self.wn.junction_name_list for p in xl[sheet[0]].iloc[1,1:]):
                    y.extend([str(i) for i in xl[sheet[0]].iloc[1,1:]])
                else:
                    messagebox.showerror("Invalid file information","The input file is invalid: Pipe name incorrect")
                    format_measurement = ""
                    break 
                for a in xl[sheet[0]].columns[1:]:
                    Real_data.extend(xl[sheet[0]][a][2:])

                if np.isnan(Real_data).any() == True:
                    messagebox.showerror("Invalid file information","The input file is invalid: Missing information in the measurements")
                    format_measurement = ""
                    break
                if format_measurement != "":
                    Elements["Type"] = x
                    Elements["ID"] = y
                    messagebox.showinfo("Created file","Measurements file has been successfully added")
                    measurements = True  
    
    # Load and exectue validation methodology                
    def validate(self):
        global network
        global temp
        global format_measurement
        if network == "" or format_measurement == "":
            messagebox.showerror("Unspecified file", "No network or measurement file has been selected")
        else:
            warning(temp, format_measurement)
            messagebox.showinfo("Validation Before Calibration", "Validation Before Calibration process has been successfully ended")
    
    
    #Load and execute data resample methodology
    def resample_data(self):
        global network
        global temp
        global data_file
        global flow_unit_factor
        global pressure_unit_factor
        global time_config
        if network == "":
            messagebox.showerror("Unspecified file", "No file has been selected")
        else:
            self.resamp_screen.destroy()
            print(data_file,flow_unit_factor.get(),pressure_unit_factor.get(),time_config.get())
            data_resample(temp, data_file, flow_unit_factor = flow_unit_factor.get(), pressure_unit_factor = pressure_unit_factor.get(), time_config= time_config.get())
            messagebox.showinfo("Resample Data", "Resample Data process has been successfully ended")
    
    # Load range and path information to the data resample methodology
    def resample_data_rang(self):
        global network
        global temp
        global data_file
        global flow_unit_factor
        global pressure_unit_factor
        global time_config
        if network == "":
            messagebox.showerror("Unspecified file", "No file has been selected")
        else:
            data_file = filedialog.askopenfilename(initialdir="/",title = "Open",
                                         filetype = (("Format files(.csv)",".csv"),("All files","*.*")))
            times= ['5T','10T','15T','30T','H','2H']
            self.resamp_screen = Toplevel(self.root)
            self.resamp_screen.iconbitmap('Logos/Icono.ico')
            self.resamp_screen.title("Resample_data_range")
            Label(self.resamp_screen,text="Flow unit factor").grid(row=0, column=0, padx=10, pady=10)
            Entry(self.resamp_screen,textvariable=flow_unit_factor).grid(row=0,column=1, padx=10, pady=10)
            Label(self.resamp_screen,text="Pressure unit factor").grid(row=1, column=0, padx=10, pady=10)
            Entry(self.resamp_screen,textvariable=pressure_unit_factor).grid(row=1,column=1, padx=10, pady=10)
            Label(self.resamp_screen,text="time_config").grid(row=2, column=0, padx=10, pady=10)
            Combobox(self.resamp_screen, state="readonly",values=times,textvariable=time_config).grid(row=2,column=1, padx=10, pady=10)
            Button(self.resamp_screen,text="OK",command=self.resample_data).grid(row=3,column=0, padx=10, pady=10)
        
    #Execute sensitivity analysis and loading bar in consecutive order
    def bar_init(self):
        global network
        if network == "":
            messagebox.showerror("Unspecified file", "No file has been selected")
        else:
            self.senscarga = Toplevel()
            self.senscarga.iconbitmap('Logos/Icono.ico')
            self.senscarga.title("Load Sensitivity Analysis")
            self.senscarga.geometry("250x30")
            self.load_bar = Progressbar(self.senscarga)
            self.load_bar.pack()
            self.root.update()
            # first layer of isolation, note var being passed along to the self.start_bar function
            # target is the function being started on a new thread, so the "bar handler" thread
            self.start_bar_thread = threading.Thread(target=self.start_bar, args=())
            # start the bar handling thread
            self.start_bar_thread.start()
    # Create loading bar
    def start_bar(self):
        # the load_bar needs to be configured for indeterminate amount of bouncing
        self.load_bar.config(mode='indeterminate', value=0, length = 200)
        # 8 here is for speed of bounce
        self.load_bar.start(8)
        # start the work-intensive thread, again a var can be passed in here too if desired
        self.work_thread = threading.Thread(target=self.work_task, args=())
        self.work_thread.start()
        # close the work thread
        self.work_thread.join()
        # stop the indeterminate bouncing
        self.load_bar.stop()
        # reconfigure the bar so it appears reset
        self.senscarga.destroy()
    # Load sensitivity information and execute methodology
    def work_task(self):
        global network
        global net_name
        sensitivity(network, net_name)
        messagebox.showinfo("Sensitivity Analysis", "Sensitivity Analysis has been successfully ended")

    # Execute calibration and loading bar in consecutive order
    def bar_init_ex(self):
        global temp
        if temp == "":
            messagebox.showerror("Unspecified file", "No file has been selected")
        else:
            self.rep.destroy()
            self.excarga = Toplevel()
            self.excarga.iconbitmap('Logos/Icono.ico')
            self.excarga.title("Load Calibration Process")
            self.excarga.geometry("250x30")
            Label(self.excarga,text="Loading...").grid(row=0, column=0, padx=10, pady=10)
            #self.load_bar = Progressbar(self.excarga)
            #self.load_bar.pack()
            self.root.update()
            # first layer of isolation, note var being passed along to the self.start_bar function
            # target is the function being started on a new thread, so the "bar handler" thread
            self.start_bar_thread = threading.Thread(target=self.start_bar_ex, args=())
            # start the bar handling thread
            self.start_bar_thread.start()
    
    # Start loadig bar
    def start_bar_ex(self):
        # the load_bar needs to be configured for indeterminate amount of bouncing
        #self.load_bar.config(mode='indeterminate', value=0, length = 200)
        # 8 here is for speed of bounce
        #self.load_bar.start(8)
        # start the work-intensive thread, again a var can be passed in here too if desired
        
        self.work_thread = threading.Thread(target=self.num_rep, args=())
        self.work_thread.start()
        # close the work thread
        self.work_thread.join()
        # stop the indeterminate bouncing
        #self.load_bar.stop()
        # reconfigure the bar so it appears reset
        self.excarga.destroy()
    
    # Calibracion Exactly: Execute calibration methodology
    def num_rep(self):
        global temp
        global net_name
        global Real_data
        global Elements
        
        global generetions_total
        global population_total
        
        global measurements
        global pipeinfo
        
        global calib_roughness
        global calib_minorloss
        global calib_demand_pattern
        global calib_demand
        global calib_leak
        global calib_valves_coeff
        global calib_valves_stat
        
        global rang_valve_coeff
        global rang_minorloss
        global rang_leak
        global rang_demand
        
        calibration(temp, net_name, Real_data, Elements,
                     calib_roughness, calib_minorloss, calib_demand, calib_leak, calib_valves_coeff, calib_valves_stat, calib_demand_pattern,
                     generations = generetions_total.get(), population= population_total.get(),
                     rang_valve_coeff = rang_valve_coeff, rang_minorloss = rang_minorloss, rang_leak = rang_leak, rang_demand = rang_demand)
        messagebox.showinfo("Calibration process", "Calibration process has been successfully ended")
    
    # Calibracion Exactly: Define number of generations and population
    def def_rang_calib(self):
        global rang_valve_coeff
        global rang_minorloss
        global rang_leak
        global rang_demand
        
        global valves_coeff_min
        global minorloss_min
        global leak_coef_min
        global dem_min
        
        global valves_coeff_max
        global minorloss_max
        global leak_coef_max
        global dem_max
        
        global generetions_total
        global population_total
        
        if calib_minorloss == True:
            rang_minorloss = [float(minorloss_min.get()),float(minorloss_max.get())]
        if calib_valves_coeff == True:
            rang_valve_coeff =[float(valves_coeff_min.get()),float(valves_coeff_max.get())]
        if calib_demand == True:
            rang_demand = [float(dem_min.get()),float(dem_max.get())]
        if calib_leak == True:                
            rang_leak = [float(leak_coef_min.get()),float(leak_coef_max.get())]
        self.rango.destroy()
        self.rep = Toplevel(self.root)
        self.rep.iconbitmap('Logos/Icono.ico')
        self.rep.title("Population settings")
        Label(self.rep,text="Number of generations").grid(row=0, column=0, padx=10, pady=10)
        Entry(self.rep,textvariable=generetions_total).grid(row=0,column=1, padx=10, pady=10)
        Label(self.rep,text="Population size").grid(row=1, column=0, padx=10, pady=10)
        Entry(self.rep,textvariable=population_total).grid(row=1,column=1, padx=10, pady=10)
        Button(self.rep,text="OK",command=self.bar_init_ex).grid(row=2,column=0, padx=10, pady=10)
    
    # Calibracion Exactly: Define ranges of the parameters
    def def_var(self):
        global calib_roughness
        global calib_minorloss
        global calib_demand_pattern
        global calib_demand
        global calib_leak
        global calib_valves_coeff
        global calib_valves_stat
        
        global roughness
        global minorloss
        global demand_pattern
        global demand
        global leak
        global valves_coeff
        global valves_stat
        
        global valves_coeff_min
        global minorloss_min
        global leak_coef_min
        global dem_min
        
        global valves_coeff_max
        global minorloss_max
        global leak_coef_max
        global dem_max       
        #roughness
        if roughness.get() == 1:
            calib_roughness = True 
        else:
            calib_roughness = False
        #minorloss
        if minorloss.get() == 1:
            calib_minorloss = True 
        else:
            calib_minorloss = False
        #demand_pattern
        if demand_pattern.get() == 1:
            calib_demand_pattern = True 
        else:
            calib_demand_pattern = False
        #Demand
        if demand.get() == 1:
            calib_demand = True 
        else:
            calib_demand = False
        #leak
        if leak.get() == 1:
            calib_leak = True 
        else:
            calib_leak = False
        #valves_coeff
        if valves_coeff.get() == 1:
            calib_valves_coeff = True 
        else:
            calib_valves_coeff = False
        #valves_status
        if valves_stat.get() == 1:
            calib_valves_stat = True 
        else:
            calib_valves_stat = False
        self.config.destroy()
        self.rango = Toplevel(self.root)
        self.rango.iconbitmap('Logos/Icono.ico')
        self.rango.title("Calibration range")
        
        Label(self.rango,text="Parameter").grid(row=0, column=0, padx=10, pady=10)
        Label(self.rango,text="Min").grid(row=0, column=1, padx=10, pady=10)
        Label(self.rango,text="Max").grid(row=0, column=2, padx=10, pady=10)
        if calib_minorloss == True:
            Label(self.rango,text="Minor loss").grid(row=1, column=0, padx=10, pady=10)
            Entry(self.rango,textvariable=minorloss_min).grid(row=1,column=1, padx=10, pady=10)
            Entry(self.rango,textvariable=minorloss_max).grid(row=1,column=2, padx=10, pady=10)
        if calib_valves_coeff == True:
            Label(self.rango,text="Valve coefficient").grid(row=2, column=0, padx=10, pady=10)
            Entry(self.rango,textvariable=valves_coeff_min).grid(row=2,column=1, padx=10, pady=10)
            Entry(self.rango,textvariable=valves_coeff_max).grid(row=2,column=2, padx=10, pady=10)
        if calib_demand == True:
            Label(self.rango,text="Demand factor").grid(row=3, column=0, padx=10, pady=10)
            Entry(self.rango,textvariable=dem_min).grid(row=3,column=1, padx=10, pady=10)
            Entry(self.rango,textvariable=dem_max).grid(row=3,column=2, padx=10, pady=10)
        if calib_leak == True:
            Label(self.rango,text="Leak coeff").grid(row=4, column=0, padx=10, pady=10)
            Entry(self.rango,textvariable=leak_coef_min).grid(row=4,column=1, padx=10, pady=10)
            Entry(self.rango,textvariable=leak_coef_max).grid(row=4,column=2, padx=10, pady=10)

        Button(self.rango,text="OK",command=self.def_rang_calib).grid(row=6,column=1, padx=10, pady=10)
    
    # Calibracion Exactly: Define calibrated parameters 
    def calibracion(self):
        global temp
        global net_name
        global roughness
        global minorloss
        global demand_pattern
        global demand
        global leak
        global valves_coeff
        global valves_stat
        
        if temp == "":
            messagebox.showerror("Unspecified file", "No input file has been selected")
        elif measurements == False:
            messagebox.showerror("Unspecified file", "No measurements file has been selected")
        elif pipeinfo == False:
            messagebox.showerror("Unspecified file", "No pipe information file has been selected")
        else:
            self.config = Toplevel(self.root)
            self.config.iconbitmap('Logos/Icono.ico')
            self.config.title("Calibration parameters")
            Label(self.config, text = "Choose the calibration parameters", width = 50, anchor = "center").pack()
            Checkbutton(self.config, text = "Roughness", variable= roughness).pack()
            Checkbutton(self.config, text = "Minor loss coeff", variable = minorloss).pack()
            Checkbutton(self.config, text = "Valves setting", variable = valves_coeff).pack()
            Checkbutton(self.config, text = "Valve state",variable = valves_stat).pack()
            Checkbutton(self.config, text = "Demand", variable=demand).pack()
            Checkbutton(self.config, text = "Emitter coeff",variable = leak).pack()
            Checkbutton(self.config, text = "Demand pattern", variable = demand_pattern).pack()
                        
            Button(self.config, text="Select", command=self.def_var).pack()

aplicacion1=Aplication()