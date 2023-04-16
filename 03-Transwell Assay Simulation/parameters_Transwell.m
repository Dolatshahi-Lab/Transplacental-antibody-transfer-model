%% Initialize parameters for the Transwell simulation.  
%% This file is called by diggdt_transwell_driver.m.
%% User can explore different experimental conditions by altering the initial
%% concentraitons of IgG subclasses (IgG10, IgG20), binding parameters (kon/koff), 
%% FcRn expression level, or time of experiment.

%% Binding parameters. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%IgG4-FcRn parameters.
p.kon1 = 1.36e5; %ml/(minmol) -> 1/(Mmin)
p.koff1 = 0.0068; %1/min
p.kd1 = p.koff1/p.kon1; %M

%IgG1-FcRn parameters.
p.kon2 = 5.4e5; %ml/(minmol) -> 1/(Mmin)
p.koff2 = 0.0068; %1/min
p.kd2 = p.koff2/p.kon2;

%% Model parameters. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.kdeg = 1e-2; %rate of lysosomal degradation (mL/min)
p.kT = 2.5e-5; %rate of transcytosis across HUVEC (mL/min)
vol = 0.1; %volume of apical chamber (mL)
IgG10 = 3.33e3; %initial apical concentration of IgG1 (nM)
IgG40 = 3.33e3; %initial apical concentration of IgG4 (nM)
FcRn0 = 2.2e3; %FcRn expression level (nM)

%% Initialize model. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [IgG10,FcRn0,0,0,IgG40,0,0]; 

