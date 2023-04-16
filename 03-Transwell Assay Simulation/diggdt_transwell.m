function diggdt_transwell = diggdt_transwell(t,x,p,x0)
%% This is the ODE model of IgG transcytosis in HUVEC in a Transwell system.
%% Parameters are initialized in 'parameters_Transwell.m'.
%% This model can be used to recreate Figure 4C-D in Erdogan (2023).

%% Initialize the solution. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diggdt_transwell = zeros(7,1);

%% Main ODE routine. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diggdt_transwell(1) = - p.kon1*x(2)*x(1) + (p.koff1)*x(3); %internalized IgG (intracellular)
diggdt_transwell(2) = - p.kon1*x(2)*x(1) - p.kon2*x(2)*x(5) + (p.koff1 + p.kT)*x(3) + (p.koff2 + p.kT)*x(6); %unbound FcRn (intracellular)
diggdt_transwell(3) = p.kon1*x(2)*x(1) - (p.koff1 + p.kT + p.kdeg)*x(3); % IgG1-FcRn complexes (intracellular)
diggdt_transwell(4) = p.kT*x(3); %basolateral IgG1

diggdt_transwell(5) = - p.kon2*x(2)*x(5) + (p.koff1)*x(6); %unbound IgG4
diggdt_transwell(6) = p.kon2*x(2)*x(5) - (p.koff2 + p.kT + p.kdeg)*x(6); % IgG4-FcRn complexes
diggdt_transwell(7) = p.kT*x(6); %basolateral IgG4

end

