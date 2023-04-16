%% This driver file runs the ODE model of IgG transcytosis in HUVEC in vitro.
%% From Erdogan, 2023. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize the parameters by running the parameter file. %%%%%%%%%%%%%%%
%% Different experimental conditions can be explored by modifying these parameters.
%% For example, to recreate the simulation from Figure 4C, you can set the 
%% initial concentration of IgG1 equal to zero.
%% To recreate Figure 4D, you can run an in silico expeirment with increasing
%% concentraitons of IgG1 while fixing IgG4 at a steady concentration.

parameters_Transwell;
tspan = [0,120]; % Time span of simulation (minutes)

%% Solve the ODE model. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt = odeset('RelTol',1e-8,'AbsTol',1e-8,'NonNegative',1);
sol = ode23s(@(t,x) diggdt_transwell(t,x,p,x0), tspan, x0,opt);

%% Plot the results. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(sol.x,sol.y(4,:)/x0(1),'k-','linewidth',1.5);title('IgG4_{basolateral}'); hold on
plot(sol.x,sol.y(7,:)/x0(1),'r-','linewidth',1.5);title('IgG1_{basolateral}'); hold on
xlabel('Minutes'); ylabel('B:A ratio (%)');
legend('IgG4','IgG1')

