%% Driver file to run compartmental placental model from Erdogan, R. 2023.
%% Updated:  04/04/2023, RE

%% Model setup - call parameters_Erdogan.m.  This file contains optimized model parameters used in (Erdogan, R. 2023).
%% You can change individual parameters in parameters_Erdogan.m to explore their affect on IgG transfer.
parameters_Erdogan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0 = [p.igg1,p.igg2,p.igg3,p.igg4,zeros(1,26)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define color palette for plotting downstream figure %%%%%%%%%%%%%%%%%%%%
ColorOrder = [0.87, 0.443, 0.3569; 0.706, 0.87, 0.286; 0.302, 0.851, 1; 0.251, 0, 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve and plot fetal IgG levels over time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solution= ode15s(@(t,x) diggdt(t,x,p), tspan, y0);

figure(1)
subplot(1,4,1)
plot(solution.x,solution.y(27,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG1'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG1 (M)')
subplot(1,4,2)
plot(solution.x,solution.y(28,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG2'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG2 (M)')
subplot(1,4,3)
plot(solution.x,solution.y(29,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG3'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG3 (M)')
subplot(1,4,4)
plot(solution.x,solution.y(30,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG4'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG4 (M)')
