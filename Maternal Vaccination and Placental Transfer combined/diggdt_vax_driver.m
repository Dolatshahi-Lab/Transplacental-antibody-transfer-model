%% Driver file to run vaccine simulation model from Erdogan, R. 2023. %%%%%
%% Updated:  04/04/2023, RE

%% Model setup - call parameters_Erdogan.m.  This file contains optimized 
%% model parameters used in (Erdogan, R. 2023).
%% You can change individual parameters in parameters_Erdogan.m to explore 
%% their impact on vaccine-induced IgG transfer.
parameters_Erdogan;

%% Initialize model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0 = [p.igg1,p.igg2,p.igg3,p.igg4,zeros(1,56),100 ,zeros(1,2)];

% Select a gestational week to simulate vaccination.
p.tvax = 25; % vaccination time (weeks gestational age)

%% Define color palette for plotting downstream figure %%%%%%%%%%%%%%%%%%%%

ColorOrder = [0.87, 0.443, 0.3569; 0.706, 0.87, 0.286; 0.302, 0.851, 1; 0.251, 0, 1];

%% Solve the model. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solution = ode23s(@(t,x) diggdt_vax(t,x,p,1), tspan, y0);

%% Plot fetal IgG subclass levels following vaccination at tvax. %%%%%%%%%%
figure(1)
subplot(1,4,1)
plot(solution.x,solution.y(57,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG1'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG1 (M)')
subplot(1,4,2)
plot(solution.x,solution.y(58,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG2'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG2 (M)')
subplot(1,4,3)
plot(solution.x,solution.y(59,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG3'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG3 (M)')
subplot(1,4,4)
plot(solution.x,solution.y(60,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG4'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG4 (M)')

%% Plot maternal IgG subclass levels following vaccination at tvax. %%%%%%%
figure(2)
subplot(1,4,1)
plot(solution.x,sum(solution.y(31:34,:)),'linewidth',2,'color',ColorOrder(1,:))
title('IgG1'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG1 (M)')
subplot(1,4,2)
plot(solution.x,solution.y(32,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG2'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG2 (M)')
subplot(1,4,3)
plot(solution.x,solution.y(33,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG3'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG3 (M)')
subplot(1,4,4)
plot(solution.x,solution.y(34,:),'linewidth',2,'color',ColorOrder(1,:))
title('IgG4'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG4 (M)')