%% Plot Data and Model from Fig 4, Erdogan et al 2023.

%% Need to go through and make sure the data is in the right place.
%% 04/02/2023
%% set up experiment to test effect of concentration on trans efficiency
% IgG0 = 0.5*logspace(-3,0,20);
% IgG0 = [0.1 0.2 0.3 0.4 0.5]; %mg/ml
parameter_setup_transwell; clear IgG0 sol_eval
IgG0 = 0:0.01:0.5; 
figure
x0(5) = 0;
for i=1:length(IgG0)
    x0(1) = IgG0(i)*vol*g2kda/(IgG_MW*avo*vol*1e-3)*1e6; %mg/ml -> mg -> kDa -> number of particles -> mmol -> mmol/L -> (nM)
    sol = ode23s(@(t,x) dx_dt_simple(t,x,p,x0), tspan, x0,opt);
%     plot(sol.x,sol.y(3,:)/x0(1),'r-','linewidth',1.5); hold on
    sol_eval(i,:,:) = deval(sol,[15,30,60,120])/p.vb/g2kda*(IgG_MW*avo*p.vb*1e-3); %ng/ml
    sol_eval_temp = deval(sol,[15,30,60,120])/p.vb/g2kda*(IgG_MW*avo*p.vb*1e-3); %ng/ml
%     subplot(2,2,3); plot([15,30,60,120],sol_eval_temp(4,:)/IgG0(i)*100*1e-6+i/500,'o-','color',[0.2*i 0 0.5],'linewidth',2); hold on
%     xlabel('Minutes'); ylabel('A:B Ratio (%)')
%     subplot(2,2,4); plot([15,30,60,120],sol_eval_temp(4,:)+i*10,'o-','color',[0.2*i 0 0.5],'linewidth',2);  hold on
%     xlabel('Minutes'); ylabel('Basolateral IgG4 (ng/ml)')
%     eta(i) = x0(1)*x0(2)/(p.kd1*1e9*(x0(1)+x0(2)));
end

figure
plot(IgG0,sol_eval(:,4,1)*2.5,'-','color',[0.5 0.5 0.5],'linewidth',2); hold on
ylabel('(ng/ml)'); xlabel('Apical Concentration IgG4 (mg/ml)');%xticks([5e-3 1e-2 5e-2 1e-1 5e-1])
xline(0.33,':','linewidth',1.5);title('Basolateral IgG4 (t = 15 min)')
load('transwell_011223_original.mat')
% igg4_interp.data_interp = igg4_interp.data_interp([1:3,7:8,6,4:5,9:end]);
% igg4_interp.data_mean = igg4_interp.data_mean([1:3,7:8,6,4:5,9:end]);

scatter(igg4_conc,igg4_interp.data_interp,'o','markerfacecolor','k','markeredgecolor','k'); hold on
scatter([0.45,0.4,0.3,0.2,0.1],igg4_interp.data_mean,150,'_','linewidth',2,'markeredgecolor','k')
errorbar([0.45,0.4,0.3,0.2,0.1],igg4_interp.data_mean,igg4_interp.data_sd,'linestyle','none','color','k')

%% 2 IgGs in competition with each other
clear sol_eval IgG2
parameter_setup_transwell; 
x0(1) = 0.2*vol*g2kda/(IgG_MW*avo*vol*1e-3)*1e6;
% IgG1 = linspace(0.05,1,10);
IgG2 = linspace(0,1,30);
for i = 1:length(IgG2)
    x0(5) = IgG2(i)*vol*g2kda/(IgG_MW*avo*vol*1e-3)*1e6;
    sol = ode23s(@(t,x) dx_dt_simple(t,x,p,x0), tspan, x0,opt);
    sol_eval(i,:,:) = deval(sol,[15,30,60,120])/vol/g2kda*(IgG_MW*avo*vol*1e-3);
end
figure
subplot(1,2,1)
plot(IgG2,sol_eval(:,4,1),'ko-','linewidth',1.5); hold on
    plot(IgG2,sol_eval(:,7,1),'bo-','linewidth',1.5); ylim([0 600])
    xlabel('Apical IgG1 (ng/ml)'); ylabel('Basolateral IgG (ng/ml)')
legend('IgG4','IgG1','IgG4+IgG1','location','best')

subplot(1,2,2)
plot(IgG2,sol_eval(:,4,4),'ko-','linewidth',1.5); hold on
    plot(IgG2,sol_eval(:,7,4),'bo-','linewidth',1.5); ylim([0 600])
    xlabel('Apical IgG1 (ng/ml)'); ylabel('Basolateral IgG (ng/ml)')
legend('IgG4','IgG1','IgG4+IgG1','location','best')

%% overlay plot against experimental data

load('C:\Users\re8sb\Box\Remziye Erdogan share\Data and Analysis\Dolatshahi lab data\Transwell\020623\transwell_020623.mat')

igg1_conc = [zeros(1,2),ones(1,3)*0.1,ones(1,3)*0.2,ones(1,3)*0.4,ones(1,3)*0.8]; %ng/ml
igg4_conc = ones(1,14)*0.2'; %ng/ml

igg4_interp.data_mean(1) = mean(igg4_interp.data_interp(2:3));
igg4_interp.data_sd(1) = std(igg4_interp.data_interp(2:3));

figure; subplot(1,3,1)
plot(IgG2,sol_eval(:,4,1)*3,'k-','linewidth',1.5); hold on
    plot(IgG2,sol_eval(:,7,1)*3,'b-','linewidth',1.5); ylim([0 700]); xlim([-0.01 0.9])
    xlabel('Apical IgG1 (ng/ml)'); ylabel('Basolateral IgG (ng/ml)')
scatter([0,0.1,0.2,0.4,0.8],igg4_interp.data_mean,150,'_','linewidth',2,'markeredgecolor','k','handlevisibility','off'); hold on
errorbar([0,0.1,0.2,0.4,0.8],igg4_interp.data_mean,igg4_interp.data_sd,'linestyle','none','color','k','handlevisibility','off')
scatter(igg1_conc,igg4_interp.data_interp(2:end),'o','markerfacecolor','k','markeredgecolor','k')
scatter([0,0.1,0.2,0.4,0.8],[0 igg1_interp.data_mean],150,'_','linewidth',2,'markeredgecolor','b','handlevisibility','off')
errorbar([0,0.1,0.2,0.4,0.8],[0 igg1_interp.data_mean],[0 igg1_interp.data_sd],'linestyle','none','color','b','handlevisibility','off')
scatter(igg1_conc,igg1_interp.data_interp(2:end),'o','markerfacecolor','b','markeredgecolor','k')

subplot(1,3,2)
plot(IgG2,sol_eval(:,4,1)*3,'k-','linewidth',1.5); hold on
    plot(IgG2,sol_eval(:,7,1)*3,'b-','linewidth',1.5); ylim([0 700]); xlim([-0.01 0.9])
    xlabel('Apical IgG1 (ng/ml)'); ylabel('Basolateral IgG (ng/ml)')

subplot(1,3,3)
scatter([0,0.1,0.2,0.4,0.8],igg4_interp.data_mean,150,'_','linewidth',2,'markeredgecolor','k','handlevisibility','off'); hold on
errorbar([0,0.1,0.2,0.4,0.8],igg4_interp.data_mean,igg4_interp.data_sd,'linestyle','none','color','k','handlevisibility','off')
plot([0,0.1,0.2,0.4,0.8],igg4_interp.data_mean,'k-','handlevisibility','off')
scatter(igg1_conc,igg4_interp.data_interp(2:end),'o','markerfacecolor','k','markeredgecolor','k')
scatter([0,0.1,0.2,0.4,0.8],[0 igg1_interp.data_mean],150,'_','linewidth',2,'markeredgecolor','b','handlevisibility','off')
plot([0,0.1,0.2,0.4,0.8],[0 igg1_interp.data_mean],'b-','handlevisibility','off')
errorbar([0,0.1,0.2,0.4,0.8],[0 igg1_interp.data_mean],[0 igg1_interp.data_sd],'linestyle','none','color','b','handlevisibility','off')
scatter(igg1_conc,igg1_interp.data_interp(2:end),'o','markerfacecolor','b','markeredgecolor','k')
xlim([-0.01 0.9]); ylim([0 700])

%% Plot supplemental figure, mixing vs. not mixing. %%%%%%%%%%%%%%%%%%%%%%%
opt = odeset('RelTol',1e-8,'AbsTol',1e-8);

%% Fig S4. - effect of subclass mixing in physiological model.
parameters_Erdogan; 
tspan = [10 40];
x0 = [p.igg1,0,0,0,zeros(1,26)];
sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,opt);
igg1_ratio = sol.y(27,end)/sol.y(1,end);

x0 = [0,p.igg2,0,0,zeros(1,26)];
sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,opt);
igg2_ratio = sol.y(28,end)/sol.y(2,end);

x0 = [0,0,p.igg3,0,zeros(1,26)];
sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,opt);
igg3_ratio = sol.y(29,end)/sol.y(3,end);

x0 = [0,0,0,p.igg4,zeros(1,26)];
sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,opt);
igg4_ratio = sol.y(30,end)/sol.y(4,end);

x0 = [p.igg1,p.igg2,p.igg3,p.igg4,zeros(1,26)];
sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,opt);
igg1_ratio_comp = sol.y(27,end)/sol.y(1,end);
igg2_ratio_comp = sol.y(28,end)/sol.y(2,end);
igg3_ratio_comp = sol.y(29,end)/sol.y(3,end);
igg4_ratio_comp = sol.y(30,end)/sol.y(4,end);

figure; 
tiledlayout(1,2)
nexttile;ax = gca;
bar([igg1_ratio igg1_ratio_comp; igg2_ratio igg2_ratio_comp;...
    igg3_ratio igg3_ratio_comp; igg4_ratio igg4_ratio_comp])
ylabel('F:M Ratio'); xticklabels({'IgG1','IgG2','IgG3','IgG4'})
legend('Isolation','Competition','location','best'); ax.FontSize = 14;
hold on; yline(1,'--','handlevisibility','off')

nexttile; ax=gca;
bar([igg1_ratio igg1_ratio_comp; igg4_ratio igg4_ratio_comp])
ylabel('F:M Ratio'); xticklabels({'IgG1','IgG4'})
legend('Isolation','Competition','location','best'); ax.FontSize = 14;
hold on; yline(1,'--','handlevisibility','off')

[igg1_ratio igg1_ratio_comp; igg2_ratio igg2_ratio_comp;...
    igg3_ratio igg3_ratio_comp; igg4_ratio igg4_ratio_comp]

%% with varying FcgR2b (Fig  S4)
ColorOrder = [0.87, 0.443, 0.3569; 0.706, 0.87, 0.286; 0.302, 0.851, 1; 0.251, 0, 1];
opt = odeset('RelTol',1e-8,'AbsTol',1e-8);

parameters_Erdogan; 
Vary_fcgr2b = logspace(-7,-4,60);
for i = 1:length(Vary_fcgr2b)
    fcgr2b = [0 Vary_fcgr2b(i)/2 Vary_fcgr2b(i)]'; %relative FcRn expression in rat placenta
    p.fcgr2b_curve = fit(ratday,fcgr2b,'poly2');

tspan = [10 40];
x0 = [p.igg1,0,0,0,zeros(1,26)];
sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,opt);
igg1_ratio(i) = sol.y(27,end)/sol.y(1,end);

x0 = [0,p.igg2,0,0,zeros(1,26)];
sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,opt);
igg2_ratio(i) = sol.y(28,end)/sol.y(2,end);

x0 = [0,0,p.igg3,0,zeros(1,26)];
sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,opt);
igg3_ratio(i) = sol.y(29,end)/sol.y(3,end);

x0 = [0,0,0,p.igg4,zeros(1,26)];
sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,opt);
igg4_ratio(i) = sol.y(30,end)/sol.y(4,end);

x0 = [p.igg1,p.igg2,p.igg3,p.igg4,zeros(1,26)];
sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,opt);
igg1_ratio_comp(i) = sol.y(27,end)/sol.y(1,end);
igg2_ratio_comp(i) = sol.y(28,end)/sol.y(2,end);
igg3_ratio_comp(i) = sol.y(29,end)/sol.y(3,end);
igg4_ratio_comp(i) = sol.y(30,end)/sol.y(4,end);
end

figure; tiledlayout(1,3); 
nexttile;ax=gca;
semilogx(Vary_fcgr2b,igg1_ratio,'color',ColorOrder(1,:),'linewidth',2); hold on
semilogx(Vary_fcgr2b,igg2_ratio,'color',ColorOrder(2,:),'linewidth',2)
semilogx(Vary_fcgr2b,igg3_ratio,'color',ColorOrder(3,:),'linewidth',2)
semilogx(Vary_fcgr2b,igg4_ratio,'color',ColorOrder(4,:),'linewidth',2)
parameters_Erdogan; ylabel('F:M Ratio'); title('Isolated IgG Subclasses'); xlabel('FcgRIIb (M)')
xline([3.096e-6,3.6078e-6],'k--'); legend('IgG1','IgG2','IgG3','IgG4'); ax.FontSize = 14;

nexttile; ax=gca;
semilogx(Vary_fcgr2b,igg1_ratio_comp,'color',ColorOrder(1,:),'linewidth',2); hold on
semilogx(Vary_fcgr2b,igg2_ratio_comp,'color',ColorOrder(2,:),'linewidth',2)
semilogx(Vary_fcgr2b,igg3_ratio_comp,'color',ColorOrder(3,:),'linewidth',2)
semilogx(Vary_fcgr2b,igg4_ratio_comp,'color',ColorOrder(4,:),'linewidth',2)
parameters_Erdogan; ylabel('F:M Ratio'); title('Mixed IgG Subclasses'); xlabel('FcgRIIb (M)')
xline([3.096e-6,3.6078e-6],'k--'); legend('IgG1','IgG2','IgG3','IgG4'); ax.FontSize = 14;

nexttile; ax=gca;
semilogx(Vary_fcgr2b,igg1_ratio-igg1_ratio_comp,'color',ColorOrder(1,:),'linewidth',2); hold on
semilogx(Vary_fcgr2b,igg2_ratio-igg2_ratio_comp,'color',ColorOrder(2,:),'linewidth',2)
semilogx(Vary_fcgr2b,igg3_ratio-igg3_ratio_comp,'color',ColorOrder(3,:),'linewidth',2)
semilogx(Vary_fcgr2b,igg4_ratio-igg4_ratio_comp,'color',ColorOrder(4,:),'linewidth',2)
parameters_Erdogan; ylabel('\Delta F:M Ratio'); title('Isolated - Mixed'); xlabel('FcgRIIb (M)')
xline([3.096e-6,3.6078e-6],'k--'); legend('IgG1','IgG2','IgG3','IgG4'); ax.FontSize = 14;