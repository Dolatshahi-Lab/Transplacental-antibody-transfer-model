%% Code that plots Figure 1 in Erdogan, 2023.

parameters_Erdogan;
ColorOrder = [0.87, 0.443, 0.3569; 0.706, 0.87, 0.286; 0.302, 0.851, 1; 0.251, 0, 1];

% sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0);
% 
% sol_eval = deval(sol,time_points)';
% ratios = sol_eval(end,27:30)'./sol_eval(end,1:4)';
% figure; plot(sol_real.x,sol_real.y(27:end,:))

%% Build LHS matrix (code adapted from Marino, S. et al)
%% Each parameter varied +/- 25% of its optimized value.
runs=200; %number of iterations (size of matrix)
% (xmin,xmean,xmax,xsd,nsample,distrib,threshold)

k_up_LHS = LHS_Call_RE(p.k_up - p.k_up*0.25, p.k_up, p.k_up + p.k_up*0.25, 1e-2 ,runs,'unif'); 
k_deg_LHS = LHS_Call_RE(p.k_deg - p.k_deg*0.25, p.k_deg, p.k_deg + p.k_deg*0.25, 1e-1,runs,'unif'); 
k_t_LHS = LHS_Call_RE(p.k_t - p.k_t*0.25, p.k_t, p.k_t + p.k_t*0.25, 1e-1, runs,'unif'); 
fcrn_LHS = LHS_Call_RE(p.fcrn - p.fcrn*0.25, p.fcrn, p.fcrn + p.fcrn*0.25,1e-7, runs,'unif'); 
fcgr2b_LHS = LHS_Call_RE(p.fcgr2b - p.fcgr2b*0.25,p.fcgr2b, p.fcgr2b + p.fcgr2b*0.25, 1e-7, runs,'unif'); 
v_m_LHS = LHS_Call_RE(p.v_m - 0.25*p.v_m, p.v_m, p.v_m + 0.25*p.v_m, 1, runs,'unif');
v_se_LHS = LHS_Call_RE(p.v_se - 0.25*p.v_se, p.v_se, p.v_se + 0.25*p.v_se, 1e-3, runs,'unif');
v_s_LHS = LHS_Call_RE(p.v_s - 0.25*p.v_s, p.v_s, p.v_s + 0.25*p.v_s, 1e-1, runs,'unif');
v_ee_LHS = LHS_Call_RE(p.v_ee - 0.25*p.v_ee, p.v_ee, p.v_ee + 0.25*p.v_ee, 1e-2, runs,'unif');
v_f_LHS = LHS_Call_RE(p.v_f - 0.25*p.v_f, p.v_f, p.v_f + 0.25*p.v_f, 1e-2, runs,'unif');

LHSmatrix=[k_up_LHS k_deg_LHS k_t_LHS fcrn_LHS fcgr2b_LHS v_m_LHS ...
    v_se_LHS v_s_LHS v_ee_LHS v_f_LHS];

%% Run global sensitivty analysis.
    ratday = [1 15 21]'*40/21; %rat days scaled to human weeks GA

for n=1:runs %Run solution x times choosing different values
    fcgr2b = [0 LHSmatrix(n,5)/2 LHSmatrix(n,5)]'; %relative FcRn expression in rat placenta
    fcrn = [0 LHSmatrix(n,4)/2 LHSmatrix(n,4)]'; %relative FcRn expression in rat placenta
    p.fcrn_curve = fit(ratday,fcrn,'poly2');
    p.fcgr2b_curve = fit(ratday,fcgr2b,'poly2');
    p.k_up = LHSmatrix(n,1);
    p.k_deg = LHSmatrix(n,2);
    p.k_t = LHSmatrix(n,3);
    p.v_m = LHSmatrix(n,6);
    p.v_se = LHSmatrix(n,7);
    p.v_s = LHSmatrix(n,8);
    p.v_ee = LHSmatrix(n,9);
    p.v_f = LHSmatrix(n,10);
    n
    sol= ode15s(@(t,x) diggdt(t,x,p), tspan, x0);
    A = deval(sol,time_points)';
    
    igg1_lhs(:,n)=A(:,27);
    igg2_lhs(:,n)=A(:,28);
    igg3_lhs(:,n)=A(:,29);
    igg4_lhs(:,n)=A(:,30);

    igg1_ratio(:,n)=A(:,27)./A(:,1);
    igg2_ratio(:,n)=A(:,28)./A(:,2);
    igg3_ratio(:,n)=A(:,29)./A(:,3);
    igg4_ratio(:,n)=A(:,30)./A(:,4);
end

% remove simulations that go below zero.
igg1_lhs(:,igg1_lhs(end,:)<0)=[];
igg2_lhs(:,igg2_lhs(end,:)<0)=[];
igg3_lhs(:,igg3_lhs(end,:)<0)=[];
igg4_lhs(:,igg4_lhs(end,:)<0)=[];

%% Plot optimized output against global sensitivty analysis runs.
parameters_Erdogan;
ColorOrder = [0.87, 0.443, 0.3569; 0.706, 0.87, 0.286; 0.302, 0.851, 1; 0.251, 0, 1];
y0 = [p.igg1,p.igg2,p.igg3,p.igg4,zeros(1,26)];

% p.k_up = 0.075;
% p.k_t = 0.07;
sol= ode15s(@(t,x) diggdt(t,x,p), tspan, y0);
sol_eval = deval(sol,time_points)';

ratios = sol_eval(end,27:30)./sol_eval(end,1:4)
% figure; plot(sol.x,sol.y(1:4,:))
% Data from Malek et al
tdata = [19.5,25,30,34.5,39]';
xCB.x1 = [0.93,2.12,3.7,5.65,10.43]/p.mass_igg; xCB.x2 = [0.31,0.74,0.93,1.19,1.56]/p.mass_igg;
xCB.x3 = [0.05,0.15,0.19,0.26,0.41]/(p.mass_igg+2e4); xCB.x4 = [0.04,0.13,0.21,0.25,0.47]/p.mass_igg;
xCB.x1sd = [0.42,0.76,0.92,1.52,3.27]/p.mass_igg; xCB.x2sd = [0.13,0.41,0.39,0.42,0.35]/p.mass_igg;
xCB.x3sd = [0.03,0.09,0.11,0.14,0.14]/(p.mass_igg+2e4); xCB.x4sd = [0.03,0.09,0.1,0.13,0.11]/p.mass_igg;
xCB.bulk = [1.44,3.32,5.57,8.47,11.98]/(p.mass_igg); xCB.bulksd = [0.67,1.29,1.1,1.71,2.18]/p.mass_igg;

figure;
subplot(1,4,1)
% plot(time_points,igg1_lhs,'color',[0.95 0.95 0.95],'linewidth',0.75); hold on
plot(time_points,sol_eval(:,27),'linewidth',2,'color',ColorOrder(1,:)); hold on
errorbar(tdata,xCB.x1,xCB.x1sd,'linestyle','none','color','k')
scatter(tdata,xCB.x1,50,'o','markerfacecolor',ColorOrder(1,:),'markeredgecolor','k','HandleVisibility','off')
xlabel('Gestational Age (weeks)'); ylabel('Cord IgG1 (M)');ylim([0 10e-5])

subplot(1,4,2)
% plot(time_points,igg2_lhs,'color',[0.95 0.95 0.95],'linewidth',0.75); hold on
plot(time_points,sol_eval(:,28),'linewidth',2,'color',ColorOrder(2,:)); hold on% 
errorbar(tdata,xCB.x2,xCB.x2sd,'linestyle','none','color','k')
scatter(tdata,xCB.x2,50,'^','markerfacecolor',ColorOrder(2,:),'markeredgecolor','k','HandleVisibility','off')
xlabel('Gestational Age (weeks)'); ylabel('Cord IgG2 (M)');ylim([0 2e-5])

subplot(1,4,3)
% plot(time_points,igg3_lhs,'color',[0.95 0.95 0.95],'linewidth',0.75); hold on
plot(time_points,sol_eval(:,29),'linewidth',2,'color',ColorOrder(3,:)); hold on
errorbar(tdata,xCB.x3,xCB.x3sd,'linestyle','none','color','k')
scatter(tdata,xCB.x3,50,'square','markerfacecolor',ColorOrder(3,:),'markeredgecolor','k','HandleVisibility','off')
xlabel('Gestational Age (weeks)'); ylabel('Cord IgG3 (M)');ylim([0 5e-6])

subplot(1,4,4)
% plot(time_points,igg4_lhs,'color',[0.95 0.95 0.95],'linewidth',0.75); hold on
plot(time_points,sol_eval(:,30),'linewidth',2,'color',ColorOrder(4,:)); hold on
errorbar(tdata,xCB.x4,xCB.x4sd,'linestyle','none','color','k')
scatter(tdata,xCB.x4,50,'diamond','markerfacecolor',ColorOrder(4,:),'markeredgecolor','k','HandleVisibility','off')
xlabel('Gestational Age (weeks)'); ylabel('Cord IgG4 (M)');ylim([0 5e-6])
% 
% figure;
% plot(time_points,sum(sol_eval(:,27:30),2),'color',[0.5 0.5 0.5],'linewidth',2); hold on
% plot(time_points,sol_eval(:,27),'linewidth',2,'color',ColorOrder(1,:)); hold on
% plot(time_points,sol_eval(:,28),'linewidth',2,'color',ColorOrder(2,:)); hold on
% plot(time_points,sol_eval(:,29),'linewidth',2,'color',ColorOrder(3,:)); hold on
% plot(time_points,sol_eval(:,30),'linewidth',2,'color',ColorOrder(4,:)); hold on
% errorbar(tdata,xCB.bulk,xCB.bulksd,'linestyle','none','color','k')
% scatter(tdata,xCB.bulk,30,'o','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k','HandleVisibility','off')
% xlabel('Gestational Age (weeks)'); ylabel('Cord IgG (M)');
% legend('Total IgG','IgG1','IgG2','IgG3','IgG4')

%% Swarmchart comparing model-predicted ratios to published literature
% load('Model_LHS.mat')
clementsdata = readmatrix('clements et al frontiers 2020 fig 1 dig.csv'); clementsdata = clementsdata(1:21,2:5);
ColorOrder = [0.87, 0.443, 0.3569; 0.706, 0.87, 0.286; 0.302, 0.851, 1; 0.251, 0, 1];

figure;
% tiledlayout(1,2)
% nexttile
ax=gca; 
yline(1,'--','HandleVisibility','off'); hold on
swarmchart(ones(1,200),igg1_ratio(end,:),'o','markeredgecolor','none',...
    'markerfacecolor',ColorOrder(1,:),'MarkerFaceAlpha',0.5,'XJitterwidth',0.5)
swarmchart(3*ones(1,200),igg2_ratio(end,:),'^','markeredgecolor','none',...
    'markerfacecolor',ColorOrder(2,:),'MarkerFaceAlpha',0.5,'XJitterwidth',0.5,'handlevisibility','off')
swarmchart(5*ones(1,200),igg3_ratio(end,:),'square','markeredgecolor','none',...
    'markerfacecolor',ColorOrder(3,:),'MarkerFaceAlpha',0.5,'XJitterwidth',0.5,'handlevisibility','off')
swarmchart(7*ones(1,200),igg4_ratio(end,:),'diamond','markeredgecolor','none',...
    'markerfacecolor',ColorOrder(4,:),'MarkerFaceAlpha',0.5,'XJitterwidth',0.5,'handlevisibility','off')
xticks(1:2:7);xticklabels({'IgG1','IgG2','IgG3','IgG4'});
ylabel('F:M Ratio'); ax.FontSize = 14; xlim([0 5]); ylim([0 3])

scatter(1,clementsdata(end,1),50,'o','markeredgecolor','k','markerfacecolor',ColorOrder(1,:),'linewidth',2)
scatter(3,clementsdata(end,2),50,'^','markeredgecolor','k','markerfacecolor',ColorOrder(2,:),'linewidth',2)
scatter(5,clementsdata(end,3),50,'square','markeredgecolor','k','markerfacecolor',ColorOrder(3,:),'linewidth',2)
scatter(7,clementsdata(end,4),50,'diamond','markeredgecolor','k','markerfacecolor',ColorOrder(4,:),'linewidth',2)
scatter([1:2:8],[mean(igg1_ratio(end,:)),mean(igg2_ratio(end,:)),mean(igg3_ratio(end,:)),...
    mean(igg4_ratio(end,:))],500,'_','markeredgecolor','k','linewidth',2.5,'handlevisibility','off')

% nexttile; ax = gca;
yline(1,'--','HandleVisibility','off'); hold on
swarmchart(2*ones(1,21),clementsdata(:,1),'o','markerfacecolor','none',...
    'markeredgecolor',ColorOrder(1,:),'MarkerFaceAlpha',0.75,'XJitterwidth',0.25)
swarmchart(4*ones(1,21),clementsdata(:,2),'^','markerfacecolor','none',...
    'markeredgecolor',ColorOrder(2,:),'MarkerFaceAlpha',0.75,'XJitterwidth',0.25,'handlevisibility','off')
swarmchart(6*ones(1,21),clementsdata(:,3),'square','markerfacecolor','none',...
    'markeredgecolor',ColorOrder(3,:),'MarkerFaceAlpha',0.75,'XJitterwidth',0.25,'handlevisibility','off')
swarmchart(8*ones(1,21),clementsdata(:,4),'diamond','markerfacecolor','none',...
    'markeredgecolor',ColorOrder(4,:),'MarkerFaceAlpha',0.75,'XJitterwidth',0.25,'handlevisibility','off')
scatter([2:2:8],[mean(clementsdata(:,1)),mean(clementsdata(:,2)),...
    mean(clementsdata(:,3)),mean(clementsdata([1:18,20:21],4))],500,'_','markeredgecolor','k','linewidth',2.5,'handlevisibility','off')


xticks(1.5:2:7.5);xticklabels({'IgG1','IgG2','IgG3','IgG4'});
ylabel('F:M Ratio (Cord:Mat Ratio)'); ax.FontSize = 14; xlim([0 9]); ylim([0 3])
% legend('Simulation','Clements et al (2020)')

%% what range of outputs emerge from varied FcRn and 2b in calipro range?
parameters_Erdogan;
fcgr2b_options = linspace(3.096e-6,3.6078e-6,20);
fcrn_options = linspace(6.591e-6,6.794e-6,20);
% fcrn_options = logspace(-7,-5,20);

figure
for i = 1:length(fcgr2b_options)
     fcrn = [0 fcgr2b_options(i)/2 fcgr2b_options(i)]'; %relative FcRn expression in rat placenta
%         p.fcrn_curve = fit(days_in_rat,fcrn,'poly2');
        p.fcgr2b_curve = fit(days_in_rat,fcrn,'poly2');
        sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0);
    plot(sol.x,sol.y(27,:)); hold on
    igg_fcrn_exp(i) = sol.y(27,end)./sol.y(1,end);
end

figure
swarmchart(ones(1,length(fcgr2b_options)),igg_fcrn_exp,'XJitterWidth',0.2)
hold on; ylabel('F:M ratio'); xticks([1 2]); xticklabels({'FcgRIIb exp.','Patient data'})
swarmchart(2*ones(1,length(clementsdata)),clementsdata(:,1),'XJitterWidth',0.2)
title('Vary FcgRIIb - Calipro range')