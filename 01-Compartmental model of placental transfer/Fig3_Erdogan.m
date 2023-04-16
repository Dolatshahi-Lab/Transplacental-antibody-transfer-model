%% Code that plots Figure 3 in Erdogan, 2023.
ColorOrder = [0.87, 0.443, 0.3569; 0.706, 0.87, 0.286; 0.302, 0.851, 1; 0.251, 0, 1];
odeopts = odeset('RelTol',1e-10,'AbsTol',1e-10);
%% Fig 3A. 
%% plot with varying FcRn. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters_Erdogan;
p.k_d_fcgr2b = p.k_d_fcrn; 

%Re-compute binding parameters:
p.kon1b = 5.4e6; %(1/Ms)
p.koff1b = p.k_d_fcgr2b(1)*p.kon1b; %(1.4 1/s)
%IgG2
p.kon2b = p.koff1b/p.k_d_fcgr2b(2); %(1/Ms)
p.koff2b = p.k_d_fcgr2b(2)*p.kon2b; %(1.4 1/s)
%IgG3
p.kon3b = p.koff1b/p.k_d_fcgr2b(3); %(1/Ms)
p.koff3b = p.k_d_fcgr2b(3)*p.kon3b; %(1/s)
%IgG4
p.kon4b = p.koff1b/p.k_d_fcgr2b(4); %(1/Ms)
p.koff4b = p.k_d_fcgr2b(4)*p.kon4b; %(1/s)

varyFcrn = logspace(-7,-4,40); 
clear ratio; 
for i = 1:length(varyFcrn)
    fcrn = [0 varyFcrn(i)/2 varyFcrn(i)]'; %relative FcRn expression in rat placenta
    p.fcgr2b_curve = fit(ratday,fcrn,'poly2');
    sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,odeopts);
    ratio(:,i) = sol.y(27:30,end)./sol.y(1:4,end);
    i
end

figure; subplot(1,2,1)
semilogx(varyFcrn,ratio(1,:),'-','color',ColorOrder(1,:),'linewidth',2); hold on
semilogx(varyFcrn,ratio(2,:),'-','color',ColorOrder(2,:),'linewidth',2); 
semilogx(varyFcrn,ratio(3,:),'-','color',ColorOrder(3,:),'linewidth',2); 
semilogx(varyFcrn,ratio(4,:),'-','color',ColorOrder(4,:),'linewidth',2);
ylabel('Transfer Ratio')
xlabel('FcRn_{EC} (M)')
legend('IgG1','IgG2','IgG3','IgG4','location','best')

%% plot with varying FcgRIIb. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters_Erdogan;

varyFcgr2b = logspace(-7,-4,40); 
% figure;
clear ratio; 
for i = 1:length(varyFcgr2b)
    fcgr2b = [0 varyFcgr2b(i)/2 varyFcgr2b(i)]'; %relative FcRn expression in rat placenta
    p.fcgr2b_curve = fit(ratday,fcgr2b,'poly2');
    sol = ode15s(@(t,x) diggdt(t,x,p), tspan, x0,odeopts);
    ratio(:,i) = sol.y(27:30,end)./sol.y(1:4,end);
    i
end

subplot(1,2,2)
semilogx(varyFcgr2b,ratio(1,:),'-','color',ColorOrder(1,:),'linewidth',2); hold on
semilogx(varyFcgr2b,ratio(2,:),'-','color',ColorOrder(2,:),'linewidth',2); 
semilogx(varyFcgr2b,ratio(3,:),'-','color',ColorOrder(3,:),'linewidth',2); 
semilogx(varyFcgr2b,ratio(4,:),'-','color',ColorOrder(4,:),'linewidth',2);
ylabel('Transfer Ratio')
xlabel('Fc\gammaRIIb_{EC} (M)')
legend('IgG1','IgG2','IgG3','IgG4','location','best')
parameters_Erdogan; xline(p.fcgr2b,'k--','handlevisibility','off')


%% scRNA-seq data violin plots (Fig 3B-C). %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('vt_fcgr2b_fcrn.mat'); load('tsang_fcgr2b_fcrn.mat');

%integrate transcript data from Vento Tormo and Tsang data sets (ratios of
%endothelial cell FCGR2B and FCGRT expression)
intx.t1 = log(violin_vt.ec_fcgr2b_rn_ratio); %integrated transcripts t1
intx.t2 = log(violin.ec_pe);
intx.t3 = log(violin.ec_healthy);
figure; ax = gca; ax.FontSize = 14;
violinplot_SD(intx',{"Tri 1","Pre-E","Term"});
hold on
yline(0,'--'); ylabel('log(FCGR2B/FCGRT) in ECs'); ylim([-10 10])
xticklabels({'<12 weeks','30 weeks','38 weeks'}); 

anova_y = [intx.t1 intx.t2 intx.t3];

figure()
g = [ones(length(intx.t1),1); 2*ones(length(intx.t2),1); 3*ones(length(intx.t3),1)];
[~,~,stats] = anova1(anova_y,g,'off');
c1 = multcompare(stats,'alpha',0.05,'ctype','bonferroni'); %compare differences between col (IgG1 v IgG4)
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

violin2.ec_pe2b = violin.ec_pe2b;
violin2.ec_healthy2b= violin.ec_healthy2b; 
violin2.ec_peRN = violin.ec_peRN; 
violin2.ec_healthyRN = violin.ec_healthyRN;
figure
violinplot(violin2,["Pre-E FcgRIIb","Healthy FcgRIIb","Pre-E FcRn","Healthy FcRn"]);
ylabel('Counts');

%% plot FcRn and FcgRIIb combined model (Fig. S3) %%%%%%%%%%%%%%%%%%%%%%%%%
parameters_Erdogan;
x0 = [x0 zeros(1,5)];
FcRn_ec_test = logspace(-7,-5,100);
for i=1:length(FcRn_ec_test)
    fcrn_ec = [0 FcRn_ec_test(i)/2 FcRn_ec_test(i)]'; %relative FcRn expression in rat placenta
    p.fcrn_ec_curve = fit(ratday,fcrn_ec,'poly2');
    sol = ode15s(@(t,x) diggdt_FcgRIIb_and_FcRn(t,x,p), tspan, x0,odeopts);
    ratio(:,i) = sol.y(27:30,end)./sol.y(1:4,end);
end
figure; ax=gca;
semilogx(FcRn_ec_test/p.fcgr2b,ratio(1,:),'-','color',ColorOrder(1,:),'linewidth',2); hold on
semilogx(FcRn_ec_test/p.fcgr2b,ratio(2,:),'-','color',ColorOrder(2,:),'linewidth',2); 
semilogx(FcRn_ec_test/p.fcgr2b,ratio(3,:),'-','color',ColorOrder(3,:),'linewidth',2); 
semilogx(FcRn_ec_test/p.fcgr2b,ratio(4,:),'-','color',ColorOrder(4,:),'linewidth',2);
xline([1 0.4],'k--'); ax.FontSize=14;% xline(p.fcgr2b,'k--')
xlabel('FcRn_{EC} / FcgRIIb_{EC}'); ylabel('F:M Ratio');
%% scatter plot to compare FCGR2B and FcRn in scRNA-seq data (Fig S4) %%%%%

figure; tiledlayout(1,5); 

nexttile; ax = gca; ax.FontSize = 14;
line([0:12],[0:12],'linestyle','--','color',[0.75 0.75 0.75]); hold on
% mdl = fitlm((violin_vt.ec_fcgr2b+0.1),(violin_vt.ec_fcrn+0.1)); hold on
% plot(0:80,[0:80]*mdl.Coefficients.Estimate(2)+mdl.Coefficients.Estimate(1),'k-','linewidth',3,'HandleVisibility','off');
scatter((violin_vt.ec_fcrn),(violin_vt.ec_fcgr2b),25,'o','markeredgecolor','k','markerfacecolor',[0.5 0.5 0.5])
first_tri_ratio = (violin_vt.ec_fcgr2b+0.1)./(violin_vt.ec_fcrn+0.1);
title(append('1st Tri. ECs, ',string(floor(length(first_tri_ratio(first_tri_ratio>1))...
    /length(first_tri_ratio)*100)),'% FcgRIIb-enriched'))
ylabel('FcgRIIb (counts)'); xlabel('FcRn (counts)')
xlim([0 12]); ylim([0 12])

nexttile; ax = gca; ax.FontSize = 14;
line([0:80],[0:80],'linestyle','--','color',[0.75 0.75 0.75]); hold on
% mdl = fitlm((violin.ec_pe2b+0.1),(violin.ec_peRN+0.1)); hold on
% plot(0:80,[0:80]*mdl.Coefficients.Estimate(2)+mdl.Coefficients.Estimate(1),'k-','linewidth',3,'HandleVisibility','off');
scatter((violin.ec_peRN+0.1),(violin.ec_pe2b+0.1),25,'o','markeredgecolor','k','markerfacecolor','r')
% title(append('Pre-E ECs, R^2 = ',string(mdl.Rsquared.Ordinary)))
preE_ratio = (violin.ec_pe2b+0.1)./(violin.ec_peRN+0.1);
title(append('Pre-E ECs, ',string(floor(length(preE_ratio(preE_ratio>1))...
    /length(preE_ratio)*100)),'% FcgRIIb-enriched'))
ylabel('FcgRIIb (counts)'); xlabel('FcRn (counts)')
xlim([0 80]); ylim([0 80])

nexttile; ax = gca; ax.FontSize = 14;
line([0:80],[0:80],'linestyle','--','color',[0.75 0.75 0.75]); hold on
% mdl = fitlm((violin.ec_healthy2b+0.1),(violin.ec_healthyRN+0.1)); hold on
% plot(0:80,[0:80]*mdl.Coefficients.Estimate(2)+mdl.Coefficients.Estimate(1),'k-','linewidth',3,'HandleVisibility','off');
scatter((violin.ec_healthyRN+0.1),(violin.ec_healthy2b+0.1),25,'o','markeredgecolor','k','markerfacecolor','b')
% title(append('Term ECs, R^2 = ',string(mdl.Rsquared.Ordinary)))
healthy_ratio = (violin.ec_healthy2b+0.1)./(violin.ec_healthyRN+0.1);
title(append('Term ECs, ',string(floor(length(healthy_ratio(healthy_ratio>1))...
    /length(healthy_ratio)*100)),'% FcgRIIb-enriched'))
ylabel('FcgRIIb (counts)'); xlabel('FcRn (counts)')
xlim([0 80]); ylim([0 80])

nexttile; ax=gca; 
plot([1:3],[0,71,84],'o-','markerfacecolor','w','markeredgecolor','k','color','k','linewidth',1.5);hold on
xticks([1:3]); xticklabels({'1st Tri.','Pre-E','Term'}); ylabel('% ECs')
title('% FcgRIIb-enriched ECs');ax.FontSize = 14;
plot([1:3],[floor(length(first_tri_ratio(first_tri_ratio<1))...
    /length(first_tri_ratio)*100), floor(length(preE_ratio(preE_ratio<1))...
    /length(preE_ratio)*100), floor(length(healthy_ratio(healthy_ratio<1))...
    /length(healthy_ratio)*100)],'^-','markerfacecolor','w','markeredgecolor','k','color','k','linewidth',1.5)
legend('FcgRIIb','FcRn')

nexttile; ax=gca;
swarmchart(ones(1,length(violin_vt.ec_fcrn)),violin_vt.ec_fcrn,'markeredgecolor','k','markerfacecolor',[0.5 0.5 0.5])
xlim([0 2]); ylim([0 12]);ax.FontSize = 14; xticks(1); xticklabels({'FcRn'}); ylabel('Counts')