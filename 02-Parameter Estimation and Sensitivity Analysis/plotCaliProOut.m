function plotCaliProOut(dotMatName)

ColorOrder = [0.87, 0.443, 0.3569; 0.706, 0.87, 0.286; 0.302, 0.851, 1; 0.251, 0, 1];

    CaliProOut = load(dotMatName);
%     pusheenPlot_pp(dotMatName);
    
% Data used for calibration from Malek et al %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tdata = [19.5,25,30,34.5,39]'; mass_igg = 1.5e5;
xCB.x1 = [0.93,2.12,3.7,5.65,10.43]/mass_igg; xCB.x2 = [0.31,0.74,0.93,1.19,1.56]/mass_igg;
xCB.x3 = [0.05,0.15,0.19,0.26,0.41]/(mass_igg+2e4); xCB.x4 = [0.04,0.13,0.21,0.25,0.47]/mass_igg;
xCB.x1sd = [0.42,0.76,0.92,1.52,3.27]/mass_igg; xCB.x2sd = [0.13,0.41,0.39,0.42,0.35]/mass_igg;
xCB.x3sd = [0.03,0.09,0.11,0.14,0.14]/(mass_igg+2e4); xCB.x4sd = [0.03,0.09,0.1,0.13,0.11]/mass_igg;
xCB.bulk = [1.44,3.32,5.57,8.47,11.98]/(mass_igg); xCB.bulksd = [0.67,1.29,1.1,1.71,2.18]/mass_igg;

    figure; subplot(1,4,1)
    plot(CaliProOut.odeSettings.tspan,CaliProOut.passedtotalModelRuns(:,:,1),'color',[0.9 0.9 0.9]); hold on
    plot(CaliProOut.odeSettings.tspan,mean(CaliProOut.passedtotalModelRuns(:,:,1),2),'color',ColorOrder(1,:),'linewidth',2)
    errorbar(CaliProOut.odeSettings.analysisTimePoints,xCB.x1,xCB.x1sd,'linestyle','none','color','k')
    scatter(CaliProOut.odeSettings.analysisTimePoints,xCB.x1,50,'markeredgecolor','k','markerfacecolor',ColorOrder(1,:))
    xlabel('Weeks'); ylabel('Fetal IgG (M)'); title('IgG1'); %ylim([0 1e-4])
        subplot(1,4,2)
        plot(CaliProOut.odeSettings.tspan,CaliProOut.passedtotalModelRuns(:,:,2),'color',[0.9 0.9 0.9]); hold on
    plot(CaliProOut.odeSettings.tspan,mean(CaliProOut.passedtotalModelRuns(:,:,2),2),'color',ColorOrder(2,:),'linewidth',2)
    errorbar(CaliProOut.odeSettings.analysisTimePoints,xCB.x2,xCB.x2sd,'linestyle','none','color','k')
    scatter(CaliProOut.odeSettings.analysisTimePoints,xCB.x2,50,'^','markeredgecolor','k','markerfacecolor',ColorOrder(2,:))
    xlabel('Weeks'); ylabel('Fetal IgG (M)'); title('IgG2'); %ylim([0 1e-4])
        subplot(1,4,3)
        plot(CaliProOut.odeSettings.tspan,CaliProOut.passedtotalModelRuns(:,:,3),'color',[0.9 0.9 0.9]); hold on
    plot(CaliProOut.odeSettings.tspan,mean(CaliProOut.passedtotalModelRuns(:,:,3),2),'color',ColorOrder(3,:),'linewidth',2)
    errorbar(CaliProOut.odeSettings.analysisTimePoints,xCB.x3,xCB.x3sd,'linestyle','none','color','k')
    scatter(CaliProOut.odeSettings.analysisTimePoints,xCB.x3,50,'diamond','markeredgecolor','k','markerfacecolor',ColorOrder(3,:))
    xlabel('Weeks'); ylabel('Fetal IgG (M)'); title('IgG3'); %ylim([0 1e-4])
        subplot(1,4,4)
        plot(CaliProOut.odeSettings.tspan,CaliProOut.passedtotalModelRuns(:,:,4),'color',[0.9 0.9 0.9]); hold on
    plot(CaliProOut.odeSettings.tspan,mean(CaliProOut.passedtotalModelRuns(:,:,4),2),'color',ColorOrder(4,:),'linewidth',2)
    errorbar(CaliProOut.odeSettings.analysisTimePoints,xCB.x4, xCB.x4sd,'linestyle','none','color','k')
    scatter(CaliProOut.odeSettings.analysisTimePoints,xCB.x4,70,'square','markeredgecolor','k','markerfacecolor',ColorOrder(4,:))
    xlabel('Weeks'); ylabel('Fetal IgG (M)'); title('IgG4'); %ylim([0 1e-4])
    
%% plot simulated transfer ratios
igg1ratio = CaliProOut.CalibrationOutputArray{1,end}.odeOut{27,1}./CaliProOut.CalibrationOutputArray{1,end}.odeOut{1,1};
igg2ratio = CaliProOut.CalibrationOutputArray{1,end}.odeOut{28,1}./CaliProOut.CalibrationOutputArray{1,end}.odeOut{2,1};
igg3ratio = CaliProOut.CalibrationOutputArray{1,end}.odeOut{29,1}./CaliProOut.CalibrationOutputArray{1,end}.odeOut{3,1};
igg4ratio = CaliProOut.CalibrationOutputArray{1,end}.odeOut{30,1}./CaliProOut.CalibrationOutputArray{1,end}.odeOut{4,1};

clementsdata = readmatrix(['C:\Users\re8sb\Box\Remziye Erdogan share\Data and Analysis\Antibody transfer model\manuscript code\Erdogan 2023\clements et al frontiers 2020 fig 1 dig.csv']); 
clementsdata = clementsdata(1:16,8:11);

figure;ax=gca; yline(1,'--','HandleVisibility','off'); hold on
swarmchart(ones(1,CaliProOut.odeSettings.NR),igg1ratio(end,:),'o','markeredgecolor','none',...
    'markerfacecolor',ColorOrder(1,:),'MarkerFaceAlpha',0.5,'XJitterwidth',0.75)
swarmchart(2*ones(1,CaliProOut.odeSettings.NR),igg2ratio(end,:),'^','markeredgecolor','none',...
    'markerfacecolor',ColorOrder(2,:),'MarkerFaceAlpha',0.5,'XJitterwidth',0.75,'handlevisibility','off')
swarmchart(3*ones(1,CaliProOut.odeSettings.NR),igg3ratio(end,:),'diamond','markeredgecolor','none',...
    'markerfacecolor',ColorOrder(3,:),'MarkerFaceAlpha',0.5,'XJitterwidth',0.75,'handlevisibility','off')
swarmchart(4*ones(1,CaliProOut.odeSettings.NR),igg4ratio(end,:),'square','markeredgecolor','none',...
    'markerfacecolor',ColorOrder(4,:),'MarkerFaceAlpha',0.5,'XJitterwidth',0.75,'handlevisibility','off')
scatter([1:4],[mean(igg1ratio(end,:)),mean(igg2ratio(end,:)),mean(igg3ratio(end,:)),...
    mean(igg4ratio(end,:))],CaliProOut.odeSettings.NR,'_','markeredgecolor','k','linewidth',3.5,'handlevisibility','off')
plot([1:4],[mean(igg1ratio(end,:)),mean(igg2ratio(end,:)),mean(igg3ratio(end,:)),...
    mean(igg4ratio(end,:))],'k-','linewidth',0.6,'handlevisibility','off')
% 
swarmchart(5*ones(1,16),clementsdata(:,1),'o','markerfacecolor','none',...
    'markeredgecolor',ColorOrder(1,:),'MarkerFaceAlpha',0.75,'XJitterwidth',0.5)
swarmchart(6*ones(1,16),clementsdata(:,2),'^','markerfacecolor','none',...
    'markeredgecolor',ColorOrder(2,:),'MarkerFaceAlpha',0.75,'XJitterwidth',0.5,'handlevisibility','off')
swarmchart(7*ones(1,15),clementsdata(1:15,3),'diamond','markerfacecolor','none',...
    'markeredgecolor',ColorOrder(3,:),'MarkerFaceAlpha',0.75,'XJitterwidth',0.5,'handlevisibility','off')
swarmchart(8*ones(1,16),clementsdata(:,4),'square','markerfacecolor','none',...
    'markeredgecolor',ColorOrder(4,:),'MarkerFaceAlpha',0.75,'XJitterwidth',0.5,'handlevisibility','off')
scatter([5:8],[mean(clementsdata(:,1)),mean(clementsdata(:,2)),...
    mean(clementsdata(1:15,3)),mean(clementsdata(:,4))],300,'_','markeredgecolor','k','linewidth',3.5,'handlevisibility','off')
plot([5:8],[mean(clementsdata(:,1)),mean(clementsdata(:,2)),...
    mean(clementsdata(1:15,3)),mean(clementsdata(:,4))],'k-','linewidth',0.6,'handlevisibility','off')

xticks(1:8);xticklabels({'IgG1','IgG2','IgG3','IgG4','IgG1','IgG2','IgG3','IgG4'});
ylabel('F:M Ratio'); ax.FontSize = 14; xlim([0 9])
legend('Simulation','Clements et al (2020)')



%% plot ratios as violin plot
% addpath(genpath('C:/Users/re8sb/PLSR-DA_MATLAB'))
% iggviolin.IgG1 = igg1ratio(end,:); iggviolin.IgG2 = igg2ratio(end,:);
% iggviolin.IgG3 = igg1ratio(end,:); iggviolin.IgG4 = igg1ratio(end,:);
% 
% clements.IgG1 = clementsdata(:,1); clements.IgG2 = clementsdata(:,2); 
% clements.IgG3 = clementsdata(:,3); clements.IgG4 = clementsdata(:,4);
% figure
% violinplot_SD(iggviolin,{'IgG1','IgG2','IgG3','IgG4'},'ShowData',true,'ShowMean',true); ylim([0 2.5])
% figure
% violinplot_SD(clements,{'IgG1','IgG2','IgG3','IgG4'},'ShowData',true,'ShowMean',true); ylim([0 2.5])
%% plot improved IgG1 fit over Calipro iterations
figure; j = 1;
for i = 1:2:14
    subplot(1,7,j)
    plot(CaliProOut.odeSettings.tspan,CaliProOut.CalibrationOutputArray{1,i}.odeOut{27,1},'color',[0.9 0.9 0.9]); hold on
    errorbar(CaliProOut.odeSettings.analysisTimePoints,xCB.x1,xCB.x1sd,'linestyle','none','color','k')
    scatter(CaliProOut.odeSettings.analysisTimePoints,xCB.x1,50,'markeredgecolor','k','markerfacecolor',ColorOrder(1,:))
    j=j+1;
    xlabel('weeks'); ylabel('fetal IgG1 (M)'); ylim([0 2e-4])
end

end