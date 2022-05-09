addpath('PPA_RA_Tests');
testName='Abbot PanBio';
 ts=4.3;
 VAC=0;
 VOC=0;
[MLE_RTPCR,MLE_Ag,MLE_PPA,Dt,totalpos,truepos,w,t,CCtestRTPCR,SymPRTPCR,CCtest,SymP] = Sensitivity_for_Plotting(testName,ts,VAC,VOC);
CCtest=CCtest(end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
subplot('Position',[0.255357142857143,0.204761904761905,0.719642857142857,0.742857142857143]);
p1=plot(t,MLE_PPA,'-','color',CCtest,'LineWidth',2); hold on;
p2=scatter(Dt(w==1),100.*truepos(w==1)./totalpos(w==1),40,SymP{1},'filled','MarkerEdgeColor',CCtest,'MarkerFaceColor',CCtest);
scatter(Dt(~isnan(w) & w<1),100.*truepos(~isnan(w) & w<1)./totalpos(~isnan(w) & w<1),40,SymP{1},'LineWidth',2,'MarkerEdgeColor',CCtest);

set(gca,'LineWidth',1.1,'tickdir','out','Fontsize',18,'XTick',[0:5:40],'xlim',[0 40],'XMinorTick','on','Yminortick','on','YTick',[0:20:100],'Ylim',[0 100]);
ytickformat('percentage')
xlabel('Days since symptom onset','Fontsize',22);
box off;
ylabel({'Percent','positive agreement'},'Fontsize',22);

rmpath('PPA_RA_Tests');
print(gcf,['PPA_PanBio.png'],'-dpng','-r600');