function [MLE_RTPCR,MLE_Ag,MLE_PPA,Dt,totalpos,truepos,w,t,CCtestRTPCR,SymPRTPCR,CCtest,SymP] = Sensitivity_for_Plotting(testName,ts,VAC,VOC)

[~,betaAg]=ParameterCOVIDTest(testName,1);
t=linspace(0,40,1001);

MLE_Ag = TestSensitivity(t,ts,betaAg,VAC,VOC);     
    
MLE_RTPCR = TestSensitivity(t,ts,[],VAC,VOC);

MLE_PPA=100.*LR(t,betaAg); % time zero is the time of symptom onset, which is what we want to be plotting from. 
    
load([testName '_LR_Parameters.mat'],'Dt','totalpos','truepos','w')


[CCtestRTPCR,SymPRTPCR]=ColourTests('RTPCR');
[CCtest,SymP]= ColourTests(testName); 
end

