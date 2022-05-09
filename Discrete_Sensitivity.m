clear;
clc;
addpath('PPA_RA_Tests');
VAC=0;
VOCN={'Original','Alpha','Delta','Omicron'};
%Delta:  LogNormal(1.249, 0.649)

% VOC=2;
% ts=lognstat(1.249, 0.649);

%Omicron: LogNormal(0.99, 0.64)
% 
% VOC=3;
% ts=lognstat(0.99, 0.64);

%Other:  LogNormal(1.434, 0.661)
% 
VOC=0;
ts=lognstat(1.434, 0.661);

Ntest=4;
Day=[1:60]';
Integration_Start=Day-1;
Integration_End=Day;

Test_Name=cell(Ntest,1);

Sensitivity_Day=zeros(Ntest,length(Day));

ii=1;

Test_Name{ii}='RT-PCR';


for dd=1:length(Day)                    
    Sensitivity_Day(ii,dd)=integral(@(x)TestSensitivity(x,ts,[],VAC,VOC),Integration_Start(dd),Integration_End(dd));
end

ii=2;

Test_Name{ii}='Abbot PanBio';
    
    [~,betaAg]=ParameterCOVIDTest(Test_Name{ii},1);
    for dd=1:length(Day)
        Sensitivity_Day(ii,dd)=integral(@(x)TestSensitivity(x,ts,betaAg,VAC,VOC),Integration_Start(dd),Integration_End(dd));
    end

ii=3;

Test_Name{ii}='BTNX Rapid Response';

% https://www.medrxiv.org/content/10.1101/2022.01.18.22269426v1.full-text
AbPB=0.645;
BTNX=0.781;

Sensitivity_Day(ii,:)=Sensitivity_Day(2,:)+(BTNX-AbPB)/(1-AbPB).*(Sensitivity_Day(1,:)-Sensitivity_Day(2,:));

ii=4;

Test_Name{ii}='Artron';

% https://dam.abbott.com/en-gb/panbio/120007883-v1-Panbio-COVID-19-Ag-Nasal-AsymptomaticSe.pdf
AbPB=102/104;
%http://www.artronlab.com/products/IFU/A03-50-422%20COVID-19AgIFU.pdf
Artron=154/168;

Sensitivity_Day(ii,:)=(Artron./AbPB).*Sensitivity_Day(2,:);


T=table(Test_Name,Sensitivity_Day);

writetable(T,['Discrete_Sensitivity_' VOCN{VOC+1} '.csv']);

rmpath('PPA_RA_Tests');