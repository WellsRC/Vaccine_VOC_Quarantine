clear;
clc;
addpath('PPA_RA_Tests');
VOCN={'Original','Alpha','Delta','Omicron'};


Ntest=4;
Day=[1:60]';
Integration_Start=Day-1;
Integration_End=Day;
Incubation_Period=[1:28]';

Test_Name=cell(length(Incubation_Period),1);
Vac_Status=cell(length(Incubation_Period),1);
Variant=cell(length(Incubation_Period),1);
T=[];
for VAC=0:1
    if(VAC==0)
        Vac_Status(:)={'Unvaccinated'};
    else
        Vac_Status(:)={'Vaccinated'};
    end
    for VOCv=3:3
        switch VOCv
            case 1
                VOC=0;
                Variant(:)={'Original'};
            case 2
                VOC=2;
                Variant(:)={'Delta'};
            case 3
                VOC=3;
                Variant(:)={'Omicron'};
        end
        %Delta:  LogNormal(1.249, 0.649)

        % VOC=2;
        % ts=lognstat(1.249, 0.649);

        %Omicron: LogNormal(0.99, 0.64)
        % 
        % VOC=3;
        % ts=lognstat(0.99, 0.64);

        %Other:  LogNormal(1.434, 0.661)
        % 
        % ts=lognstat(1.434, 0.661);
    for test=1:4
        
        Sensitivity_Day=zeros(length(Incubation_Period),length(Day));
        switch test
            case 1
                Test_Name(:)={'RT-PCR'};
                for ii=1:length(Incubation_Period)
                    for dd=1:length(Day)   
                        Sensitivity_Day(ii,dd)=integral(@(x)TestSensitivity(x,Incubation_Period(ii),[],VAC,VOC),Integration_Start(dd),Integration_End(dd));
                    end
                end
            case 2
                Test_Name(:)={'Abbot PanBio'};
                [~,betaAg]=ParameterCOVIDTest(Test_Name{1},1);
                for ii=1:length(Incubation_Period)
                    for dd=1:length(Day)
                        Sensitivity_Day(ii,dd)=integral(@(x)TestSensitivity(x,Incubation_Period(ii),betaAg,VAC,VOC),Integration_Start(dd),Integration_End(dd));
                    end
                end
            case 3
                Test_Name(:)={'BTNX Rapid Response'};

                    % https://www.medrxiv.org/content/10.1101/2022.01.18.22269426v1.full-text
                    AbPB=0.645;
                    BTNX=0.781;
                    
                    [~,betaAg]=ParameterCOVIDTest('Abbot PanBio',1);
                    for ii=1:length(Incubation_Period)
                        for dd=1:length(Day)
                            Sensitivity_Day_Ag=integral(@(x)TestSensitivity(x,Incubation_Period(ii),betaAg,VAC,VOC),Integration_Start(dd),Integration_End(dd));
                            Sensitivity_Day_R=integral(@(x)TestSensitivity(x,Incubation_Period(ii),[],VAC,VOC),Integration_Start(dd),Integration_End(dd));
                            Sensitivity_Day(ii,dd)=Sensitivity_Day_Ag+(BTNX-AbPB)/(1-AbPB).*(Sensitivity_Day_R-Sensitivity_Day_Ag);
                        end
                    end
            case 4
                Test_Name(:)={'Artron'};
                % https://dam.abbott.com/en-gb/panbio/120007883-v1-Panbio-COVID-19-Ag-Nasal-AsymptomaticSe.pdf
                AbPB=102/104;
                %http://www.artronlab.com/products/IFU/A03-50-422%20COVID-19AgIFU.pdf
                Artron=154/168;
                
                
                [~,betaAg]=ParameterCOVIDTest('Abbot PanBio',1);
                for ii=1:length(Incubation_Period)
                    for dd=1:length(Day)
                        Sensitivity_Day_Ag=integral(@(x)TestSensitivity(x,Incubation_Period(ii),betaAg,VAC,VOC),Integration_Start(dd),Integration_End(dd));
                        Sensitivity_Day(ii,dd)=(Artron./AbPB).*Sensitivity_Day_Ag;
                    end
                end
        end
        if(isempty(T))
            T=table(Variant,Test_Name,Vac_Status,Incubation_Period,Sensitivity_Day);
        else
            T=[T;table(Variant,Test_Name,Vac_Status,Incubation_Period,Sensitivity_Day)];
        end
    end
    end
end



writetable(T,['Discrete_Sensitivity_Table.csv']);

rmpath('PPA_RA_Tests');