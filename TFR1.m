%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency testing: BD Veritor 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parpool(20); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW Sensitivity curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('RAgTest_Name.mat','testName');
NumTests=length(testName);


tsv=[6.3 5 4.3 3.2];
tdv=tsv+20;
R0=[2.79 4.19 5.08 6.57];


for TestN=1:NumTests
    for VAC=0:1
        for VOC=0:3
            ts=tsv(VOC+1);
            td=tdv(VOC+1);
            SelfIsolate=1; %If sympmatics self-isolate
            R0S=R0(VOC+1);
            R0A=R0(VOC+1);
            RTotA=zeros(14,1);
            RTotS=zeros(14,1);


            [~,betaAg]=ParameterCOVIDTest(testName{TestN},1);

            testtype=cell(14,1);
            timetoff=cell(14,1);
            NT=zeros(14,1);
            for dT=1:14
                timetoff{dT}=[0:dT:(floor(td))];
                NT(dT)=length([timetoff{dT}]);
                temp=cell(NT(dT),1);

                for ii=1:length([timetoff{dT}])
                    temp{ii}=betaAg;
                end
                testtype{dT}=temp;
            end

            parfor dT=1:14  

                [RS,RA] = SerialTesting(testtype{dT},[timetoff{dT}],R0S,R0A,ts,td,SelfIsolate,NT(dT),dT,VAC,VOC);

                RTotS(dT)=sum(RS);
                RTotA(dT)=sum(RA);
            end
            save(['Testing_Frequency_' testName{TestN} '_VOC= ' num2str(VOC) '_VAC= ' num2str(VAC) '.mat']);
        end 
    end
end