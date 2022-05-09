%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency testing: RTPCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parpool(20); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW Sensitivity curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts=[6.3 5 4.3 3.2];
td=ts+20;
R0=[2.79 4.19 5.08 6.57];

SelfIsolate=1; %If sympmatics self-isolate
R0S=R0;
R0A=R0;


delayTR=1;


for VAC=0:1
    for VOC=0:3
        RTotA=zeros(14,1);
        RTotS=zeros(14,1);
        testtype=cell(14,1);
        timetoff=cell(14,1);
        NT=zeros(14,1);
        for dT=1:14
            timetoff{dT}=[0:dT:(floor(td(VOC+1)))];
            NT(dT)=length([timetoff{dT}]);
            temp=cell(NT(dT),1);
            testtype{dT}=temp;
        end
        parfor dT=1:14  

            [RS,RA] = SerialTestingDelay(testtype{dT},[timetoff{dT}],R0S(VOC+1),R0A(VOC+1),ts(VOC+1),td(VOC+1),SelfIsolate,NT(dT),dT,VAC,VOC,delayTR);

            RTotS(dT)=sum(RS);
            RTotA(dT)=sum(RA);
        end
        save(['VOC=' num2str(VOC) '-day_24h_Delay_Testing_Frequency_RTPCR_Vaccinated=' num2str(VAC) '.mat']);
    end
end
