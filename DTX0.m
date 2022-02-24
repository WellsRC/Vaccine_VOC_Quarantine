% Random entry in qiaratine with testing on exit
clear;

pobj=parpool(20); % Parallel pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT-PCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=[0:14]; % Quarantine durations consideredd
SelfIsolate=1; % Self-isolation

IncP=[3 4.3 5 6.3];
R0=[6.93 5.08 4.19 2.79];
VOC=[3 2 1 0];
% Allcoate memory for output
IDSLS=zeros(length(IncP),length(q)); 
IDSLA=zeros(length(IncP),length(q)); 

[IncPv,qv]=meshgrid(IncP,q);
[R0v,~]=meshgrid(R0,q);

[VOCv,~]=meshgrid(VOC,q);


IDSLS=IDSLS(:); % Vectorize the matrix
IDSLA=IDSLA(:); % Vectorize the matrix
IncPv=IncPv(:);
qv=qv(:);
R0S=R0v(:); % Set R0 for symptomatic
R0A=R0v(:); % Set R0 for asymptomatic
VOCv=VOCv(:);

testtype=cell(1,1);
testtype{1}=[];

PT24=1;
for VAC=0:1
    parfor jj=1:60

        if(qv(jj)>0)
            IDSLS(jj)=((1./IncPv(jj)).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,max(qv(jj)-1,0),testtype,R0S(jj),R0A(jj),0,IncPv(jj),(IncPv(jj)+20),SelfIsolate,VAC,VOCv(jj)),0,IncPv(jj),qv(jj),inf));
            IDSLA(jj)=((1./(IncPv(jj)+20)).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t+u,u,max(qv(jj)-1,0),testtype,R0S(jj),R0A(jj),1,IncPv(jj),(IncPv(jj)+20),0,VAC,VOCv(jj)),0,(IncPv(jj)+20),qv(jj),inf));  
        else
            % Infected after the test was adminsted 
            Posttest=((1./PT24).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,[],testtype,R0S(jj),R0A(jj),0,IncPv(jj),(IncPv(jj)+20),SelfIsolate,VAC,VOCv(jj)),0,PT24,@(u)(PT24-u),inf));
            Pretest=((1./(IncPv(jj)-PT24)).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,0,testtype,R0S(jj),R0A(jj),0,IncPv(jj),(IncPv(jj)+20),SelfIsolate,VAC,VOCv(jj)),0,IncPv(jj)-PT24,@(u)(u+PT24),inf));
            IDSLS(jj)=(PT24./IncPv(jj)).*Posttest+((IncPv(jj)-PT24)./IncPv(jj)).*Pretest;

            Posttest=((1./PT24).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,[],testtype,R0S(jj),R0A(jj),1,IncPv(jj),(IncPv(jj)+20),0,VAC,VOCv(jj)),0,PT24,@(u)(PT24-u),inf));
            Pretest=((1./((IncPv(jj)+20)-PT24)).*integral2(@(u,t)InfectiousnessfromInfectionTesting(t,u,0,testtype,R0S(jj),R0A(jj),1,IncPv(jj),(IncPv(jj)+20),0,VAC,VOCv(jj)),0,(IncPv(jj)+20)-PT24,@(u)(u+PT24),inf));
            IDSLA(jj)=(PT24./(IncPv(jj)+20)).*Posttest+(((IncPv(jj)+20)-PT24)./(IncPv(jj)+20)).*Pretest;            
        end
    end

save(['TestingonExit_RTPCR_24hrDelay_Vaccinated=' num2str(VAC) '.mat']);
end