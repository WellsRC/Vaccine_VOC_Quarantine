function [RS,RA] = SerialTestingDelay(testtype,timetoff,R0S,R0A,ts,td,SelfIsolate,NT,dT,VAC,VOC,delayTR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


RA=zeros(length(timetoff)+1,1);
RS=zeros(length(timetoff)+1,1);

RS(1)=integral2(@(u,tx)InfectiousnessfromInfectionTesting(tx,u,[],testtype,R0S,R0A,0,ts,td,SelfIsolate,VAC,VOC),0,dT,0,@(u)(u+delayTR))./dT;
RA(1)=integral2(@(u,tx)InfectiousnessfromInfectionTesting(tx,u,[],testtype,R0S,R0A,1,ts,td,0,VAC,VOC),0,dT,0,@(u)(u+delayTR))./dT;

parfor ii=0:(NT-2)
    dayT=[timetoff(1:(ii+1))]; % have to add +1 to ii as ii starts at zero
    % have to add 2 to ii for RS and RA as the ii starts at zero and
    % need index to start at 2
    RS(ii+2)=integral2(@(u,tx)InfectiousnessfromInfectionTesting(tx,u,dayT,testtype,R0S,R0A,0,ts,td,SelfIsolate,VAC,VOC),0,dT,@(u)(u+dT.*(ii)+delayTR),@(u)(u+dT.*(ii+1)+delayTR))./dT;
    RA(ii+2)=integral2(@(u,tx)InfectiousnessfromInfectionTesting(tx,u,dayT,testtype,R0S,R0A,1,ts,td,0,VAC,VOC),0,dT,@(u)(u+dT.*(ii)+delayTR),@(u)(u+dT.*(ii+1)+delayTR))./dT;
end

dayT=[timetoff];
RS(end)=integral2(@(u,tx)InfectiousnessfromInfectionTesting(tx,u,dayT,testtype,R0S,R0A,0,ts,td,SelfIsolate,VAC,VOC),0,dT,@(u)(u+dT.*(NT-1)+delayTR),td)./dT;
RA(end)=integral2(@(u,tx)InfectiousnessfromInfectionTesting(tx,u,dayT,testtype,R0S,R0A,1,ts,td,0,VAC,VOC),0,dT,@(u)(u+dT.*(NT-1)+delayTR),td)./dT;
        
end

