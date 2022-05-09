function V=Transform_Vac_VOC(t,VAC,VOC,ts,td)

if(VAC==0) && (VOC==0)
    V=Infectivity_Profile(t,td,ts);
else
    
    
    V=zeros(size(t)); 
    load('Peak_Infection.mat','mmv','tsv');
    mm=pchip(tsv,mmv,ts);
    if(VAC==1) && (VOC==0)   
       % Prior to the peak
       m_nVAC=(40-20.7)/(-3.5-0);
       m_VAC=(40-20.5)/(-3.2-0);
       xxp=t(t<=mm)-mm;
       y_VAC=m_VAC.*xxp; % The average of the peak
       t_transform=(y_VAC)./m_nVAC+mm; % The average of the peak
       
       V(t<=mm)=Infectivity_Profile(t_transform,td,ts);
       % After the peak

       m_nVAC=(20.7-40)/(0-7.5);
       m_VAC=(20.5-40)/(0-5.5);
       xxa=t(t>mm)-mm;
       y_VAC=m_VAC.*xxa; % The average of the peak   
       t_transform=(y_VAC)./m_nVAC+mm; % The average of the peak
       V(t>mm)=Infectivity_Profile(t_transform,td,ts);

    elseif(VAC==1)
       y_VAC=zeros(size(t));   
       y_nVAC=zeros(size(t));
       % Prior to the peak
       m_nVAC=(40-20.7)/(-3.5-0);
       m_VAC=(40-20.5)/(-3.2-0);
       xxp=t(t<=mm)-mm;
       y_nVAC(t<=mm)=m_nVAC; % The average of the peak
       y_VAC(t<=mm)=m_VAC; % The average of the peak

       % After the peak

       m_nVAC=(20.7-40)/(0-7.5);
       m_VAC=(20.5-40)/(0-5.5);
       xxa=t(t>mm)-mm;

       y_nVAC(t>mm)=m_nVAC; % The average of the peak   
       y_VAC(t>mm)=m_VAC; % The average of the peak   
       
    elseif(VAC==0)   
       y_VAC=ones(size(t));   
       y_nVAC=ones(size(t)); 
       
    end

    if(VOC==1)

        % Prior to the peak
       m_nVOC=(40-20.1)/(-4.2-0);
       m_VOC=(40-21.0)/(-3.4-0);
       xxp=t(t<=mm)-mm;

       y_VOC=m_VOC.*xxp; % The average of the peak

       Ct=y_VOC.*y_VAC(t<=mm)./y_nVAC(t<=mm);
       t_transform=(Ct)./m_nVOC+mm; % The average of the peak
       
       
       V(t<=mm)=Infectivity_Profile(t_transform,td,ts);
       
       % After the peak

       m_nVOC=(20.1-40)/(0-7.3);
       m_VOC=(21.0-40)/(0-6.2);
       xxa=t(t>mm)-mm;

       y_VOC=m_VOC.*xxa; % The average of the peak

       Ct=y_VOC.*y_VAC(t>mm)./y_nVAC(t>mm);
       t_transform=(Ct)./m_nVOC+mm; % The average of the peak
       
       V(t>mm)=Infectivity_Profile(t_transform,td,ts);
       

    elseif(VOC==2)

       % Prior to the peak
       m_nVOC=(40-20.1)/(-4.2-0);
       m_VOC=(40-19.8)/(-3.0-0);
       xxp=t(t<=mm)-mm;

       y_VOC=m_VOC.*xxp; % The average of the peak

       Ct=y_VOC.*y_VAC(t<=mm)./y_nVAC(t<=mm);
       t_transform=(Ct)./m_nVOC+mm; % The average of the peak
       
       V(t<=mm)=Infectivity_Profile(t_transform,td,ts);
       
       % After the peak

       m_nVOC=(20.1-40)/(0-7.3);
       m_VOC=(19.8-40)/(0-6.2);
       xxa=t(t>mm)-mm;

       y_VOC=m_VOC.*xxa; % The average of the peak

       Ct=y_VOC.*y_VAC(t>mm)./y_nVAC(t>mm);
       t_transform=(Ct)./m_nVOC+mm; % The average of the peak
       
       
       V(t>mm)=Infectivity_Profile(t_transform,td,ts);

    elseif(VOC==3)
      % Delta First
       VD=zeros(1,2001);
       
       % Prior to the peak
       m_nVOC=(40-20.1)/(-4.2-0);
       m_VOC=(40-19.8)/(-3.0-0);
       xxp=linspace(-mm-20,0,1001);
       y_VOC=m_VOC.*xxp; % The average of the peak

       Ct=y_VOC;
       t_transform=(Ct)./m_nVOC+mm; % The average of the peak
       
       VD(1:1001)=Infectivity_Profile(t_transform,td,ts);
       
       % After the peak

       m_nVOC=(20.1-40)/(0-7.3);
       m_VOC=(19.8-40)/(0-6.2);
       xxa=linspace(0,ts+40,1001);
       xxa=xxa(2:end); 
       
       t_temp=[xxp xxa]+mm;
       y_VOC=m_VOC.*xxa; % The average of the peak

       Ct=y_VOC;
       t_transform=(Ct)./m_nVOC+mm; % The average of the peak
       
       
       VD(1002:2001)=Infectivity_Profile(t_transform,td,ts);
       
       
      % Prior to the peak
       m_delta=(40-20.5)/(-4.67-0);
       m_omicron=(40-23.3)/(-4.52-0);

       xxp=t(t<=mm)-mm;
       
       y_omicron=m_omicron.*xxp; % The average of the peak


       Ct=y_omicron.*y_VAC(t<=mm)./y_nVAC(t<=mm);
       t_transform=(Ct)./m_delta+mm; % The average of the peak
       
       
       V(t<=mm)=pchip(t_temp(t_temp<=mm),VD(t_temp<=mm),t_transform);

       % After the peak
       m_delta=(20.5-40)/(0-6.23);
       m_omicron=(23.3-40)/(0-5.35);

       xxa=t(t>mm)-mm;
       y_omicron=m_omicron.*xxa; % The average of the peak


       Ct=y_omicron.*y_VAC(t>mm)./y_nVAC(t>mm);
       t_transform=(Ct)./m_delta+mm; % The average of the peak
       
       V(t>mm)=pchip(t_temp(t_temp>mm),VD(t_temp>mm),t_transform);
    end
end
end