function V=Infectivity_Profile(t,td,ts)
    V=zeros(size(t));
    m=-4;
    s=1.85;
    a=5.85;
    tau=5.42;

    rt=tau/ts;


    t_temp=t-ts;

    V(t_temp>=0)=exp(-(t_temp(t_temp>=0)-m)./s)./((1+exp(-(t_temp(t_temp>=0)-m)./s)).^(a+1));

    V(t_temp<0)=exp(-(t_temp(t_temp<0).*rt-m)./s)./((1+exp(-(t_temp(t_temp<0).*rt-m)./s)).^(a+1));

    V(t>td)=0;

end