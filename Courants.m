%CAlcul du courant de derive diffusion

function [XJ,J,Jn,Jp]=Courants(V,X,n,p)


    material_properties;
    physical_constants;
    simulation_parameters;

    XJ=X+(dx/2)*ones(size(X)); % creation d'un nouveau maillage pour les courants
    XJ(end)=[];

    Jn=zeros(size(XJ));
    Jp=zeros(size(XJ));
    for i=1:length(XJ)
        Jn(i)=(q*muSi*kT/dx)*(n(i+1)*B((V(i+1)-V(i))/Vt)-n(i)*B((V(i)-V(i+1))/Vt));
        Jp(i)=(-q*muSi*kT/dx)*(p(i+1)*B(-(V(i+1)-V(i))/Vt)-p(i)*B(-(V(i)-V(i+1))/Vt));
    end

    J=Jn+Jp; %densite de courant totale 

end