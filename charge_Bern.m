%Tjarriba dyal Bernoulli 

function [n,p,rho]=charge_Bern(V_sol,X)
    
    material_properties;
    physical_constants;
    simulation_parameters;


    n=zeros(size(X));p=zeros(size(X));
    
    N=length(X);
    V=(1/Vt)*V_sol;
    n(1) = (-Na+sqrt(Na^2+4*ni^2))/2; % Left boundary (p-side)
    n(2)=n(1);
    n(N) = (Nd+sqrt(Na^2+4*ni^2))/2;         % Right boundary (n-side)
    p(1) = (Na+sqrt(Na^2+4*ni^2))/2;% Left boundary (p-side)
    p(2)=p(1);
    p(N) = (-Nd+sqrt(Na^2+4*ni^2))/2;  % Right boundary (n-side)
    
    for i=2:N-1
        n(i+1)=(1/Bern(V(i+1)-V(i)))*(n(i-1)*(Bern(V(i)-V(i+1))+Bern(V(i)-V(i-1)))-n(i-1)*Bern(V(i-1)-V(i)));
    
        p(i+1)=(1/Bern(-V(i+1)+V(i)))*(p(i-1)*(Bern(-V(i)+V(i+1))+Bern(-V(i)+V(i-1)))-p(i-1)*Bern(-V(i-1)+V(i)));
    
    end
    
    
    left_indices = (X <= 0);
    right_indices = (X > 0);
    dopage=zeros(size(X));
    dopage(right_indices)=Nd;
    dopage(left_indices)=-Na;
    
    % mise a jour de la densite de charges
    rho = q*(p-n+dopage);


end