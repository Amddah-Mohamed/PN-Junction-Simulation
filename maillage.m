function [X,Nx]=maillage(Wtot,dx)


% Maillage fixe, ne pas modifer !

X= (-Wtot/2):dx:(Wtot/2);

Nx=length(X) ;



end