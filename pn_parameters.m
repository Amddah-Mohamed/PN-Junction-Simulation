% Paramèter de la capa MOS%
%-------------------------%


% Tension interne de la diode
Vbi = 0.56 ;%%% A COMPLETER %%% ;


% Epaisseur des cours de desertion a Vapp = 0
Wn = 0.492e-6; %%% A COMPLETER %%% ; % Epaisseur cote n
Wp = 0.492e-6;%%% A COMPLETER %%% ; % Epaisseur cote p
Wdep = 2*Wn ;%%% A COMPLETER %%% ;

% Warning
if(Wdep>Wtot)
    disp('attention : epaisseur de film inférieur à la couche de désertion max');
    pause
end



