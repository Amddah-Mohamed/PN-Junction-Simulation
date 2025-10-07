%Epaisseur totale de la jonction
Wtot = 2e-6 ; %[m]

% Dopage di Si
Na = 5e20 ; %[/m3] - type p
Nd = 5e20 ; %[/m3] - type n
ni = 1e16 ; %[/m3]

% Nombre de points de simulation (definition utilisateur)
Nx=1000;

% deltaX dans tout le futur maillage
dx = Wtot/Nx  ; 

% Tension appliquee %
Vapp_plus = 0:0.1:1.0;
Vapp_moins = -1.0:0.05:0.1;
Vapp=[Vapp_moins,Vapp_plus];
% Critère de convergence

Tau_list=[1e-3,1e-6];
SetErr =1e-8 ;