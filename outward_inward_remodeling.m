%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Outward_Inward_Remodeling                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Outward_Inward_Remodeling: the module decides between inward/outward
% remodeling direction by minimizing the energy across the graft. Remember
% that matrix S in principle kept track of where the energy has been
% UNbalanced either by an SMC or by an ECM.

% We here use a mechanical model based on the circumferential simmetry
% assumption. We indeed want to minimize the energy of all the strucutre
% and not only locally.

Ksupporting=4000; % Stiffness of the supporting tissue.

r10=r1; % radius of the lumen
r20=r2; % radius of the graft
Aire0=pi*(r20^2-r10^2); % area of the wall (excluding the lumen)

v0=1; % inlet velocity supposed unitary
flux0=pi*r10*r10*v0; % v0 exact value does not matter as long as long as we keep the flux constant.

% Philosophy: what we need to decide is basically how to move the membranes
% in order to minimize the total energy of the structure. Whether we need
% to ove them inward or to move them outward of one step. In order to do
% that, we need to know:
% 1) the current energy recorded at r1=r10
% 2) the energy with r1plus = 1.0001*r1 that would move the membrane outward
% 3) the energy with r1minus = 0.0009*r1 that would move the membrane inward
% Once we know in which condition we will have the lower value of energy,
% we will decide whether to move inward or outward


% We compute the total mechanical energy for the current position r1;
% we will compute the new energy provided a 5 % variation inward or outward
% and decide on where remodeling go in order to pick the smallest energy.
rmoins=0.9999*r1;
rplus=1.0001*r1;

%% 1) Energy at r1=r10
for j=1:Ncell,
    rj=r1+(r2-r1)*(j-1)/Ncell; % unloaded position
    sigma_r(j)=cA*(1-r2^2/rj^2)-cB*(1-r1^2/rj^2); % radial tension
    sigma_a(j)=cA*(1+r2^2/rj^2)-cB*(1+r1^2/rj^2); % axial tension
    Tension(j)=sqrt(sigma_r(j)^2+sigma_a(j)^2); % normal tension
end
Et0=sum(Tension.*Tension)/Ncell; % Energy recorded at r1=r10
clear rj sigma_r sigma_a Tension


%% 2) Energy at r1minus = 0.0009*r1 --> inward remodeling assumed
r1moins=rmoins;
r2moins=sqrt((pi*r1moins^2+Aire0)/pi); % accordingly would move the external membrane too

% Re-compute all the mechanical quantities across the wall
% (velocity/constants/pressures) with the new radii we are considering to
% assume inward remodeling 
vmoins=flux0/(pi*r1moins*r1moins); % is the max of the paraboloide for example, all dependence are linear. 
                                   % The constant factor is eliminated in the ratio of flux.
% assume Bernouilli: v^2/2+p/rho=Ct
p1moins=1000/133*(v0^2/2+P_internal*133/1000-vmoins^2/2); % (1 mm Hg = 133 Pa, Pascal is for USI, i.e kg/m/s, rho is 1000 kg per cube meter)
p2moins=P_external+Ksupporting*(r2moins-r20); % neglect displacement
    
% constant for the stress distribution in the wall
cAmoins=p1moins*r1moins^2/(r2moins^2-r1moins^2);
cBmoins=p2moins*r2moins^2/(r2moins^2-r1moins^2);

% Compute the energy with r=rmoins
for j=1:Ncell,
    rj=r1moins+(r2moins-r1moins)*(j-1)/Ncell; % unloaded position
    sigma_r(j)=cAmoins*(1-r2moins^2/rj^2)-cBmoins*(1-r1moins^2/rj^2); % radial
    sigma_a(j)=cAmoins*(1+r2moins^2/rj^2)-cBmoins*(1+r1moins^2/rj^2); % axial
    Tension(j)=sqrt(sigma_r(j)^2+sigma_a(j)^2); % normal
end
if p1moins<=0,
    display('problem p1 negative')
    pause
end
Etmoins=sum(Tension.*Tension)/Ncell; % Energy 
clear rj sigma_r sigma_a Tension


%% 3) Energy at r1plus = 1.0001*r1 --> outward remodeling assumed
r1plus=rplus;
r2plus=sqrt((pi*r1plus^2+Aire0)/pi);
% conservation of flux gives:
vplus=flux0/(pi*r1plus*r1plus); % is the max of the paraboloide for example, all dependence are linear. 
                        % The constant factor is eliminated in the ratio of flux.
% assume Bernouilli: v^2/2+p/rho=Ct
p1plus=1000/133*(v0^2/2+P_internal*133/1000-vplus^2/2); % (1 mm Hg = 133 Pa, Pascal is for USI, i.e kg/m/s)
p2plus=P_external+Ksupporting*(r2plus-r20); % neglect displacement
    
% constant for the stress distribution in the wall
cAplus=p1plus*r1plus^2/(r2plus^2-r1plus^2);
cBplus=p2plus*r2plus^2/(r2plus^2-r1plus^2);

for j=1:Ncell
    rj=r1plus+(r2plus-r1plus)*(j-1)/Ncell; % unloaded position
    sigma_r(j)=cAplus*(1-r2plus^2/rj^2)-cBplus*(1-r1plus^2/rj^2);
    sigma_a(j)=cAplus*(1+r2plus^2/rj^2)-cBplus*(1+r1plus^2/rj^2);
    Tension(j)=sqrt(sigma_r(j)^2+sigma_a(j)^2);
end
if p1plus<=0,
    display('problem p1 negative')
    pause
end
Etplus=sum(Tension.*Tension)/Ncell;
clear rj sigma_r sigma_a Tension

%% What he have here
% Et0 --> energy recordedin absence of remodeling
% Etmoins = energy recorded in presence of inward remodeling
% Etplus = energy recorded in presence of outward remodeling
% We will choose the direction to move toward basing on the minimum energy
% that the structure assumes

% EE --> Value of minimum energy
% EI --> Index of minimum energy in the vector of energies
[EE,EI]=min([Etmoins Et0 Etplus]);

% Store the energies recorded in the history matrix. It is a good indicator
% of the direction of the remodeling 
histoire(hour,8)=Etmoins; 
histoire(hour,9)=Et0; 
histoire(hour,10)=Etplus; 




