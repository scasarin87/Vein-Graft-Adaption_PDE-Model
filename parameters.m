% Parameters

analytical_shear=0; % No radial simmetry assumption! Otherwise it would be 1

%% Initial geometrical measures
Rlumen=1.0; % Initial lumen radius [mm]
REEL=1.5; % Initial External Elastic Lamina radius [mm]
R_SMC=0.04; % Radius of a SMC [mm]

% Physical Scale constant for normalization 
Cs=4; % 1 unit of ABM = 4mm

% Geometrical measured normalized on the dimension of a single unit (Cs=4mm)
Rb=Rlumen/Cs; % Lumen radius
Ra=REEL/Cs; % Normalized EEL radius
Rs=R_SMC/Cs; % Normalized radius of a SMC

N_SMC=3; % thickness of intima expressed in number of SMC (~3 cells)
Rab=Rlumen/Cs+N_SMC*Rs; % IEL radius already normalized


%% Potential field to drive cells motion
% Motion is led by the region of steepest descent + noise.
% in other words the force derive from a potential:
% f=-grad P
% displacement will be proportional to f, i.e pure viscous flow with no
% inertia, and instantateneous limit speed.

% Stiffness parameters associated to graft membranes:
% Ca0(stiffness of the EEL), Cab0(stiffness of the IEL), and Cb0(stiffness of the lumen)
% are proportional to the repulsion that each membrane exploits on a single
% cell. The IEL is porous (cells can pass through) IF the chemotaxis
% driving the cell is big enough to win the resistance of the IEL.
Ca0=1; % strenght EEL.
Cab0=0.5; % strenght IEL.
Cb0=1; % strenght Lumen.

% relative strenght of potential forces between cell repulsion, chemotaxis
% and membrane

amplitude(2)=0.1; % chemotaxis has a 1/10 gradient of the cell-cell repulsion force.
amplitude(3)=1.; % motility toward lower density of SMC

chemo_depth=300; % related to the slope of chemoattractant concentration expressed in number of SMC diameter.
coef_chemo=0; % 0 no chemotaxis, 1 yes chemotaxis.

% Relative velocity for cell motility
time_scale(1)=0.3; % expressed in cell(diameter)/hour

% physical parameter in the Navier-Stokes equation :
epsilon=1.;

% elasticity parameter of the moving boundary :
sigma=1000;
sigma_ab=200.;

% relaxation of the membrane force to zero when X go to steady to avoid internal stress:
relaxation_force=0;

% number of SMC:
% Proportion of SMC versus ECM in tissue
theta=0.25; % ECM occupy 1-theta the total area of the wall.
% Ncell=60;

% parameter vector of the ABM:
bruit=0.1; % factor of random walk in smooth cell motility

% Base line for SMC division/apoptosis and ECM production/degradation 

% ATTENTION
Ap(1)=0.1; % Base level of cell division/apoptosis probability;
Ap(2)=0.1/12; % probability of matrix production/ probability of matrix degradation 

 
Ap(3)=24*14; Ap(4)=24*60; % activity factor related to macrophage =exp((hour-A[3])^2/A[4]);

% probability for SMC to go through the internal elastic lamina (IEL) to the intima;
% = A(5)*activity*(1+A(12)*tau(IEL))*(1+A(13)*sigma_r(IEL));
Ap(5)=0; Ap(12)=0; Ap(13)=0;
% -------------------------------------------------------------------------

Ap(6)=0;  % promotion of cell division by shear stress next to the lumen
        % moderated by promotion of cell division globaly inversely proportionel to shear stress
        % at the wall:
        % (1+A[6]* tau(j))/(1+A[7]*tau(1)) in the intima;
        
Ap(7)=0;

Ap(8)=40; % rate of decay of shear stress effect inside the wall expressed in SMC cell dimension number.

% -------------------------------------------------------------------------

Ap(9)=0; % promotion of cell apoptosis by Tension: (1+A[9]*Tension) in the media
        % where tension is in fact the radial mechanical stress;



%%%%%%%%%%%%%%%%%


Ap(10)=0; % matrix degeneration regulated by shear stress: (1+A[10]*tau) in the intima;

Ap(11)=0; % matrix degeneration regulated by tension in the media: (1+A[11]* Tension);