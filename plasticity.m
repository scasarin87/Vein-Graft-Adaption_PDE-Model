%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Tissue Plasticity (TP)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time allowed to the structure to apply to plasiticity moduel and minimize
% the energy
maxT_plasticity=5;

% unknowns:
u=zeros(N,N-1); 
v=zeros(N-1,N); 
p=zeros(N-1,N-1);

hist_Aire = zeros(1,maxT_plasticity); % initialize the vector to keep track of the wall area variation
                                    % during the relaxation of the
                                    % strucutre
for kt=1:maxT_plasticity, % repeat the relaxation for N times, where N is defined
                          % at the beginning of the code
    
    oldXa=Xa; % coordinates of EEL
    oldXb=Xb; % coordinates of lumen border
    
    % Get the mean distance between the 360 points belonging to the lumen
    % border(Hb), EEL(Ha), and IEL(Hab)
    distancea=zeros(Ma,1); % Initialize the vector of distances
    for k2=1:Ma-1,
        distancea(k2)=sqrt((Xa(k2,1)-Xa(k2+1,1))^2+(Xa(k2,2)-Xa(k2+1,2))^2);
    end
    distancea(k2)=sqrt((Xa(1,1)-Xa(Ma,1))^2+(Xa(1,2)-Xa(Ma,2))^2);
    Ha=mean(distancea); % mean distance between coordinates of EEL
    
    distanceb=zeros(Mb,1);
    for k2=1:Mb-1,
        distanceb(k2)=sqrt((Xb(k2,1)-Xb(k2+1,1))^2+(Xb(k2,2)-Xb(k2+1,2))^2);
    end
    distanceb(k2)=sqrt((Xb(1,1)-Xb(Mb,1))^2+(Xb(1,2)-Xb(Mb,2))^2);
    Hb=mean(distanceb); % mean distance between coordinates of lumen border
    
    distanceab=zeros(Mab,1);
    for k2=1:Mab-1,
        distanceab(k2)=sqrt((Xab(k2,1)-Xab(k2+1,1))^2+(Xab(k2,2)-Xab(k2+1,2))^2);
    end
    distanceab(k2)=sqrt((Xab(1,1)-Xab(Mab,1))^2+(Xab(1,2)-Xab(Mab,2))^2);
    Hab=mean(distanceab); % mean distance between coordinates of IEL
    clear distancea distanceb distanceab
    
    Aire=0; % Initialize the area of the wall
    
    % Awall = Agraft - Alumen
    for k1=1:Ma-1, 
        Aire=Aire+0.5*abs((Xa(k1,1)-0.5)*(Xa(k1+1,2)-0.5)-(Xa(k1+1,1)-0.5)*(Xa(k1,2)-0.5));
    end
    for k1=1:Mb-1,
        Aire=Aire-0.5*abs((Xb(k1,1)-0.5)*(Xb(k1+1,2)-0.5)-(Xb(k1+1,1)-0.5)*(Xb(k1,2)-0.5));
    end
    hist_Aire(kt)=Aire; % keep track of the wall area variation during relaxation
    
    %% Evaluation of V(*)
    if semi_imp==0,
        % Fully explicit 
        prediction_fully_explicit
    else    
        % Semi-implicit
        prediction_semi_implicit
    end
    
    %% Evaluation of V(n+1), P(n+1)
    correction
    
    % Evaluation of X(n+1)
    move_membrane
    
    % Move SMC according to the wall displacement
    move_smc
    
    %% Evaluation of f(n+1)
    force_eval
    %% Evaluation of F(n+1)
    force_distrib
   
    NormXt=(norm(Xa-oldXa)/Ma+norm(Xb-oldXb)/Mb)/dt;
    if kt==1,
        NormXt0=NormXt;
    end
end
 if print_cfl==1,
    cfl=dt*max(max(max(u)),max(max(v)))/h;
end

clear hist_Aire
Interface_polar_coordinate

