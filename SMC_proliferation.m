% SMC_dynamic:

% Generate a random number that will be compared with the probability of
% mitosis/apoptosis. If such probability will be higher, the cell will
% divide/die, otherwise nothing will happen and the algorithm will switch
% to the next cell
test=rand(1); 

% General note: the probability densities associated to the different
% cellular events are function of the environmental conditions that are
% specific of the current SMC. Furthermore mitosis/apoptosis in intima
% should be function of shear stress, and at the same way the
% mitosis/apoptosis in media should be function of wall tension.

% Apoptosis (forget about this for now)
if intima_wall(j1,k1)==1, %inside the intima 
    p_apoptosis=probability_SMC(1)*activity;
else
    p_apoptosis=probability_SMC(1)*activity*(1+Ap(9)*sigma_r_local);
end

%% Mitosis: this formulation needs to be reviewed
if intima_wall(j1,k1)==1, %inside the intima 
    p_division=probability_SMC(2)*activity*(1-Ap(6)*min(tau_local,0));
else % in the media
    p_division=probability_SMC(2)*activity;
end

if test<p_apoptosis,
    change_cell(k)=-1; % SMC go to apoptosis.
else
    if (1-test)<p_division;
        change_cell(k)=1; % SMC go to mitosis.
    end
end

% If the cell has undergone mitosis, create a new SMC and place it close to
% the original one with random orientation
if change_cell(k)==1, 
    teta=rand(1)*2*pi; % pick random oriention on new cell location 
    
    % X2smc: coordinates of the daughter.
    % The phylosophy for SMC placement is that the daughter stays close to the mother
    X2smc(k,1)=Xsmc(k,1)+1.5*cos(teta)*Rs; % x coordinate
    X2smc(k,2)=Xsmc(k,2)+1.5*sin(teta)*Rs; % y coordinate
    X2smc(k,3)=Xsmc(k,3); % SMC is in the same layer of the mother
    
    XXsmc=X2smc(k,1:3);
    SMC_nearwall_correction % makes sure that the cell foes not go out of bound
    
    % Here X2smc contains the coordinates of the daughter cell (just in
    % case corrected) and in the third column just if the cell is in the
    % media(3) or in the intima(2)
    X2smc(k,1:2)=XXsmc(1,1:2); 
end


if change_cell(k)~=0, % manage mass balance
% create the source or sink term corresponding to mitosis/apopotosis 
    for i=1:N-1,
        for j=1:N-1,
            deltau1=0; 
            deltau2=0;
            % Cartesian distance between current site and the cell that
            % divided or died
            dx=abs(x(i)-Xsmc(k,1)); 
            dy=abs(x(j)-Xsmc(k,2));
            
            % If the current site is close enough to the cell thst has
            % divided, either died (close enough = in the range of the resolution of the grid)
            if dx<2*h,
                if dy<2*h,
                   deltau1=C2*(1.+cos(C4*dx)); 
                   deltau2=C3*(1.+cos(C5*dy));
                   % create an imbalance that is positive if the cell
                   % divided and negative if the cell died. The figurative
                   % idea is the following: the medium where the cells are
                   % immersed, works as a sort of bubble. IF a cell is
                   % generated, we will have a spike that will drive the
                   % matrix to re-organize around that spike. IF a cell
                   % dies, leaving a hole in the bubble, the matrix will
                   % re-organize accordingly.
                   S(i,j)=S(i,j)+Asmc*deltau1*deltau2*change_cell(k);
                end
            end
        end    
    end
end

