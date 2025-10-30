% keep track of where the cells are and eventually correct its position
% NB: current cell's coordainte in XXsmc = |x|y|where|

% Absolute Position of the current cell (radius and angle)
r_smc=sqrt((XXsmc(1,1)-0.5)^2+(XXsmc(1,2)-0.5)^2);
Teta_smc=atan2(XXsmc(1,1)-0.5,XXsmc(1,2)-0.5);

% Relative radiuses respect lumen/intima/media border
r_a=interp1(Theta_a,R_a,Teta_smc);
r_b=interp1(Theta_b,R_b,Teta_smc);
r_ab=interp1(Theta_ab,R_ab,Teta_smc);

if r_smc>r_a, % if the cell is outside the external wall
    where_smc(k)=1; % we need to set the cell back inside the wall
    
    % We set it back setting it just an SMC radius (Rs) away from the EEL 
    XXsmc(1,1)=(r_a-Rs)*sin(Teta_smc)+0.5; 
    XXsmc(1,2)=(r_a-Rs)*cos(Teta_smc)+0.5;
    
else
    
    if r_smc>r_ab, % here the cell is in media
        where_smc(k)=2;
    else
        if r_smc>r_b, % here is in intima
            where_smc(k)=3;
        else  % if the cell is inside the lumen
            where_smc(k)=4; % we need to set the cell back inside the wall
            
            % We set it back setting it just an SMC radius (Rs) away from the lumen 
            XXsmc(1,1)=(r_b+Rs)*sin(Teta_smc)+0.5; 
            XXsmc(1,2)=(r_b+Rs)*cos(Teta_smc)+0.5;
        end
    end
end

% If we need to move the SMC back in intima or media the phylosophy is the
% same. We decided to set it one radius far from the border

if where_smc(k)==2 && XXsmc(1,3)==3, % ned to move the smooth cell from the media back to the intima:
    XXsmc(1,1)=(r_ab-Rs)*sin(Teta_smc)+0.5; 
    XXsmc(1,2)=(r_ab-Rs)*cos(Teta_smc)+0.5; 
end
if where_smc(k)==3 && XXsmc(1,3)==2, % ned to move the smooth cell from the intima back to the media:
    XXsmc(1,1)=(r_ab+Rs)*sin(Teta_smc)+0.5; 
    XXsmc(1,2)=(r_ab+Rs)*cos(Teta_smc)+0.5; 
end
    
    
clear r_a r_b r_ab r_smc Teta_smc