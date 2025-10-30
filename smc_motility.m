% smc motility

% first get the influence of other smc:
% this potential express cell-cell interaction.
% we set just a localized repelant force to start:

% cells move of 1/3 of a space step at the time (less than cell's size)
deltat=h/3;  

% Re-initialize the matrix of SMCs coordinates
newXsmc=Xsmc;

% Randomize the wya we consider the various cells
nouvel_ordre=rand(Ncell,1); 
[~, newK]=sort(nouvel_ordre); clear nouvel_ordre ordre

% Initialize the matrix that keeps track of where cells are
where_smc=zeros(Ncell,1);  

for kr=1:Ncell, % Explore all the cells in the wall
    
    % k is the index of a random cell
    k=newK(kr); 
    
    FT=zeros(N,N);
    FT=coef_chemo*FC; % Motion by chemotaxis.
    
    i1=floor((Xsmc(k,1))/h)+1; 
    j1=floor((Xsmc(k,2))/h)+1;
    
    % Define the potential driven by chemotaxis (not using it right now)
    FCx=zeros(2,2); FCy=zeros(2,2);
    FCx(1,1)=(FT(i1+1,j1)-FT(i1-1,j1))/2/h; 
    FCy(1,1)=(FT(i1,j1+1)-FT(i1,j1-1))/2/h;
    FCx(2,1)=(FT(i1+2,j1)-FT(i1,j1))/2/h; FCy(2,1)=(FT(i1+1,j1+1)-FT(i1+1,j1-1))/2/h;
    FCx(1,2)=(FT(i1+1,j1+1)-FT(i1-1,j1+1))/2/h; FCy(1,2)=(FT(i1,j1+2)-FT(i1,j1))/2/h;
    FCx(2,2)=(FT(i1+2,j1+1)-FT(i1,j1+1))/2/h; FCy(2,2)=(FT(i1+1,j1+2)-FT(i1+1,j1))/2/h;
    Fx=deltat.*FCx; Fy=deltat.*FCy;
    
    % Gradient of the potential obtained thorugh linear combination
    ax=(x1(i1+1)-Xsmc(k,1))/h; bx=(Xsmc(k,1)-x1(i1))/h;
    ay=(x1(j1+1)-Xsmc(k,2))/h; by=(Xsmc(k,2)-x1(j1))/h;
    fx=ay*(ax*Fx(1,1)+bx*Fx(2,1))+by*(ax*Fx(1,2)+bx*Fx(2,2)); % potential toward x coordinate
    fy=ay*(ax*Fy(1,1)+bx*Fy(2,1))+by*(ax*Fy(1,2)+bx*Fy(2,2)); % potential toward y coordinate
    
    % Move the cell according to chemotaxis + natural noise given by the
    % particles crawling suspended in the medium
    newXsmc(k,1)=Xsmc(k,1)+deltat/h*bruit*(rand(1)-0.5)-deltat*fx; % make sure the noise is centered.
    newXsmc(k,2)=Xsmc(k,2)+deltat/h*bruit*(rand(1)-0.5)-deltat*fy;
    
    % Correct the position of the cell to make sure it stays in the right
    % layer
    
    % keep track of where the cells are:
    r_smc=sqrt((Xsmc(k,1)-0.5)^2+(Xsmc(k,2)-0.5)^2); % ditance of the cell from the center of the graft
    Teta_smc=atan2(Xsmc(k,1)-0.5,Xsmc(k,2)-0.5); % angle of the cell 
    r_a=interp1(Theta_a,R_a,Teta_smc); % correspondent graft radius
    r_b=interp1(Theta_b,R_b,Teta_smc); % correspondent lumen border radius
    r_ab=interp1(Theta_ab,R_ab,Teta_smc); % correspondent IEL radius
    
    if r_smc>r_a, % if cell ended up outside the wall ...
        where_smc(k)=1; % we need to set the cell back inside the wall
        Xsmc(k,1)=(r_a-Rs)*sin(Teta_smc)+0.5; Xsmc(k,2)=(r_a-Rs)*cos(Teta_smc)+0.5; newXsmc(k,:)=Xsmc(k,:);
    else
        if r_smc>r_ab, % if cell is beyond IEL
            where_smc(k)=2;
        else
            if r_smc>r_b, % if cell is beyond the lumen wall
                where_smc(k)=3;
            else
                % last chance ... the cell ended up inside the lumen
                where_smc(k)=4; % and obviously it requires correction
                % need to set the cell back inside the wall (in the intima)
                Xsmc(k,1)=(r_b+Rs)*sin(Teta_smc)+0.5; Xsmc(k,2)=(r_b+Rs)*cos(Teta_smc)+0.5; newXsmc(k,:)=Xsmc(k,:);
            end
        end
    end
    
    % If the cell leaves its domain, correct the position
    r_smc=sqrt((newXsmc(k,1)-0.5)^2+(newXsmc(k,2)-0.5)^2);
    Teta_smc=atan2(newXsmc(k,1)-0.5,newXsmc(k,2)-0.5);
    r_a=interp1(Theta_a,R_a,Teta_smc);
    r_b=interp1(Theta_b,R_b,Teta_smc);
    r_ab=interp1(Theta_ab,R_ab,Teta_smc);
        
    if where_smc(k)==2, 
        if r_smc<r_ab| r_smc>r_a, % correction needed if the SMC leaves its tissue layer
            newXsmc(k,:)=Xsmc(k,:); 
        end
    end
    if where_smc(k)==3,
        if r_smc>r_ab | r_smc<r_b, % correction needed if the SMC leaves its tissue layer
            newXsmc(k,:)=Xsmc(k,:); 
        end
    end
            
    clear Fx Fy ax ay bx by fx fy xx yy
end
Xsmc=newXsmc;

% motion to invade ECM cluster in the intima:

analysis_space_distribution
newXsmc=Xsmc;

for k=1:Ncell,
    if Xsmc(k,3)==3,
        i1=floor((Xsmc(k,1))/h)+1; j1=floor((Xsmc(k,2))/h)+1;
        % computation of the distance to ECM:
        distance_ECM=1;
        for i=1:N-1,
            for j=1:N-1,
                if ECM_intima(i,j)==1,
                    test=sqrt((x(i)-Xsmc(k,1))^2+(x(j)-Xsmc(k,2))^2);
                    if test<distance_ECM,
                        distance_ECM=test; xnext=x(i);ynext=x(j);
                    end
                end
            end
        end
        if distance_ECM<=(Ndiffusif+1)*h,
            normalisation=sqrt((xnext-Xsmc(k,1))^2+(ynext-Xsmc(k,2))^2);
            newXsmc(k,1)=Xsmc(k,1)+deltat*(xnext-Xsmc(k,1))/normalisation; newXsmc(k,2)=Xsmc(k,2)+h*(ynext-Xsmc(k,2))/normalisation;
        end 
    end
    if Xsmc(k,3)==2,
        i1=floor((Xsmc(k,1))/h)+1; j1=floor((Xsmc(k,2))/h)+1;
        % computation of the distance to ECM:
        distance_ECM=1;
        for i=1:N-1,
            for j=1:N-1,
                if ECM_media(i,j)==1,
                    test=sqrt((x(i)-Xsmc(k,1))^2+(x(j)-Xsmc(k,2))^2);
                    if test<distance_ECM,
                        distance_ECM=test; xnext=x(i);ynext=x(j);
                    end
                end
            end
        end
        if distance_ECM<=(Ndiffusif+1)*h,
            normalisation=sqrt((xnext-Xsmc(k,1))^2+(ynext-Xsmc(k,2))^2);
            newXsmc(k,1)=Xsmc(k,1)+deltat*(xnext-Xsmc(k,1))/normalisation; newXsmc(k,2)=Xsmc(k,2)+h*(ynext-Xsmc(k,2))/normalisation;
        end 
    end
end
oldXsmc=Xsmc;
Xsmc=newXsmc;
        
    


