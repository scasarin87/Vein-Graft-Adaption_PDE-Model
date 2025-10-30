% meca_solver: compute the mechanical environment conditions

% Parabolic profile of velocity
if mod(hour,24)==0, % Re-evaluate it every day (~24hours)
    h2=1/h/h;
    % calcul du residus:
    nsa=zeros(N,N); residus=zeros(N,N);
    for ks=1:N*N,
        for j=2:N-1,
            for k=2:N-1,
                if lumen_hole(j,k)==1,
                    res=h2*(-nsa(j,k)+1/4*(nsa(j+1,k)+nsa(j-1,k)+nsa(j,k+1)+nsa(j,k-1)))-1;
                    nsa(j,k)=nsa(j,k)+1/h2/5*res;
                end
            end
        end
    end
    v_ns=-1/max(max(abs(nsa))).*nsa; % normalization of the basic solution
    clear LF UF PF Y RHS
    

    % Flow
    flux1=sum(sum(v_ns));
    v_ns=Flux0/flux1.*v_ns;
    
    % new we computer the normal direction at the wound edge using
    % a level set technic:
    wa=zeros(N,N); % wounded area
    wa=lumen_hole; tau=wa;
    newwa=zeros(N,N);
    for ks=1:round(max(N,N)),
        for j=2:N-1,
            for k=2:N-1,
                if lumen_hole(j,k)==1,
                    delta=wa(j+1,k)+wa(j-1,k)+wa(j,k+1)+wa(j,k-1);
                    newwa(j,k)=(0.5*delta/4+0.5*wa(j,k));
                end
            end
        end
        wa(2:N-1,2:N-1)=newwa(2:N-1,2:N-1);
    end
    clear newwa
    
    % get normal direction and gradient of flow:
    for j=2:N-1,
        for k=2:N-1,
            wax=wa(j+1,k)-wa(j-1,k);
            way=wa(j,k+1)-wa(j,k-1);
            delta=sqrt(wax^2+way^2);
            if delta>0,
                wax=wax/delta; way=way/delta;
            end
            v_nsx=v_ns(j+1,k)-v_ns(j-1,k);
            v_nsy=v_ns(j,k+1)-v_ns(j,k-1);
            tau(j,k)=wax*v_nsx+way*v_nsy;
        end
    end
    clear ks
    tau=tau.*lumen_hole; % tau is zero inside the wall.
    oldtau=tau; % keep track of tau vaule for some times.
else
    tau=oldtau; % tau is recomputed only every 24 hours.
end



% from the shear sress at the wall, we derive the growth factor inside the
% wall, thanks to a diffusion process. Note that shear stress is relative
% to a target value. As a matter of fact, we made the hypothesis that metabolism of cell is sensitive to
% the deviation from the basic solution, i.e stable vein.

% construction of the boundary of the lumen_hole (inside the hole):
newlumen=zeros(N,N); wall=zeros(N,N);
for j=2:N-1,
    for k=2:N-1,
        if lumen_hole(j,k)==1,
            newlumen(j,k)=(lumen_hole(j+1,k)+lumen_hole(j-1,k)+lumen_hole(j,k+1)+lumen_hole(j,k-1))/4;
        end
    end
end
for j=2:N-1,
    for k=2:N-1,
        if newlumen(j,k)<1 & newlumen(j,k)>0,
            wall(j,k)=1;
        end
    end
end
clear newlumen

% construction of Tau at the wall:
Tau_wall=tau.*wall;

% computation of the shear stress mean at the wall.
tau_wall_mean=0; n_mean=0;
for j=2:N-1,
    for k=2:N-1,
        if Tau_wall(j,k)~=0,
            tau_wall_mean=tau_wall_mean+Tau_wall(j,k);
            n_mean=n_mean+1;
        end
    end
end
tau_wall_mean=tau_wall_mean/n_mean;

% tau:Shear stress in the sites belonging to membrane lumen/intima
tau=wall.*1./mean_tau_wall_0.*Tau_wall.*(tau_wall_mean-mean_target_0)/mean_tau_wall_0;


% ATTENTION: we look at the negative delta of shear stress versus its
% target now.
newtau=tau;
for ks=1:tau_depth,
    for j=2:N-1,
        for k=2:N-1,
            if lumen_hole(j,k)==0,
                delta=tau(j+1,k)+tau(j-1,k)+tau(j,k+1)+tau(j,k-1);
                newtau(j,k)=(0.5*delta/4+0.5*tau(j,k));
            end
        end
    end
    tau(2:N-1,2:N-1)=newtau(2:N-1,2:N-1);
end
clear newtau

%% Shear stress inside the wall
Tau=(tau.*(1-lumen_hole));
hhhh(hour)=max(max(abs(Tau))); % keep track of the max value of shear stress in the wall



%% Useful mechanical quantities to calculate later the strain energy
p1=P_internal; p2=P_external; % should be invariant in fact

% r1 = lumen radius; r2 = graft radius
cA=p1*r1^2/(r2^2-r1^2);
cB=p2*r2^2/(r2^2-r1^2);

% constant for the displacement assuming no dilatation in longitudinal direction:
C10=(cA-cB)/(lambda+mu);
C20=(p1-p2)/(r2^2-r1^2)*r1^2*r2^2/(2*mu);

u1=C10/2*r1+C20/r1;
u2=C10/2*r2+C20/r2;




