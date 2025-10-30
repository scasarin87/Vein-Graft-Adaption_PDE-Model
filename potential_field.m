% potential_field: build the potential field that drives the re-arrangement
% of SMCs within the viscous flow (~matrix)

% Motion is determined by 3 independent components:
% 1) Gradient of concentration
% 2) Random walk
% 3) Repulsion from more concentrated neighborhoods

% potential_field due to a chemoattractant that would diffuse uniformly
% from the wall.

% wall(i,j) = 1 if wall(i,j) belongs to lumen border
% wall(i,j) = 0 if wall(i,j) does NOT belong to lumen border
wa=wall; clear chemo_wall;
chemo_wall=lumen_hole; % initialize chemo_wall to overlap the lumen
for i=1:N,
    for j=1:N,
        if wa(i,j)~=0,
            wa(i,j)=1; % this seems to be useless
            %chemo_wall(j,k)=1;
            chemo_wall(i,j)=1;
        end
    end
end
newwa=wa;

for ks=1:chemo_depth,
    for j=2:N-1,
        for k=2:N-1,
            if chemo_wall(j,k)==0,
                delta=wa(j+1,k)+wa(j-1,k)+wa(j,k+1)+wa(j,k-1);
                newwa(j,k)=(0.5*delta/4+0.5*wa(j,k));
            end
        end
    end
    % explicit psuedo homogeneous Neuman boundary conditions:
    newwa(1,:)=newwa(2,:); newwa(N,:)=newwa(N-1,:); newwa(:,1)=newwa(:,2); newwa(:,N)=newwa(:,N-1);
    wa=newwa;
end
FC=(1-wa).*(1-lumen_hole);

% however the gradent of FC is much lower than the gradient of other local potential forces
% due to cell-cell interaction and membrane. Those are short interaction,
% i.e are esentially smooth dirac of the order of few mesh points.
% we compute below a correction factor, to make sure that amplitude give
% the weight of chemotaxis in the right space scale:
FCx=zeros(N,N); FCy=FCx;
for j=2:N-1,
    for k=2:N-1,
        FCx(j,k)=(FC(j+1,k)-FC(j-1,k))/2/h;
        FCy(j,k)=(FC(j,k+1)-FC(j,k-1))/2/h;
    end
end
Factor_correctif=0.5*0.5*C4*0.5/max(max(max(abs(FCx))),max(max(abs(FCy))));
%amplitude(2)=0; % do not use this potential field
FC=amplitude(2)*Factor_correctif.*FC;


