function [U]=poisson(Nx,Ny,UBL,RHS,L,UU,P)
% Nx dimension in x.
% Ny dimension in y.
% UBL for boundary conditions
% RHS for rhs of Helmotz problem.

U=UBL;
% construction of linear system:
% boundary condition:
b=RHS;
%b(:,1)=UBL(:,1); b(:,Ny+1)=UBL(:,Ny+1);
%b(1,:)=UBL(1,:); b(Nx+1,:)=UBL(Nx+1,:);

b=reshape(b,(Nx+1)*(Ny+1),1); b((Nx+1)*(Ny+1)+1,1)=1;
Y=L\(P*b);
U=UU\Y;
U=reshape(U,Nx+1,Ny+1);