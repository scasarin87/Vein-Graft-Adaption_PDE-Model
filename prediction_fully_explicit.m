% prediction.m
% prediction on the 2D incompressible equations
% convection and diffusion semi-implicit
% Gauss-Seidel

U=zeros(N,N+1);
V=zeros(N+1,N);
U(2:N-1,2:N)=u(2:N-1,1:N-1);
U(2:N-1,1)=-U(2:N-1,2);
U(2:N-1,N+1)=-U(2:N-1,N);
V(2:N,2:N-1)=v(1:N-1,2:N-1);  
V(1,2:N-1)=-V(2,2:N-1); 
V(N+1,2:N-1)=-V(N,2:N-1);

VenU=zeros(N,N+1); 
UenV=zeros(N+1,N);
px=zeros(N,N+1); 
py=zeros(N+1,N);
convU=zeros(N,N+1);
convV=zeros(N+1,N);
tauU=zeros(N,N+1);
tauV=zeros(N+1,N);

% approximation of explicit convective term:
% centre=0 for upwind of order 1 
% centre=1 for second order central differences.
centre=1;

for i=2:N-1,
    for j=2:N,
        VenU(i,j)=0.25*(V(i,j)+V(i+1,j)+V(i,j-1)+V(i+1,j-1));
        if centre==0,
              convU(i,j)=max(U(i,j),0)*h1*(U(i,j)-U(i-1,j))+min(U(i,j),0)*h1*(U(i+1,j)-U(i,j));
              convU(i,j)=convU(i,j)+max(VenU(i,j),0)*h1*(U(i,j)-U(i,j-1))+min(VenU(i,j),0)*h1*(U(i,j+1)-U(i,j));
        else
              convU(i,j)=0.5*h1*U(i,j)*(U(i+1,j)-U(i-1,j));
              convU(i,j)=convU(i,j)+0.5*h1*VenU(i,j)*(U(i,j+1)-U(i,j-1));
        end
        px(i,j)=h1*(p(i,j-1)-p(i-1,j-1));
        tauU(i,j)=h2*(U(i+1,j)-2.*U(i,j)+U(i-1,j))+h2*(U(i,j+1)-2.*U(i,j)+U(i,j-1));
        UenV(j,i)=0.25*(U(j,i)+U(j-1,i)+U(j,i+1)+U(j-1,i+1));
        if centre==0,
              convV(j,i)=max(UenV(j,i),0)*h1*(V(j,i)-V(j-1,i))+min(UenV(j,i),0)*h1*(V(j+1,i)-V(j,i));
              convV(j,i)=convV(j,i)+max(V(j,i),0)*h1*(V(j,i)-V(j,i-1))+min(V(j,i),0)*h1*(V(j,i+1)-V(j,i));
        else
              convV(j,i)=0.5*h1*UenV(j,i)*(V(j+1,i)-V(j-1,i));
              convV(j,i)=convV(j,i)+0.5*h1*V(j,i)*(V(j,i+1)-V(j,i-1));
        end
        py(j,i)=h1*(p(j-1,i)-p(j-1,i-1));
        tauV(j,i)=h2*(V(j+1,i)-2.*V(j,i)+V(j-1,i))+h2*(V(j,i+1)-2.*V(j,i)+V(j,i-1));
    end
end

Up=U+dt*(epsilon*tauU-convU-px+ForceU);
Vp=V+dt*(epsilon*tauV-convV-py+ForceV);

clear px py UenV VenU tauU tauV convU convV