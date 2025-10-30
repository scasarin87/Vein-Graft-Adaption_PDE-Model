% prediction3.m
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

for i=2:N-1,
    for j=2:N,
        VenU(i,j)=0.25*(V(i,j)+V(i+1,j)+V(i,j-1)+V(i+1,j-1));
        px(i,j)=h1*(p(i,j-1)-p(i-1,j-1));
        UenV(j,i)=0.25*(U(j,i)+U(j-1,i)+U(j,i+1)+U(j-1,i+1));
        py(j,i)=h1*(p(j-1,i)-p(j-1,i-1));
    end
end

Up=U;
Vp=V;

% remark:
% approximation of explicit convective term:
% second order central differences

diff_gs1=1000;
diff_gs2=1000;
k_gs=0;
while (0.5*(diff_gs1+diff_gs2)>1.e-14) && k_gs<1000,
    k_gs=k_gs+1;
    old_Up=Up;
    old_Vp=Vp;
    for i=2:N-1,
        Up(i,2)=const_gs*(U(i,2)-dt*px(i,2)+epsilon*dt*h2*(Up(i+1,2)+Up(i,3)+Up(i-1,2)-Up(i,2))+dt*ForceU(i,2));
        Vp(2,i)=const_gs*(V(2,i)-dt*py(2,i)+epsilon*dt*h2*(Vp(3,i)+Vp(2,i+1)-Vp(2,i)+Vp(2,i-1))+dt*ForceV(2,i));
        for j=3:N-1,
            Up(i,j)=const_gs*(U(i,j)-dt*px(i,j)+epsilon*dt*h2*(Up(i+1,j)+Up(i,j+1)+Up(i-1,j)+Up(i,j-1))+dt*ForceU(i,j));
            Vp(j,i)=const_gs*(V(j,i)-dt*py(j,i)+epsilon*dt*h2*(Vp(j+1,i)+Vp(j,i+1)+Vp(j-1,i)+Vp(j,i-1))+dt*ForceV(j,i));
        end
        Up(i,N)=const_gs*(U(i,N)-dt*px(i,N)+epsilon*dt*h2*(Up(i+1,N)-Up(i,N)+Up(i-1,N)+Up(i,N-1))+dt*ForceU(i,N));
        Vp(N,i)=const_gs*(V(N,i)-dt*py(N,i)+epsilon*dt*h2*(-Vp(N,i)+Vp(N,i+1)+Vp(N-1,i)+Vp(N,i-1))+dt*ForceV(N,i));
    end
    Up(2:N-1,1)=-Up(2:N-1,2);
    Up(2:N-1,N+1)=-Up(2:N-1,N);
    Vp(1,2:N-1)=-Vp(2,2:N-1); 
    Vp(N+1,2:N-1)=-Vp(N,2:N-1);
    diff_gs1=l2norm(old_Up-Up)/l2norm(Up);
    diff_gs2=l2norm(old_Vp-Vp)/l2norm(Vp);
end
clear old_Up old_Vp diff_gs1 diff_gs2 px py UenV VenU