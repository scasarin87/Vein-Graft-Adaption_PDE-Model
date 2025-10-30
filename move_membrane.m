% motion of membrane and feed back:
% move_membrane.m

dXa=zeros(Ma,2);

for k=1:Ma,
    for i=1:N,
        for j=1:N-1,
            deltau1=0.; 
            deltau2=0.;
            dx=abs(x1(i)-Xa(k,1)); 
            dy=abs(x(j)-Xa(k,2));
            if dx<2.*h,
                if dy<2.*h,
                   deltau1=C2*(1.+cos(C4*dx)); 
                   deltau2=C3*(1.+cos(C5*dy));
                   dXa(k,1)=dXa(k,1)+dt*u(i,j)*deltau1*deltau2*h*h;
                end
            end
            deltau1=0.; 
            deltau2=0.;
            dx=abs(x(j)-Xa(k,1)); 
            dy=abs(x1(i)-Xa(k,2));
            if dx<2.*h,
                if dy<2.*h,
                   deltau1=C2*(1.+cos(C4*dx)); 
                   deltau2=C3*(1.+cos(C5*dy));
                   dXa(k,2)=dXa(k,2)+dt*v(j,i)*deltau1*deltau2*h*h;
                end
            end
        end    
    end
end
  
Xa=Xa+dXa;

dXb=zeros(Mb,2);

for k=1:Mb,
    for i=1:N,
        for j=1:N-1,
            deltau1=0.; 
            deltau2=0.;
            dx=abs(x1(i)-Xb(k,1)); 
            dy=abs(x(j)-Xb(k,2));
            if dx<2.*h,
                if dy<2.*h,
                   deltau1=C2*(1.+cos(C4*dx)); 
                   deltau2=C3*(1.+cos(C5*dy));
                   dXb(k,1)=dXb(k,1)+dt*u(i,j)*deltau1*deltau2*h*h;
                end
            end
            deltau1=0.; 
            deltau2=0.;
            dx=abs(x(j)-Xb(k,1)); 
            dy=abs(x1(i)-Xb(k,2));
            if dx<2.*h,
                if dy<2.*h,
                   deltau1=C2*(1.+cos(C4*dx)); 
                   deltau2=C3*(1.+cos(C5*dy));
                   dXb(k,2)=dXb(k,2)+dt*v(j,i)*deltau1*deltau2*h*h;
                end
            end
        end    
    end
end
  
Xb=Xb+dXb;

dXab=zeros(Mab,2);

for k=1:Mab,
    for i=1:N,
        for j=1:N-1,
            deltau1=0.; 
            deltau2=0.;
            dx=abs(x1(i)-Xab(k,1)); 
            dy=abs(x(j)-Xab(k,2));
            if dx<2.*h,
                if dy<2.*h,
                   deltau1=C2*(1.+cos(C4*dx)); 
                   deltau2=C3*(1.+cos(C5*dy));
                   dXab(k,1)=dXab(k,1)+dt*u(i,j)*deltau1*deltau2*h*h;
                end
            end
            deltau1=0.; 
            deltau2=0.;
            dx=abs(x(j)-Xab(k,1)); 
            dy=abs(x1(i)-Xab(k,2));
            if dx<2.*h,
                if dy<2.*h,
                   deltau1=C2*(1.+cos(C4*dx)); 
                   deltau2=C3*(1.+cos(C5*dy));
                   dXab(k,2)=dXab(k,2)+dt*v(j,i)*deltau1*deltau2*h*h;
                end
            end
        end    
    end
end
  
Xab=Xab+dXab;