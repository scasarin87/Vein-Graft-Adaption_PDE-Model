% using the same interpolation method on velocity,
% we move SMC with the medium: extremely viscous flow .

dXsmc=zeros(Ncell,2);

for k=1:Ncell,
    for i=1:N,
        for j=1:N-1,
            deltau1=0.; 
            deltau2=0.;
            dx=abs(x1(i)-Xsmc(k,1)); 
            dy=abs(x(j)-Xsmc(k,2));
            if dx<2.*h,
                if dy<2.*h,
                   deltau1=C2*(1.+cos(C4*dx)); 
                   deltau2=C3*(1.+cos(C5*dy));
                   dXsmc(k,1)=dXsmc(k,1)+dt*u(i,j)*deltau1*deltau2*h*h;
                end
            end
            deltau1=0.; 
            deltau2=0.;
            dx=abs(x(j)-Xsmc(k,1)); 
            dy=abs(x1(i)-Xsmc(k,2));
            if dx<2.*h,
                if dy<2.*h,
                   deltau1=C2*(1.+cos(C4*dx)); 
                   deltau2=C3*(1.+cos(C5*dy));
                   dXsmc(k,2)=dXsmc(k,2)+dt*v(j,i)*deltau1*deltau2*h*h;
                end
            end
        end    
    end
end
  
Xsmc(:,1:2)=Xsmc(:,1:2)+dXsmc(:,1:2);

