%% Problema 4.1 - Método de shooting — determinação da frequência do primeiro modo normal de vibração
% Feito com RK4 em vez de Euler-Cromer
clc;
clear all;
close all;

T= 1000;
mu= 10^-3; % kg/m
L= 1; %m

w1 = (pi/L)*sqrt(T/mu);
w(1)=2000;
w(2)=3500;

tolera=1e-12;
xf=1;
h=0.001;
x=0:h:xf;
N=numel(x);

Dy0= 1e-3;

figure(1)
xlabel('x')
ylabel('y')

for iw=1:150
   y=zeros(1,N);
   Dy=zeros(1,N);
   Dy(1)=Dy0;

   C= (w(iw)^2)*mu/T;
   fy = @(DY) DY;   
   fDy = @(Y) -C*Y;

   for k=1:N-1
       r1Dy=fDy(y(k));
       r1y=fy(Dy(k));

       r2Dy=fDy(y(k)+r1y*h/2);
       r2y=fy(Dy(k)+r1Dy*h/2);

       r3Dy=fDy(y(k)+r2y*h/2);
       r3y=fy(Dy(k)+r2Dy*h/2);

       r4Dy=fDy(y(k)+r3y*h);
       r4y=fy(Dy(k)+r3Dy*h);

       Dy(k+1)= Dy(k)+(r1Dy+2*r2Dy+2*r3Dy+r4Dy)*h/6;
       y(k+1)=y(k)+(r1y+2*r2y+2*r3y+r4y)*h/6;

   end
    
   plot(x,y);
   xlabel('x'); ylabel('y');
   pause(0.5); % Para ir vendo como nos aproximamos da solução

   yf(iw)=y(end);

   if(iw>1)
       dif=(yf(iw)-yf(iw-1))/(w(iw)-w(iw-1));
       w(iw+1)=w(iw)-yf(iw)/dif; % método da secante
       %fprintf('Ciclo %i : %d \n',iw, abs(w(iw+1)-w(iw)))

       if(abs(w(iw+1)-w(iw)) < tolera)
           fprintf('Convergido no ciclo %i: %d < %d \n',iw, abs(w(iw+1)-w(iw)),tolera)
           %fprintf('Omega = %f s^-1\n', w(iw)
           break
       end
   end
end

fprintf('Omega %d s^-1\n',w(end)) 
fprintf('Omega/Omega_i %d/ \n',w(end)/w1)