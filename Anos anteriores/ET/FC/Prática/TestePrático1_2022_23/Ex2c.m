%%  Ex2c

clear,clc,close all

%CONSTANTES

a=0.2;
b=0.2;
c=1.8;

tf=100;
h=0.05;

t=0:h:tf;

N=length(t);
x=zeros(1,N);
y=zeros(1,N);
z=zeros(1,N);

x(1)=1;
y(1)=1;
z(1)=1;

%% MÉTODO DE RUNGE-KUTTA

%Equações das derivadas:  (Só é necessário colocar as equações)
fx = @(t,X,Y,Z) -Y-Z;        %derivada em ordem ao tempo de x
fy = @(t,X,Y,Z) X+a*Y;   %derivada em ordem ao tempo de v
fz = @(t,X,Y,Z) b+(X-c)*Z;

for i=1:N-1
                %VER PARTE DE CIMA ( DE CIMA PARA BAIXO)
    r1x= fx(t(i),x(i),y(i),z(i));    %primeiro valor =0 ou seja,
    r1y= fy(t(i),x(i),y(i),z(i));    %r1x = fx(t(i) +0*h, x(i) +0*h*r0x, ...)    
    r1z= fz(t(i),x(i),y(i),z(i));
                                %segundo valor =1/2 ou seja, h*1/2
    r2x= fx(t(i)+h/2 , x(i)+r1x*h/2 , y(i)+r1y*h/2, z(i)+r1z*h/2);
    r2y= fy(t(i)+h/2 , x(i)+r1x*h/2 , y(i)+r1y*h/2, z(i)+r1z*h/2);
    r2z= fz(t(i)+h/2 , x(i)+r1x*h/2 , y(i)+r1y*h/2, z(i)+r1z*h/2);

    x(i+1)=x(i)+ r2x*h; %x(i) + h* ( 0*r1x + 1*r2x)
    y(i+1)=y(i)+ r2y*h; %v(i) + h* ( 0*r1v + 1*r2v)
    z(i+1)=z(i)+ r2z*h;
end

%%GRAFIVOS

figure(1)
plot(t,x)
title('x em função do tempo')
xlabel('tempo (s)')
ylabel('x')

figure(2)
plot(t,y)
title('y em função do tempo')
xlabel('tempo (s)')
ylabel('Corrente (A)')

figure(3)
plot(t,z)
title('z em função do tempo')
xlabel('tempo (s)')
ylabel('z')

figure(4)
plot(x,y)
title('Y em função de x')
xlabel('x')
ylabel('y')
