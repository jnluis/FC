%%  Ex2d

clear,clc,close all

%CONSTANTES

a=0.2;
b=0.2;
c=1.8;

tf=100;
h=0.05;

t=0:h:tf;

x0=1;
y0=1;
z0=1;

reltol = 3e-14;
abstol_1 = 1e-13;
abstol_2 = 1e-13;
abstol_3 = 1e-13;

options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2 abstol_3]);
[t,solucao]=ode45(@f,0:h:tf,[x0,y0,z0],options);

tt=t';
xx=solucao(:,1)';
yy=solucao(:,2)';
zz=solucao(:,3)';

%%GRAFIVOS

figure(1)
plot(tt,xx)
title('x em função do tempo')
xlabel('tempo (s)')
ylabel('x')

figure(2)
plot(tt,yy)
title('y em função do tempo')
xlabel('tempo (s)')
ylabel('Corrente (A)')

figure(3)
plot(tt,zz)
title('z em função do tempo')
xlabel('tempo (s)')
ylabel('z')

figure(4)
plot(xx,yy)
title('Y em função de x')
xlabel('x')
ylabel('y')
