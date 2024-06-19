%%  Ex2b

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

options=optimset('Display','off','Tolx',1e-10,'TolFun',1e-10);

for k=1:N-1
    func = @(xv) fcr(xv,x(k),y(k),z(k),h);
    holder = func;
    xv0 = [x(k),y(k),z(k)];

    aux = fsolve(holder,xv0,options);

    x(k+1) = aux(1);
    y(k+1) = aux(2);
    z(k+1) = aux(3);
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

