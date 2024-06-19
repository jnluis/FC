%Teste Pratico 2023
t0=0;
h=0.001;
tf=2;
t=t0:h:tf;
N=length(t);

%Dados
L=0.25;
C=1e-3;
alfa=0.2;
Ck=C/alfa;

dIa=zeros(1,N);
dIb=zeros(1,N);
Ia=zeros(1,N);
Ib=zeros(1,N);

%Condições iniciais
Ia0=0.2;
Ib0=0;
dIa0=0;
dIb0=0;

Ia(1)=Ia0;
Ib(1)=Ib0;
dIa(1)=dIa0;
dIb(1)=dIa0;

for k=1:N-1
    % Euler-Cromer
    a=(-1/(L*C))*Ia(k)-(1/(L*Ck))*(Ia(k)-Ib(k));
    dIa(k+1)=dIa(k)+a*h;
    Ia(k+1)=Ia(k)+dIa(k+1)*h;

    b=(-1/(L*C))*Ib(k)+(1/(L*Ck))*(Ia(k)-Ib(k));
    dIb(k+1)=dIb(k)+b*h;
    Ib(k+1)=Ib(k)+dIb(k+1)*h;
end
%Solução analítica
w1=(1/sqrt(L*C));
w2=(1/sqrt((L*((1/C)+(2/Ck)))^(-1)));
IA=0.5*(Ia0+Ib0)*cos(w1*t)+0.5*(Ia0-Ib0)*cos(w2*t);
IB=0.5*(Ia0+Ib0)*cos(w1*t)-0.5*(Ia0-Ib0)*cos(w2*t);

figure()
subplot(2,1,1)
title('Corrente Ia em função do Tempo')
plot(t,Ia)
xlabel('Tempo');ylabel('Corrente Ia')
subplot(2,1,2)
title('Corrente Ib em função do Tempo')
plot(t,Ib)
xlabel('Tempo');ylabel('Corrente Ib')

figure()
subplot(2,1,1)
title('Corrente Ia em função do Tempo')
plot(t,IA)
xlabel('Tempo');ylabel('Corrente Ia')
subplot(2,1,2)
title('Corrente Ib em função do Tempo')
plot(t,IB)
xlabel('Tempo');ylabel('Corrente Ib')
%%
%2)
%b)
t0=0;
h=0.05;
tf=100;
t=t0:h:tf;
N=length(t);

%Dados
x0=1;
y0=1;
z0=1;
a=0.2;
b=0.2;
c=1.8;

x=zeros(1,N);
y=zeros(1,N);
z=zeros(1,N);

x(1)=x0;
y(1)=y0;
z(1)=z0;

for k=1:N-1
    %Euler-Implicito
   % x(k+1)=x(k)-(y(k+1)+z(k+1))*h;
   % y(k+1)=y(k)+(x(k+1)+a*y(k+1))*h;
   % z(k+1)=z(k)+(b+(x(k+1)-c)*z(k+1))*h;

   func = @(xv) fcr(xv,x(k),y(k),z(k),const);
xv0 = [x(k),y(k),z(k)];
aux = fsolve(func,xv0,options);
x(k+1) = aux(1);
y(k+1) = aux(2);
z(k+1) = aux(3);

end



function F = fcr(xv,xold,vold,zold,const)
% const(1), const(2) e const(3) estão definidas no programa
%principal.
% xold é x(k) e vold é vx(k).
% xv(1) é x(k+1) e xv(2) é vx(k+1).
F(1)=xv(1)+const(1)*xv(2)+const(1)*xv(3);
F(2)=const(1)*xv(1)+const(2)*xv(2);
F(3)=xv(3)-const(3)-xv(1)*xv(3)+const(4)*xv(3);
end

%%
clc
clear all
close all

t0=0;
h=0.05;
tf=100;
t=t0:h:tf;
N=length(t);

%Dados
x0=1;
y0=1;
z0=1;
a=0.2;
b=0.2;
c=1.8;

x=zeros(1,N);
y=zeros(1,N);
z=zeros(1,N);

x(1)=x0;
y(1)=y0;
z(1)=z0;

fx=@(Y,Z)-y-z;
fy=@(X,Y)X+a.*Y;
fz=@(X,Z)b+(X-c).*Z;

for k=1:N-1
    r1x=fx(y(k),z(k));
    r1y=fy(x(k),y(k));
    r1z=fz(x(k),z(k));

    r2x=fx((y(k)+r1y*h/2),(z(k)+r1z*h/2));
    r2y=fy((x(k)+r1x*h/2),(y(k)+r1y*h/2));
    r2z=fz((x(k)+r1x*h/2),(z(k)+r1z*h/2));

    x(k+1)= x(k) + r2x.*h;
    y(k+1)= y(k) + r2y.*h;
    z(k+1)= z(k) + r2z.*h;
end

plot(x,y,z)


