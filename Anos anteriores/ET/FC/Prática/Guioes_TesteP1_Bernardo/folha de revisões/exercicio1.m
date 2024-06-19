%exercicio 1
%alínea a)
clc;
clear all;
close all;

%constantes
%h pode variar
h=0.02;
t=0:h:20;

N=length(t);

alfa=-0.1;
m=1;
K=1;

x=zeros(1,N);
x(1)=1;

vx = zeros(1,N);
vx(1)= 1;

const = [h/2, K*h/(2*m), 2*alfa];

%opções da função fsolve
options=optimset('Display','off','Tolx',1e-10,'TolFun',1e-10);

%ciclo do método de crank-nicolson
for k=1:N-1

    func = @(xv) fcr(xv,x(k),vx(k),const);
    holder = func;
    xv0 = [x(k),vx(k)];

    aux = fsolve(holder,xv0,options);

    x(k+1) = aux(1);

    vx(k+1) = aux(2);
end

%figuras
figure(1)
plot(t,x);
figure(2)
plot(t,vx);

figure(3)
plot(x,vx);
%%
%alínea b)
xmax = [];
tmax = [];
imax = 0;

for k = 2:N-1
    if and(x(k+1)-x(k)<=0,x(k)-x(k-1)>=0) % condição de máximo
        imax = imax+1;
        POL = lagr(t(k-1:k+1),x(k-1:k+1));
        tmax(imax) = POL(1);
        xmax(imax) = POL(2);
    end
end

A = mean(xmax);
NI = length(tmax);

PLF = polyfit(1:NI,tmax,1);
figure()
plot(1:NI,tmax)
T = PLF(1);