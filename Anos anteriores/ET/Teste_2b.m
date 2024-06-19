%%
clc
clear all
close all
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
const = [h, 1-h*a, b,c];
options = optimset('Display','off','Tolx',1e-10,'TolFun',1e-10);

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

plot(x,y,z)



function F = fcr(xv,xold,yold,zold,const)
% const(1), const(2) e const(3) estão definidas no programa
%principal.
% xold é x(k) e vold é vx(k).
% xv(1) é x(k+1) e xv(2) é vx(k+1).
F(1)=xv(1)+const(1)*xv(2)+const(1)*xv(3)-xold;
F(2)=const(1)*xv(1)+const(2)*xv(2)-yold;
F(3)=xv(3)-const(3)-xv(1)*xv(3)+const(4)*xv(3)-zold;
end