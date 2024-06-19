clc;
clear all;
close all; 

tf=0.5;
L=1; % [cm]

dx=0.01;
dt=0.00005;

t=0:dt:tf;
x=0:dx:L;
nx=numel(x);
nt=numel(t);

T=zeros(nx,nt);
T(:,1)=x.^2;
T(1,:)= 0; 
T(nx,:)=exp(-t); 

C1 = dt/(dx^2);
C2 =0.5*(dt/dx);
C3= dt;

f=@(x,t) (x^2+(2*x)-2)*exp(-t);

for n=1:nt-1
    for i=2:nx-1
        T(i,n+1)=(C1+C2)*T(i-1,n) + (1-2*C1-2*C3)*T(i,n) + (C1-C2)*T(i+1,n) + C3 * f(i,n);
    end
end

Ta=@(x,t) (x.^2)*exp(-t);

figure(1);
plot(x,T(:,1:1000:end))
xlabel('Posicao X)')
ylim([0 L])
ylabel('Temperatura ºC')

figure(2)
plot(x,T(:,nt),'r',x,Ta(x,tf),'b')
xlabel('x')
ylabel('T')
legend('solucao analitica','solucao numerica')

figure(3);
mesh(t(1:10:end),x,T(:,1:10:end))
xlim([0 tf])
xlabel('Tempo (s)')
ylim([0 L])
ylabel('Posicao x')
zlabel('Temperatura ºC')