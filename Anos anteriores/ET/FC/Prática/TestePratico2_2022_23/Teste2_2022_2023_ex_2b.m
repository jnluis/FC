clc
clear all
close all

L=1;
beta=1;
epsilon=1;
dt=0.001;
dx=0.1;
tf=1;
x=0:dx:L;
t=0:dt:tf;
Nx=length(x);
Nt=length(t);
T=zeros(Nx,Nt);
T(:,1)=exp((-1/2) *x);
T(1,:)=exp((3/4) *t);
T(Nx,:)=exp((-2+3*t)/4);


C1=(beta*dt)/(2*dx);
C2=(epsilon*dt)/(dx^(2));
for n=1:Nt-1
    for i=2:Nx-1
        T(i,n+1)=(1-2*C2)*T(i,n) + (C2-C1)*T(i+1,n) + (C1+C2)*T(i-1,n);
    end
end

%Solução Analítica

Ta=@(x,t)exp((-2*x + 3*t)/4);

figure()
mesh(t(:,1:100:end),x,T(:,1:100:end))
xlabel('t')
ylabel('x')
zlabel('T')

figure()
subplot(1,2,1)
plot(x,T(:,100),'r',x,Ta(x,0.1),'b')
xlabel('x')
ylabel('T')
xlabel('espessura')
ylabel('Temperatura')
legend('solucao analitica','solucao numerica')

Erro=abs(Ta(x,0.1)-T(:,100));
figure()
plot(x,Erro)
