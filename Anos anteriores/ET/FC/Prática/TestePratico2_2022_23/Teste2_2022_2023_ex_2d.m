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

T(1,:)=exp((3/4) *t);
T(Nx,:)=exp((-2+3*t)/4);
T(:,1)=exp((-1/2) *x);


C1=(beta*dt)/(2*dx);
C2=(epsilon*dt)/(dx^(2));

n_mat = Nx-2;
aux_a= (1+C2);
aux_b= -0.5*(C2+C1);
aux_c= -0.5*(C2-C1);
diagonal = aux_a*ones(n_mat,1);
diagonal_inf = aux_b *ones(n_mat-1,1);
diagonal_sup = aux_c *ones(n_mat-1,1);
A = diag(diagonal)+diag(diagonal_inf,-1)+diag(diagonal_sup,+1);

b = zeros(Nx-2,1);
tic

aux_1=0.5*(C1+C2);
aux_2=(1-C2);
aux_3=0.5*(C2-C1);
for n = 1:Nt-1
    b = aux_1*T(1:Nx-2,n)+aux_2*T(2:Nx-1,n)+aux_3*T(3:Nx,n);
    b(1)=b(1)+T(1,n+1);
    b(n_mat)=b(n_mat)+T(Nx,n+1);

     T(2:Nx-1,n+1)=linsolve(A,b);

end
toc
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
legend('solucao numerica','solucao analitica')

Erro=abs(Ta(x,0.1)-T(:,100));
figure()
plot(x,Erro)

