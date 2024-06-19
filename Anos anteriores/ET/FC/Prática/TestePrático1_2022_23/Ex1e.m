%%  Ex1e

clear,clc,close all

%CONSTANTES

L=0.25;
C=10e-3;
a=0.2;
Ck=C/a;

h=0.001;
t=0:h:2;

N=length(t);
Ia=zeros(1,N);
Ib=zeros(1,N);
DIa=zeros(1,N);
DIb=zeros(1,N);

Ia(1)=0.2;
DIa(1)=0;
Ib(1)=0;
dIb(1)=0;

%EULER-CROMER

for k=1:N-1
    DIa(k+1) = DIa(k)+(-1/(C*L)*Ia(k)-1/(Ck*L)*(Ia(k)-Ib(k)))*h;
    DIb(k+1) = DIb(k)+(-1/(C*L)*Ib(k)+1/(Ck*L)*(Ia(k)-Ib(k)))*h;
    
    Ia(k+1) = Ia(k)+DIa(k+1)*h;   
    Ib(k+1) = Ib(k)+DIb(k+1)*h;
end

figure(1)
plot(t,Ia+Ib)
title('Soma das Correntes')
xlabel('tempo (s)')
ylabel('Corrente (A)')

%FRQUENCIA ANGULAR
ind=find(islocalmax(Ia+Ib));
tt=(t(ind));
for i=2:length(tt)
    Ppra(i-1)=tt(i)-tt(i-1);    %calculo de todos os períodos possiveis
end

Tpratico=mean(Ppra);
fpratico=1/Tpratico;

wpratico = 2*pi*fpratico

w1=1/sqrt(L*C)
w2=1/sqrt(L*(1/C+2/Ck)^-1)
