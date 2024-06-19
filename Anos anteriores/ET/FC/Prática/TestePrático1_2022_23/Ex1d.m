%%  Ex1d

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

%SOLUÇAO ANALITICA

w1=1/sqrt(L*C);
w2=1/sqrt(L*(1/C+2/Ck)^-1);

IA = 1/2*(Ia(1)+Ib(1))*cos(w1*t) + 1/2*(Ia(1)-Ib(1))*cos(w2*t);
IB = 1/2*(Ia(1)+Ib(1))*cos(w1*t) - 1/2*(Ia(1)-Ib(1))*cos(w2*t);

%%

figure(1)
plot(t,Ia,'b',t,IA,'r')
title('Corrente Ia em função do tempo')
xlabel('tempo (s)')
ylabel('Corrente (A)')
legend('Método de Euler','Solução Analitica')

figure(2)
plot(t,Ib,'b',t,IB,'r')
title('Corrente Ib em função do tempo')
xlabel('tempo (s)')
ylabel('Corrente (A)')
legend('Método de Euler','Solução Analitica')