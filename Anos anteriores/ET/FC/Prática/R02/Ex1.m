%%  Ex1

clear,clc,close all

%%  CONSTANTES

beta=18;
h=0.01;
x=0:h:7;

N=length(x);

T=zeros(N,1);
DT=zeros(N,1);

T(1)=10e-4;

v(1)=1.9;
v(2)=2;

tol=10e-4;
B=0;

%%  METODO DE SHOOTING COM EULER-CROMER

for i=1:1000
    
    DT(1)=10e-4*v(i);

    for k=1:N-1     %APLICAÇÃO DO EULER-CROMER
        DT(k+1)=DT(k)+(v(i)*DT(k)-beta*exp(-1/T(k))*(1+DT(k)-v(i)*T(k)))*h;
        T(k+1)=T(k) + DT(k+1)*h;
    end

    Tend(i)=DT(end); %ver o último valor de y

    if(i>1)     %Para aplicar o metodo das secantes precisa-se de 2 valores no míni
    
        %APLICAÇÃO DO MÉTODO DA SECANTE
        m=(Tend(i)-Tend(i-1))/(v(i)-v(i-1));    %formula do declive normal
                                                %result=y    e    guess=w

        v(i+1)=v(i)+(B-Tend(i))/m;      %guess(3) = guess(2) + (B-result(2))/m                                
    
        if(abs(v(i+1)-v(i))<tol)   %repetir a integração até que w esteja muito próximo da tolerancia
          break;
        end
    end
end

%% GRÁFICOS

plot(x,T)
title('y em função de x')
xlabel('x')
ylabel('y')