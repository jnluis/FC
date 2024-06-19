%%  Ex4.1 pontos  

clear,clc,close all

%% CONSTANTES

u = 10e-3;
L= 1;
T = 10e3;
n=1;
h=0.0001;

x=0:h:L;
N=length(x);

Dy=zeros(1,N);
y=zeros(1,N);

Dy(1)=0.01;
y(1)=0;

%% METODO DE SHOOTING com EULER-CROMER

%solução analitica de w para agora nao a iremos utilizar
w1 = n*pi/L * sqrt(T/u);  %w1 = 3141.6

%dois valores iniciais de w (dado no enunciado) DIFERENTES DO ANALITICO 
w(1)=2000;
w(2)=3500;

%tolerancia defenida
tol=1e-10;

%para o metodo secante
B=0; %é o resultado pretendido para o problema

for i=1:100

    for k=1:N-1     %APLICAÇÃO DO EULER-CROMER
        Dy(k+1)=Dy(k)-w(i)^2*u/T*y(k)*h;
        y(k+1)=y(k) + Dy(k+1)*h;
    end

    yend(i)=y(end); %ver o último valor de y

    if(i>1)     %Para aplicar o metodo das secantes precisa-se de 2 valores no míni
    
        %APLICAÇÃO DO MÉTODO DA SECANTE
        m=(yend(i)-yend(i-1))/(w(i)-w(i-1));    %formula do declive normal
                                                %result=y    e    guess=w

        w(i+1)=w(i)+(B-yend(i))/m;      %guess(3) = guess(2) + (B-result(2))/m                                
    
        if(abs(w(i+1)-w(i))<tol)   %repetir a integração até que w esteja muito próximo da tolerancia
          break;
        end
    end
end

%% GRÁFICOS

plot(x,y)
title('y em função de x')
xlabel('x')
ylabel('y')

%Trata-se de um método fundamental pois começa em 0 e acaba em 0

%% CALCULO DE Y(L)

yL=y(L);   %y(L) = 0