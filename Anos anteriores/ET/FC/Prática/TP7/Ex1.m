%%  Ex1

clear,clc,close all

%%  Variaveis

m=1;    %massa
a=1;

%solução analitica
n=1;
En=n^2*pi^2/(8*a^2);    %valores exato

tol=10^-2;      %tolerancia
B=0;            %valor pretendido para a energia

%vetor x (dominio entre -a e a)
h=0.001;    
x=-a:h:a;   
N=numel(x);

y=zeros(1,N);
%condiçoes iniciais e fronteira
y(1)=0;
y(2)=h;  
%condiçoes para o shooting  (1.2337 valor base para n=1->valor exato)
E(1)=1;
E(2)=1.5;

%%  Método de shooting e Numerov

for i=1:100

    g=2*E(i);       %energia altera ao longo do tempo

    %Metodo de Numerov
    for k=2:N-1
        C1 = 1+h^2*g/12;
        C2=1-5*h^2*g/12;

        y(k+1)= C1^-1*(-C1*y(k-1)+2*C2*y(k));   %S=0
    end

    yf(i)=y(end);

    %Aplicação do metodo da secante  (E-> guess ; yf-> result)
    if(i>1)
        m=(yf(i)-yf(i-1))/(E(i)-E(i-1));
        E(i+1)=E(i)+(B-yf(i))/m;      %guess(3) = guess(2) + (B-result(2))/m 
    end

    %condição de paragem
    if(abs(E(i+1)-E(i))<tol)
        break;
    end

end

%%  Renormalização

c=trapz(y.^2);      %calculo do integral
ynor=y./sqrt(c);

%%  Gráficos

figure(1)
plot(x,ynor)
xlabel('x')
ylabel('Ynormalizado')
grid on