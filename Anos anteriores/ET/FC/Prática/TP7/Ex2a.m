%%  Ex2a)

clear,clc,close all

%%  Variaveis

l=0;

%solução analitica
n=1;
En=-1/2*n^-2;    %valores exato

tol=10*-7;
B=0;

%vetor r (dominio entre -a e a)
h=0.0001;
rmax=50;
r=0:h:rmax;   
N=length(r);

u=zeros(1,N);
%condiçoes iniciais e fronteira
u(1)=0;
u(N)=0;
u(N-1)=h;
comp=1;
%condiçoes para o shooting  (-0.5 valor base para n=1->valor exato)
E(1)=-0.45;
E(2)=-0.55;

%%  Resultado analítico

Ranalitica = 2*exp(-r);

%%  Método de shooting e Numerov

for i=1:1000
    g(1)=1;
    g(2:N)=-2*((l*(l+1))./(2*r(2:N).^2)-1./r(2:N)-E(i));   

    %Metodo de Numerov
    for k=N-1:-1:3  %ciclo ao contrario (Metodo de Numerov Regressivo)
        C1 = 1+h^2*g(k-1)/12;
        C2=1-5*h^2*g(k)/12;
        C3= 1+h^2*g(k+1)/12;

        u(k-1)= C1.^-1.*(2*C2*u(k)-C3*u(k+1));   %S=0
    end

    u(1)=interp1(r(2:5),u(2:5),0,'spline');
    comp=u(1);
    u0(i)=comp;
     %Aplicação do metodo da secante  (E-> guess ; u-> result)
    if(i>1)
        m=(u0(i)-u0(i-1))/(E(i)-E(i-1));
        E(i+1)=E(i)+(B-u0(i))/m;      %guess(3) = guess(2) + (B-result(2))/m 
    end

    %condição de paragem
    if(abs(E(i+1)-E(i))<tol)
        break;
    end
end