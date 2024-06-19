%1)
%d)
close all;clc;clear all
%constantes
k=16;
m=1;
%valores iniciais
vo=0;
xo=1;

%variável independente
h=0.01;
t=0:h:10;
N=length(t);

%vari+aveis dependentes
x=zeros(1,N);
v=zeros(1,N);
x(1)=xo;
v(1)=vo;

%solução analítica
w  = sqrt(k/m);
xx = xo*cos(w.*t);
vx = -xo*w*sin(w.*t);

%solução numérica
%equações das derviadas
fx = @(t,X,V) V;           %derivada de x em função do tempo
fv = @(t,X,V) -k*X/m;      %derivada da v em função do tempo

for k=1:N-1
    %neste casoo podiamos só por um dos argumentos pq a formula da derivada
    %só usa um também, teria de mudar a função anónima 

    r1x = fx(t(k),x(k),v(k));
    r1v = fv(t(k),x(k),v(k));

    r2x = fx(t(k)+h/2, x(k)+r1x*h/2 , v(k)+r1v*h/2);
    r2v = fv(t(k)+h/2, x(k)+r1x*h/2 , v(k)+r1v*h/2);

    x(k+1)=x(k)+ r2x*h;
    v(k+1)=v(k)+ r2v*h;
end

%plots
figure()
plot(t,xx,t,x);
xlabel('tempo(s)');
ylabel('distâcia(m)');
legend('sol.analítica','sol.numérica');
%escolhi um h tão bom que os valores são literalmente são iguais sou bué
%slay
figure()
plot(t,vx,t,v);
xlabel('tempo(s)');
ylabel('velocidade(m/s)');
legend('sol.analítica','sol.numérica');

%%
%e)
%simplificar seria tirar as variáveis que não são necessárias para as
%formulas 

fx = @(V) V;           %derivada de x em função do tempo
fv = @(X) -k*X/m;      %derivada da v em função do tempo

    r1x = fx(v(k));
    r1v = fv(x(k));

    r2x = fx(v(k)+r1v*h/2);
    r2v = fv(x(k)+r1x*h/2);

    x(k+1)=x(k)+ r2x*h;
    v(k+1)=v(k)+ r2v*h;
%%
%f)
close all;clc;clear all
%---------------------------------método de Euler--------------------------
%constantes
k=16;
m=1;
alfa = -k/m;

%valores iniciais
vo=0;
xo=1;

%variável independente
h=0.010;
teuler=0:h:15;
N=length(teuler);

%vari+aveis dependentes
xeuler=zeros(1,N);
veuler=zeros(1,N);
xeuler(1)=xo;
veuler(1)=vo;

%solução analítica
w  = sqrt(k/m);
xxeuler = xo*cos(w.*teuler);
vxeuler = -xo*w*sin(w.*teuler);

%solução numérica
for k=1:N-1
    veuler(k+1)=veuler(k)+alfa*xeuler(k)*h;
    xeuler(k+1)=xeuler(k)+veuler(k)*h;  
end

%-------------------------método de runge-kutta----------------------------
%voltei a por pq os valores que pediam de h eram diferentes
%variável independente
h=0.01;
t=0:h:15;
N=length(t);

%solução analítica
w  = sqrt(k/m);
xx = xo*cos(w.*t);
vx = -xo*w*sin(w.*t);

%vari+aveis dependentes
xrk=zeros(1,N);
vrk=zeros(1,N);
xrk(1)=xo;
vrk(1)=vo;

%solução numérica
%equações das derviadas
fx = @(t,X,V) V;           %derivada de x em função do tempo
fv = @(t,X,V) -k*X/m;      %derivada da v em função do tempo

for k=1:N-1
    %neste casoo podiamos só por um dos argumentos pq a formula da derivada
    %só usa um também, teria de mudar a função anónima 

    r1x = fx(t(k),xrk(k),vrk(k));
    r1v = fv(t(k),xrk(k),vrk(k));

    r2x = fx(t(k)+h/2, xrk(k)+r1x*h/2 , vrk(k)+r1v*h/2);
    r2v = fv(t(k)+h/2, xrk(k)+r1x*h/2 , vrk(k)+r1v*h/2);

    xrk(k+1)=xrk(k)+ r2x*h;
    vrk(k+1)=vrk(k)+ r2v*h;
end

%plots
figure()
subplot(1,2,1);
plot(teuler,xxeuler,teuler,xeuler);
xlabel('tempo(s)');
ylabel('distâcia(m)');
legend('sol.analítica','sol.euler');

subplot(1,2,2);
plot(t,xx,t,xrk);
xlabel('tempo(s)');
ylabel('distâcia(m)');
legend('sol.analítica','sol.rk');

%escolhi um h tão bom que os valores são literalmente são iguais sou bué
%slay

figure()
subplot(1,2,1);
plot(teuler,vxeuler,teuler,veuler);
xlabel('tempo(s)');
ylabel('velocidade(m/s)');
legend('sol.analítica','sol.euler');

subplot(1,2,2)
plot(t,vx,t,vrk);
xlabel('tempo(s)');
ylabel('velocidade(m/s)');
legend('sol.analítica','sol.rk');

figure()
subplot(1,2,1);
plot(veuler,xeuler);
xlabel('posição(m)');
ylabel('velocidade(m/s)');
legend('sol.euler');

subplot(1,2,2)
plot(vrk,xrk);
xlabel('posição(m)');
ylabel('velocidade(m/s)');
legend('sol.rk');

%calculo da energia
Em_euler = 0.5 * m * w.^2 * xeuler.^2 + 0-5 * m * veuler.^2;
Em_rk = 0.5 * m * w^2 * xrk.^2 + 0-5 * m * vrk.^2;
Em_ark = 0.5 * m * w^2 * xx.^2 + 0-5 * m * vx.^2;
Em_aeuler = 0.5 * m * w^2 * xxeuler.^2 + 0-5 * m * vxeuler.^2;
figure()
subplot(1,2,1)
plot(teuler,Em_aeuler,teuler,Em_euler)
xlabel('tempo(s)');
ylabel('energia(joule)');
legend('sol.analítica','sol.euler');

subplot(1,2,2)
plot(t,Em_ark,t,Em_rk)
xlabel('tempo(s)');
ylabel('energia(joule)');
legend('sol.analítica','sol.rk');