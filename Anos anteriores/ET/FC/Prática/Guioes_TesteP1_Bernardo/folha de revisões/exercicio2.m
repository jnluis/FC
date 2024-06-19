%exercicio2)
% clc;
% clear all;
% close all;

%variáveis
m=1;
K=1;
h=0.01;
vo=1;
xo=1;
tfin=50;
alfa=-0.1;
w=sqrt(K/m);

%vetores
t=0:h:tfin;
N=length(t);
x=zeros(1,N);
x(1)=xo;
v=zeros(1,N);
v(1)=vo;

%funções anónimas
fx=@(V) V;
fv=@(X) (-K/m)*(X+2*alfa*X^3);

%método de Runge-Kutta
for k=1:N-1
    r1x = fx(v(k));
    r1v = fv(x(k));
    %e aqui tbm , usas os rv no x pq o x leva com v
    r2x = fx(v(k)+r1v*h/2);
    r2v = fv(x(k)+r1x*h/2);

    r3x = fx(v(k)+r2v*h*3/4);
    r3v = fv(x(k)+r2x*h*3/4);

    x(k+1)=x(k)+h*(2*r1x/9 + r2x/3 + (4*r3x)/9);
    v(k+1)=v(k)+h*(2*r1v/9 + r2v/3 + (4*r3v)/9);
end

figure()
plot(t,x);
figure()
plot(t,v);

%o erro estava nesta formula smh
E = (1/2)*m.*v.^2 + (1/2)*K.*x.^2.*(1+alfa.*x.^2);

figure()
plot(t,E);