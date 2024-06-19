clc
clear all
close all
ww=1:0.1:1.6;
Nw=length(ww);

for iw=1:Nw
   
h=0.01;
tf=100;
t=0:h:tf;
N=length(t);
y=zeros(N,1);
v=zeros(N,1);

for is=1:150

guess(1)=0.3;
guess(2)=0.5;
alfa=guess(is);
tol=1e-4;    
ro=0.2;
w=ww(iw);
y(1)=0;
v(1)=0.2;

for k=1:N-1

    v(k+1) = v(k) + ( -((y(k).^(3)) - y(k))- alfa*v(k)+ro*cos(w*t(k)) )*h;

    y(k+1)= y(k) + v(k+1)*h;
end

   %Codigo para calcular amplitudes minimas
imin=0;
for j=2:N-1
    if and(v(j+1)>v(j),v(j-1)>=v(j))
        imin=imin+1;
        aux=lagr(t(j-1:j+1),-v(j-1:j+1));
        vmin(imin)=-aux(2);
    end
end
    Velocidade_minima=mean(vmin);

    result(is)= Velocidade_minima; % Condição Fronteira
        B=-0.2;

 if(is > 1)
     m=(result(is)-result(is-1))/(guess(is)-guess(is-1));
     guess(is+1)=guess(is)+ (B-result(is))/m;

     if(abs(B-result(is))<tol)
         break
     end
 end
end

 alfa_final=guess(end);

alfa_saida(iw)=alfa_final;





plot(v,y)




end

figure()
plot(alfa_saida,ww)
