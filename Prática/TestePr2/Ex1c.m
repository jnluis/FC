clc;
clear all;
close all;

T= 1000;
L= 0.35; %m
vsom=345; % velocidade do som

xf=L;
h=0.0001;
x=0:h:xf;
N=numel(x);

y=zeros(N,1);
Dy=zeros(N,1);

y(1)=5*10^-6;

tolera=1e-10;

omega = (pi/L)*vsom;
w4= ((4*pi)/L) * vsom;

omega(1)= w4 + (w4*0.1);
omega(2)= w4 - (w4*0.1);

for iw=1:150
   y=zeros(1,N);
   Dy=zeros(1,N);
   Dy(1)=0;
   
   aux= omega(iw)^2
    for k=1:N-1
        Dy(k+1)=Dy(k)+ ((aux / vsom) * y(k)) *h;
        y(k+1) = y(k)+Dy(k+1)*h;
    end
    
   plot(x,y);
   xlabel('x'); ylabel('y');

   yf(iw)=Dy(end); % CondiÁ„o fronteira

   if(iw>1)
       dif=(yf(iw)-yf(iw-1))/(omega(iw)-omega(iw-1));
       omega(iw+1)=omega(iw)-yf(iw)/dif; % m√©todo da secante

       if(abs(omega(iw+1)-omega(iw)) < tolera)
           fprintf('Convergido no ciclo %i: %d < %d \n',iw, abs(omega(iw+1)-omega(iw)),tolera)
           break
       end
   end
end

plot(x,Dy);
xlabel('x'); ylabel('y');