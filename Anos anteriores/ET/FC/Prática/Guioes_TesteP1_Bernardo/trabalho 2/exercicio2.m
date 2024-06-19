%% Orbita de mercúrio
clc;clear;close all;
%% a)
m = 1; %3.285e23; %Kg
GM = 4*pi^2;
tf = 3; h = 0.0001;t = 0:h:tf; %ano

x0 = 0.47;y0 = 0;
vy0 = 8.2;vx0 = 0;

N = length(t);
x = zeros(1,N);x(1) = x0;
y = zeros(1,N);y(1) = y0;
vx = zeros(1,N);vx(1) = vx0;
vy = zeros(1,N);vy(1) = vy0;
r = zeros(1,N);


for k = 1:N-1
    r(k) = sqrt(x(k)^2+y(k)^2);
    
    vx(k+1) = vx(k) - (GM*x(k)*h)/(r(k)^3);
    vy(k+1) = vy(k) -(GM*y(k)*h)/(r(k)^3);
    x(k+1)= x(k) + vx(k+1)*h;
    y(k+1)= y(k) + vy(k+1)*h;
end

figure()
plot(t,x,t,y)
%plot(t,vx,t,vy)

%%
%determinar o periodo
imax = 0;
for k = 2:N-1
    if and(y(k+1)-y(k)<=0,y(k)-y(k-1)>=0)
        imax = imax + 1;
        POL = lagr(t(k-1:k+1),y(k-1:k+1));
        ymax(imax) = POL(2);
        tmax(imax) = POL(1);
    end
end

POLY = polyfit(1:imax,tmax,1);
T = POLY(1)
%%
a=zeros(1,round(T/(2*h))-1);

for k = 1:round(T/(2*h))-1
    t1 = mod(atan2(y(k),x(k)),2*pi)
    t2 = mod(atan2(y(k+1),x(k+1)),2*pi)
    a(k) = r(k+1)^2 * (t2-t1)/2;
end
figure()
plot(t(1:round(T/(2*h)-1)),a)