%% Problema 1.2 Bola de ténis - topSpin e backSpin
close
clear
clc

% Dados:
g=9.8;
ro=1.225; %densidade
m=0.057;  %massa (kg)
d=0.067;  %diâmetro (m)
R=d/2;
alt=0.7; %metros
B=0.5*ro*pi*d^2/4;
v0=20;
ang=5; %graus
t0=0;
tf=2;
h=0.001;

%% a) sem rotação

t=t0:h:tf;
N=numel(t);
r=zeros(N,3);
v_vec=zeros(N,3);
r(1,:)=[0 0 alt];
v_vec(1,:)=[v0*cosd(ang) 0 v0*sind(ang)];
wy=0; % 0 rotações
w=abs(wy);
w_vers=[0 sign(wy) 0];

for i=1:N-1
    v=norm(v_vec(i,:));
    S=R*w/v;
    CD=0.508+(22.503+4.196*S^-2.5)^-0.4;
    FD=-B*CD*v*v_vec(i,:);
    CL=(2.022+0.981*S^-1)^-1;
    FL=B*CL*v*cross(w_vers,v_vec(i,:));
    a=(FD+FL)/m+[0 0 -g];
    v_vec(i+1,:)=v_vec(i,:)+a*h;
    r(i+1,:)=r(i,:)+v_vec(i,:)*h;
    
    if r(i+1,3)<0
       break
    end
end

t=t(1:i+1);
r=r(1:i+1,:);
v_vec=v_vec(1:i+1,:);

alcance_a=interp1(r(:,3),r(:,1),0,'linear')
alt_max_a=max(r(:,3))

%% b) topspin de 3000 rpm

t=t0:h:tf;
N=numel(t);
r=zeros(N,3);
v_vec=zeros(N,3);
r(1,:)=[0 0 alt];
v_vec(1,:)=[v0*cosd(ang) 0 v0*sind(ang)];
wy=3000*2*pi/60; % 3000 rpm em rad/s
w=abs(wy);
w_vers=[0 sign(wy) 0];

for i=1:N-1
    v=norm(v_vec(i,:));
    S=R*w/v;
    CD=0.508+(22.503+4.196*S^-2.5)^-0.4;
    FD=-B*CD*v*v_vec(i,:);
    CL=(2.022+0.981*S^-1)^-1;
    FL=B*CL*v*cross(w_vers,v_vec(i,:));
    a=(FD+FL)/m+[0 0 -g];
    v_vec(i+1,:)=v_vec(i,:)+a*h;
    r(i+1,:)=r(i,:)+v_vec(i,:)*h;
    
    if r(i+1,3)<0
       break
    end
end

t=t(1:i+1);
r=r(1:i+1,:);
v_vec=v_vec(1:i+1,:);

alcance_b=interp1(r(:,3),r(:,1),0,'linear')
alt_max_b=max(r(:,3))

%% c) backspin de 3000 rpm

t=t0:h:tf;
N=numel(t);
r=zeros(N,3);
v_vec=zeros(N,3);
r(1,:)=[0 0 alt];
v_vec(1,:)=[v0*cosd(ang) 0 v0*sind(ang)];
wy=-3000*2*pi/60; % -3000 rpm em rad/s
w=abs(wy);
w_vers=[0 sign(wy) 0];

for i=1:N-1
    v=norm(v_vec(i,:));
    S=R*w/v;
    CD=0.508+(22.503+4.196*S^-2.5)^-0.4;
    FD=-B*CD*v*v_vec(i,:);
    CL=(2.022+0.981*S^-1)^-1;
    FL=B*CL*v*cross(w_vers,v_vec(i,:));
    a=(FD+FL)/m+[0 0 -g];
    v_vec(i+1,:)=v_vec(i,:)+a*h;
    r(i+1,:)=r(i,:)+v_vec(i,:)*h;
    
    if r(i+1,3)<0
       break
    end
end

t=t(1:i+1);
r=r(1:i+1,:);
v_vec=v_vec(1:i+1,:);

alcance_c=interp1(r(:,3),r(:,1),0,'linear')
alt_max_c=max(r(:,3))