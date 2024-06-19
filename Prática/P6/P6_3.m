%% Problema 6.3 - Resolução da equação não linear de Schrödinger
clc;
clear all;
close all;

h=0.05;
N=256;
x=-(N-1)/2*h:h:(N-1)/2*h;

zi=0;
zf=4;
hz=0.02;
z=zi:hz:zf;

Nz=length(z);

% Espcaço dos k
hk=2*pi/(N*h);
kmax=(N/2-1)*hk;
kmin=(-N/2)*hk;
k=kmin:hk:kmax;

q=zeros(Nz,N);
q(1,:)=sech(x); % Forma do Feixe em z=0

qt=zeros(Nz,N);

abstol=ones(1,N);
abstol=1e-9.*abstol;
options=odeset('RelTol',1e-9,'AbsTol',abstol);

qt0=fftshift(fft(q(1,:)));
qtexp0=qt0.*exp(1i.*k.^2*zi/2);
[z,qtexp]= ode45(@f6_3,zi:0.02:zf,qtexp0,options,N,k);
Nz=numel(z);

for n=1:Nz
    qt(n,:)=qtexp(n,:).*exp(-1i.*k.^2.*z(n)/2);
    q(n,:)=ifft(ifftshift(qt(n,:)));
end

figure(1)
perfil=abs(q).^2;
mesh(x,z,perfil)
axis([-Inf, Inf, -Inf, Inf, -Inf, Inf])
set(gca,'YDir','reverse');
xlabel('x')
ylabel('z')
zlabel('perfil')