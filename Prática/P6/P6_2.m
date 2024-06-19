%% Problema 6.2 - Resolução da equação paraxial
clc;
clear all;
close all;

h=0.05;
N=1024; % 2^10
x=-(N-1)/2*h:h:(N-1)/2*h;

zi=0;
zf=4;
hz=0.02;
z=zi:hz:zf;

Nz=length(z);

q=zeros(Nz,N);
q(1,:)=exp(-x.^2./2); % Forma do Feixe em z=0
%q(1,:)=c1./cosh(x); % Forma do Feixe em z=0

% Espcaço dos k
hk=2*pi/(N*h);
kmax=(N/2-1)*hk;
kmin=(-N/2)*hk;
k=kmin:hk:kmax;

qt0=fftshift(fft(q(1,:))); % Transformada de Fourier
qt=zeros(Nz,N);

for n=1:Nz
    qt=qt0.*exp(-i.*k.^2.*z(n)/2);
    q(n,:)=ifft(ifftshift(qt));
% ifft -> Transformada inversa de Fourier
% para converter de volta para o espaço real
end

figure(1)
intensidade = abs(q).^2;
mesh(x,z, intensidade)
axis([-Inf, Inf, -Inf, Inf, -Inf, Inf])
set(gca,'YDir','reverse');
xlabel('x')
ylabel('z')