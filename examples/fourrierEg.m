clc,clearvars
dt=0.001;
time=0:dt:7;

f=exp(1i*time);


%fourrier transform in codomain
fourrierT=fftshift(fft(f));

%working on the dommain


omega=2*pi*time/size(time,2)/dt^2;
dOmega=omega(2)-omega(1)
omega=omega-omega(end)/2;


%square modulous(must be one) normalization

modulous=sqrt(abs(fourrierT)*abs(fourrierT')*dOmega)
fourrierT=fourrierT/modulous


%plotting

plot(omega,abs(fourrierT).^2)
