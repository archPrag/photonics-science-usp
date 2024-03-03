clc,clearvars
time=0:0.1:7;

f=exp(1i*time);

omega=0:0.1:2;

fourrierT=fourier(f,time,omega);

plot(omega,fourrierT)
