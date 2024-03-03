clc, clearvars
%problem parameters
initE=1;
omega=1;
gamma=10^(-3)*omega;
couplingR=0.5;

%my assumptions:
mass=1;

%time constant
tau=2*pi/couplingR/gamma;

%simulating parameters(in seconds)
dt=0.00001;
plotted=10;
final=0.16*tau;


%important arrays:
t=0:(dt*plotted):final;
phase=zeros(4,size(t,2));

%simulating the phase
phase(1,1)=sqrt(2*initE/omega^2/mass);
phase(2,1)=0;
phase(3,1)=0;
phase(4,1)=0;
matrix=[0 1 0 0;-omega^2 -gamma couplingR*omega*gamma 0;0 0 0 1;couplingR*omega*gamma 0 -omega^2 -gamma];
itteration=(dt*matrix+eye(4))^plotted;
for index=2:(size(t,2))
  phase(:,index)=itteration*phase(:,index-1);
end


%frequency power
fourrier1=fftshift(fft(phase(1,:)))*plotted*dt/2/pi;
fourrier2=fftshift(fft(phase(3,:)))*plotted*dt/2/pi;
omegaVar=2*pi*t/size(t,2)/(plotted*dt)^2;
omegaVar=omegaVar-omegaVar(1,end)/2;
power1=omega^2*gamma*abs(fourrier1).^2;
power2=omega^2*gamma*abs(fourrier2).^2;
plot(omegaVar,power1+power2)
xlim([0.99 1.01])




