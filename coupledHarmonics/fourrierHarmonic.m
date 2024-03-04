clc, clearvars
%problem parameters
initE=1;
omega=1;
gamma=10^(-3)*omega;
couplingR=3;

%my assumptions:
mass=1;


%simulating parameters(in rad/s)
dOmega=0.0001;
finalOmega=50;
timeDivisions=10000;

%fullOmega
omegaVar=0:dOmega:2*finalOmega;

%Hence time:
dt=2*pi/size(omegaVar,2)/dOmega;
finalTime=finalOmega/dOmega*dt;


%important arrays:
t=2*pi/size(omegaVar,2)*omegaVar/dOmega^2;
phase=zeros(4,size(t,2));

%simulating the phase
phase(1,1)=sqrt(2*initE/omega^2/mass);
phase(2,1)=0;
phase(3,1)=0;
phase(4,1)=0;
matrix=[0 1 0 0;-omega^2 -gamma couplingR*omega*gamma 0;0 0 0 1;couplingR*omega*gamma 0 -omega^2 -gamma];
itteration=(dt/timeDivisions*matrix+eye(4))^timeDivisions;
for index=2:(size(t,2))
  phase(:,index)=itteration*phase(:,index-1);
end


%fourrier Transform
fourrier1=fftshift(fft(phase(1,:)));
fourrier2=fftshift(fft(phase(3,:)));
%normalization
modulous=abs(fourrier1*fourrier1')*dOmega;
fourrier1=fourrier1/sqrt(abs(fourrier1*fourrier1')*dOmega);
fourrier2=fourrier2/sqrt(abs(fourrier2*fourrier2')*dOmega);

%power:
omegaVar=omegaVar-omegaVar(1,end)/2;
power1=omega^2*gamma*abs(fourrier1).^2;
power2=omega^2*gamma*abs(fourrier2).^2;
plot(omegaVar,power1+power2)
xlim([0.99 1.01]);
ylabel("power coefficient")
xlabel("angular frequency")




