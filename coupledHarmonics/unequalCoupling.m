clc, clearvars
%problem parameters
initE=1;
omega1=1;
omega2=1;
gamma1=22*10^(-3);
gamma2=10^(-3);
omegaCoupling=0.1;

%my assumptions:
mass=1;

%time constant
tau=2*pi*omega1/omegaCoupling^2;

%simulating parameters(in seconds)
dt=0.00001;
plotted=50000;
final=0.8*tau;


%important arrays:
t=0:(dt*plotted):final;
phase=zeros(4,size(t,2));

%simulating the phase
phase(1,1)=sqrt(2*initE/omega1^2/mass);
phase(2,1)=0;
phase(3,1)=0;
phase(4,1)=0;
matrix=[0 1 0 0;-omega1^2 -gamma1 omegaCoupling^2 0;0 0 0 1;omegaCoupling^2 0 -omega2^2 -gamma2];
itteration=(dt*matrix+eye(4))^plotted;
for index=2:(size(t,2))
  phase(:,index)=itteration*phase(:,index-1);
end

  %energy:
energy=zeros(2,size(t,2));
energy(1,:)=mass*(phase(2,:).^2)./2+mass*omega1^2*phase(1,:).^2/2;
energy(2,:)=mass*(phase(4,:).^2)./2+mass*omega2^2*phase(3,:).^2/2;
regularFall=energy(1,1)*exp(-(gamma1+gamma2)*t/2);
plot(t/tau,regularFall,"DisplayName","uncoupled")

xlabel("t*Omega²/2/pi/omega")
ylabel("T+V")
hold on

plot(t/tau,energy(1,:),"DisplayName","particle 1")

plot(t/tau,energy(2,:),"DisplayName","particle 2")

hold off

