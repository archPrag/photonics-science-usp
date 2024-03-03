clc, clearvars
%problem parameters
initE=1;
omega=1;
gamma=10^(-3)*omega;
couplingR=3;

%my assumptions:
mass=1;

%time constant
tau=2*pi/couplingR/gamma;

%simulating parameters(in seconds)
dt=0.00001;
plotted=5000000;
final=4*tau;


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

  %energy:
energy=zeros(2,size(t,2));
energy(1,:)=mass*(phase(2,:).^2)./2+mass*omega^2*phase(1,:).^2/2;
energy(2,:)=mass*(phase(4,:).^2)./2+mass*omega^2*phase(3,:).^2/2;
regularFall=energy(1,1)*exp(-gamma*t);
plot(t/tau,regularFall,"DisplayName","uncoupled")

title('energy as a function of time')
xlabel("t*Omega²/2/pi/omega")
ylabel("T+V")
hold on

plot(t/tau,energy(1,:),"DisplayName","particle 1")

plot(t/tau,energy(2,:),"DisplayName","particle 2")

hold off

