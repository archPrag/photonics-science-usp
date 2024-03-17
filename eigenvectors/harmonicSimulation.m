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
dt=0.1;
final=0.16*tau;


%important arrays:
t=0:(dt):final;
phase=zeros(4,size(t,2));
eigenPhase=zeros(4,1);

%simulating the phase
matrix=[0 1 0 0;-omega^2 -gamma couplingR*omega*gamma 0;0 0 0 1;couplingR*omega*gamma 0 -omega^2 -gamma];
[eigenVectors,diagonal]=eig(matrix);
eigenVectors*eigenVectors'
eigenPhase(:,1)=eigenVectors\[sqrt(2*initE/omega^2/mass);0;0;0];
for index=1:size(t,2)
  phase(:,index)=real(eigenVectors*expm(diagonal*t(index))*eigenPhase);
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
