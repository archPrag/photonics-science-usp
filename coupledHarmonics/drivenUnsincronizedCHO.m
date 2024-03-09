
clc, clearvars
couplingRatio=2;
omega1=1;%my assumption
gamma=10^(-3);
forceAmplitude=[1;1];

%parameters
dOmega=0.0001;
ddelta=0.01;
%Variables
drivingOmega=0.996:dOmega:1.004;
delta=(-2:ddelta:2)*gamma;

%simulation
matrix=zeros(2,2);
pow=zeros(size(drivingOmega,2),size(gamma,2));
x=zeros(2,1);
for line=1:size(drivingOmega,2)
  for row=1:size(delta,2)
    matrix(1,1)=omega1^2-drivingOmega(line)-1i*gamma*drivingOmega(line);
    matrix(1,2)=-couplingRatio*omega1*gamma;
    matrix(2,1)=-couplingRatio*omega1*gamma;
    matrix(2,2)=(delta(row)+omega1)^2-drivingOmega(line)-1i*gamma*drivingOmega(line);
    x=matrix^(-1)*forceAmplitude;
    pow(line,row)=omega1^2*gamma*abs(x(1,1))^2+(omega1+delta(row))^2*gamma*abs(x(2,1))^2;
  end
end
pow=pow/(sum(pow,"all")*ddelta*dOmega);
%plotting
heatmap(delta/gamma,drivingOmega/omega1,pow)
xlabel("delta/gamma")
ylabel("omega/omega0")


