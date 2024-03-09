clc, clearvars

dDelta=0.03;
deltaReal=-1:dDelta:1;
deltaImaginary=-1:dDelta:1;
omega=zeros(size(deltaReal,2),size(deltaImaginary,2));
g=-.5;

for line=1:size(deltaReal,2)
  for row=1:size(deltaImaginary,2)
    omega(line,row)=(deltaReal(line)+deltaImaginary(row)*1i)*movedSqrt((2*g/(deltaReal(line)+1i*deltaImaginary(row)))^2)/2;
  end

end

surf(deltaReal,deltaImaginary,real(omega)/2/g)
hold on

surf(deltaReal,deltaImaginary,-real(omega)/2/g)

hold off

%I will aproximate omega to omega~ \pm delta/2*(4th(2/2) degree padéAproximand)
function root=movedSqrt(argument)
  root=(16+20*argument+5*argument.^2)/(16+12*argument+argument.^2);
end

