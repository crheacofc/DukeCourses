%Problem 4
%matlab program to solve two-dimensional N-R problem
clear
close all
d = zeros(2,1);
Ntsteps = 100;
t = linspace(1/Ntsteps,1,Ntsteps);
etol = 1.0e-16;

for titer = 1:Ntsteps
   Fext = [ 56*t(titer); 132*t(titer) ];
   enorm = 1000;
   iter = 0;
   while (abs(enorm) > etol)
      iter = iter+1;
      d1 = d(1);  d2 = d(2);
      Fint = [ d1 + 4*d2^3 + 5*d2^2; d1^3 + 3*d1^2 + 10*d2];
      Resid = Fext-Fint;
      Ktan = [1 12*d2^2 + 10*d2; 3*d1^2 + 6*d1 10 ];
      %Ktan = [10 0; 0 10];
      deltad = inv(Ktan)*Resid;
      if (iter == 1)
          deltad0 = deltad;
          Resid0 = Resid;
      end
      enorm = (deltad'*Resid)/(deltad0'*Resid0);
      disp(enorm)
      d = d + deltad;
   end
   d1t(titer) = d(1);  d2t(titer) = d(2);
end

plot(t,d1t,t,d2t)