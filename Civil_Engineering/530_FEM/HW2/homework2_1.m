%%%
% Homework2 for FEM (CEE 530) problem 1. We have to use the given psi
% functions as a basis for our test space. Then we are solving the FEM
% problem which is crazy simple since we have only 2 functions here. The
% formulas for K and F are derived from the weak form of the equation
% (Eu')'+b(x)=0. E is youngs modulus and b has a functional form given
% later.
%%%
function homework2_1()
%set up with symbols
x = sym('x');
E = sym('E');
L = sym('L');
Psi = [sin((pi*x)/(4*L));sin((pi*x)/(2*L))];
%snag derivative and then find K,F, and a
dPsi = diff(Psi,x);
K = int(dPsi*E*transpose(dPsi),x,0,1);
F = int((10*sin(pi*x)*Psi),x,0,1);
a = K\F;
%now I must clean up a to just get my values
a_sub1 = subs(a,L,1);
a_sub2 = subs(a_sub1,E,1);
a_fin = a_sub2;
a = double(a_fin);
disp(a)
%Ok now to apply a to our test functions
Psi_corrected = dot(a,Psi);
Psi_corrected = subs(Psi_corrected,L,1);
x_val = linspace(0,1,100);
%making the curves to plot
U_h = double(subs(Psi_corrected,x,x_val)); %h is approx
U_r = ((10)/(pi^2))*(sin(pi*(x_val))+pi*(x_val));
plot(x_val,U_r,x_val,U_h);

title('Approximate vs exact solutions');
legend('exact','approximate')
savefig('exact_vs_approximated.fig')

sigma_h = diff(Psi_corrected);
sigma_h_val = double(subs(sigma_h,x,x_val));
sigma_r = (10/pi)*(cos(pi*x_val)+1);

plot(x_val,sigma_r,x_val,sigma_h_val);
title('Approximate vs exact sigma');
legend('exact','approximate')
savefig('exact_vs_approximated_sigma.fig')


% ADDED SEPTEMBER 25 for problem 4!
function E = L2_norm(f,fn)
    x_val = linspace(0,1,100);
    E = trapz(x_val,(f-fn).^2);
    E = E^(1/2);
end

fileid = fopen('L2_norm_p2.txt','w');
fprintf(fileid,'The  L^2 norm value is %6f',L2_norm(U_r,U_h));
end
