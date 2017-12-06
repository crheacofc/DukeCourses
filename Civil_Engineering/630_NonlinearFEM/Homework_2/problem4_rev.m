%Problem 4
%matlab program to solve two-dimensional N-R problem
clear 
close all
clc
d = zeros(2,1);
%d1 = d(1);
%d2 = d(2);
%final_load = [56 ; 132];
% N = [d1 + 4*d2^3+5*d2^2 ; d1^3+3*d1^2 + 10*d2];
% K = [1 12*d2^2+10*d2; 3*d1^2+6*d1 10];

%Temporal Loading
tol = 1e-12;
n_load_steps = 10000;
pseudo = linspace(1/n_load_steps,1,n_load_steps);

d1_total = zeros(n_load_steps,1);
d2_total = zeros(n_load_steps,1);
for pt = 1:n_load_steps
    Fext = [56*pseudo(pt);132*pseudo(pt)];
    %NR iteration
    %disp(k)
    NR_iterate_counter = 0;
    %Residual_init = load_increment;
    res_norm = 1000000;
    while (abs(res_norm)>tol)
 
        NR_iterate_counter = NR_iterate_counter +  1;
        d1 = d(1); d2 = d(2);

        Fint = [ d1 + 4*d2^3 + 5*d2^2; d1^3 + 3*d1^2 + 10*d2];
        K = [1 12*d2^2 + 10*d2; 3*d1^2 + 6*d1 10 ];
        Residual = Fext-Fint;
        del_d = inv(K)*Residual;
       
        if (NR_iterate_counter ==1)
           %need info about first attempt to compare
           del_init = del_d;
           Residual_init = Residual;
        end
        d = d+del_d;
       res_norm = norm(del_d'*Residual)/norm(del_init'*Residual_init);

    end
    d1_total(pt) = d(1) ; d2_total(pt) = d(2);
    %disp(k)

end

h = figure;
plot(pseudo,d1_total,pseudo,d2_total)
title('displacement vs time')
xlabel('time')
ylabel('displacement')
legend('disp 1','disp 2')
saveas(h,'problem4_10000.png')