%%%%%%%%%
% Homework 3 Problem 1: Solve the following problem using a newton raphson
% solver and incremental loading.
%
%%%%%%%%%
clear all
close all
clc
%Initialization
ToL = 1.e-4;
d = zeros(2,1);
d1 = d(1);d2=d(2);

%Incremental Loading
loading_steps = 100;
d1_total = zeros(loading_steps);
%d1_total(1) = 0.0;
d2_total = zeros(loading_steps);
%d2_total(1) = 0.0;
t = linspace(1/loading_steps,1,loading_steps);
pt = 1;%pseudo time
nr_iters_final = zeros(100,1);

while (pt<loading_steps)
    Fext = [t(pt);t(pt)+2*t(pt)^2+2*t(pt)^3];
    %Newton Raphson Solver
    rnorm = 100;
    nriter = 0;
    while (rnorm>ToL)
        %disp(nriter)
        Fint = [d2-d1;d2+d2^2-d1^2];
        Constan = [-1 1; -2*d1 1+2*d2];
        Residual = Fext - Fint;
        %disp(Residual)
        if (nriter == 0)
            deltad0 = inv(Constan)*Residual;
            Residual0 = Residual;
        end
        deltad = inv(Constan)*Residual;
        d = d + deltad;
        nriter = nriter + 1;
        rnorm = norm(deltad'*Residual)/norm(deltad0'*Residual0);
        d1 = d(1);d2 = d(2);
    end
    %disp(d1);disp(d2);
    pt = pt + 1;
    d1_total(pt) = d1;
    d2_total(pt) = d2;
    nr_iters_final(pt) = nriter;
end  

h = figure;
hold;
plot(t,d1_total,'--')
plot(t,d2_total,'--')
title('displacement vs time')
xlabel('time')
ylabel('displacement')
xlim([0,1])
legend('disp 1','disp 2')
saveas(h,'problem1.png')
        
        