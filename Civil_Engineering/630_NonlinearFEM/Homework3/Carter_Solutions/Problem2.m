%%%%%%%%%
% Homework 3 Problem 2: Solve the following problem using a newton raphson
% solver and incremental loading.
% Add line search to problem 1!
%
%%%%%%%%%
clear all
%Initialization
ToL = 1.e-4;
STOL = 1.e-3;
d = zeros(2,1);
d1 = d(1);d2=d(2);

%Incremental Loading
loading_steps = 100;
d1_total = zeros(loading_steps,1);
d1_total(1) = 0.0;
d2_total = zeros(loading_steps,1);
d2_total(1) = 0.0;
t = linspace(1/loading_steps,1,loading_steps);
pt = 1;%pseudo time
nr_iters_final = zeros(loading_steps,1);
Fext_history1 = zeros(loading_steps,1);
Fext_history2 = zeros(loading_steps,2);
while (pt<loading_steps)
    Fext = [t(pt);t(pt)+2*t(pt)^2+2*t(pt)^3];
    d1 = d(1);d2 = d(2);
    %Newton Raphson Solver
    rnorm = 10000;
    nriter = 0;
    while (rnorm>ToL)
        %disp(nriter)
        Fint = [d2-d1;d2+d2^2-d1^2];
        Constan = [-1 1; -2*d1 1+2*d2];
        Residual = Fext - Fint;
        %disp(Residual)
        if (nriter == 0)
            Constan = [-1 1; -2*d1 1+2*d2];
            deltad0 = Constan\Residual;
            Residual0 = Residual;
        end
        
        deltad = Constan\Residual;
        %s solve
        %s_func = @(s)deltad'*(Fext-Fint'*(d+s*deltad));
        s_i = linesearch(d,deltad,Residual,Fext,ToL);
        d = d + s_i*deltad;
        %d = d + deltad;
        nriter = nriter + 1;
        rnorm = norm(deltad'*Residual)/norm(deltad0'*Residual0);
        d1 = d(1);d2 = d(2);
    end
    %disp(d1);disp(d2);
    Fext_history1(pt) =  Fext(1);
    Fext_history2(pt) =  Fext(2);
    pt = pt + 1;
    d1_total(pt) = d1;
    d2_total(pt) = d2;
    nr_iters_final(pt) = nriter;
end  


h = figure;
plot(t,d1_total,t,d2_total)
title('displacement vs time')
xlabel('time')
ylabel('displacement')
xlim([0,1.])
legend('disp 1','disp 2')
saveas(h,'displacements_vs_time_P2.png')



disp(length(d1_total));
h2 = figure;
len = length(d2_total);
plot(d1_total(1:len-1),Fext_history1(1:len-1))
title('Force Displacement Curve')
xlabel('displacement 1')
ylabel('Force')
saveas(h2,'displacements_vs_Force_d1_P2.png')

h3 = figure;
plot(d2_total(1:len-1),Fext_history2(1:len-1))
title('Force Displacement Curve')
xlabel('displacement 2')
ylabel('Force')
saveas(h3,'displacements_vs_Force_d2_P2.png')

function val = bisection(f,a,b,ToL,MAXIT)
    c=(a+b)/2;
    k=0;
    while (abs(c-a)>=ToL) && (k<=MAXIT)
        if f(a)*f(c)<0
            b=c;
        else
            a=c;
        end
        c=(b+a)/2;
        k = k +1;
    end
    val = c;
end

function s_i = linesearch(d,deltad,Residual,Fext,ToL)
%initialize
s_i = 1;
G_init = deltad'*Residual;
d1_step = d(1) + deltad(1);
d2_step = d(2) + deltad(2);
Fint = [d2_step-d1_step; d2_step+d2_step^2-d1_step^2];
G_si = deltad'*(Fext - Fint);%straight from calculation
    if(G_init*G_si<0)%check if there is even a root
      while(abs(G_si/G_init)>ToL && G_si<0)
        d1_step = d(1) + s_i*deltad(1);
        d2_step = d(2) + s_i*deltad(2);
        Fint_s = [d2_step-d1_step; d2_step+d2_step^2-d1_step^2];
        G_si = deltad'*(Fext - Fint_s);
        s_i = (s_i*G_init-0*G_si)/(G_init-G_si);
      end
    else  
    end
end
