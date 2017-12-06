%%%%%%%%%
% Homework 3 Problem 3: Solve the following problem using a newton raphson
% solver and incremental loading.
% BFGS method added to problem 2!
%
%%%%%%%%%
clear all
%Initialization
ToL = 1.e-12;
STOL = 0.5;
bfgs_max = 30;
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
nr_iters_final = zeros(100,1);
%del_prev = [0;0];
Residual_Prev = [0;0];
s_i_prev = 1;
while (pt<loading_steps)
    Fext = [t(pt);t(pt)+2*t(pt)^2+2*t(pt)^3];
    %Newton Raphson Solver
    rnorm = 100;
    nriter = 0;
    while (rnorm>ToL)
       
        %disp(nriter)
        Fint = [d2-d1;d2+d2^2-d1^2];
        %Constan = [-1 1; -2*d1 1+2*d2];
        Residual = Fext - Fint;
        %disp(Residual)
        if (nriter == 0)
            Constan = [-1 1; -2*d1 1+2*d2];
            deltad0 = inv(Constan)*Residual;
            del_prev = deltad0;
            Residual0 = Residual;
        end
        
        deltad = inv(Constan)*Residual;
        %s solve
        s_func = @(s)deltad'*(Fext-Fint'*(d+s*deltad));
        s_i = bisection(s_func,0,1,STOL,100);
        d = d + s_i*deltad;
        %d = d + deltad;
        nriter = nriter + 1;
        rnorm = norm(deltad'*Residual)/norm(deltad0'*Residual0);
        %check if we can end immediately without continuing to BFGS
        if (rnorm<ToL)
           break
        %if not go into BFGS solve
        else
           %do for either max or when rnorm is small enough
           bfgs_it = 0;
           %del_prev = deltad;
           %Residual_Prev = Residual0;
           %calculate BFGS vectors
           v = del_prev/(del_prev'*(Residual-Residual_Prev));
           denom = (Residual_Prev'*del_prev);
           alpha = sqrt((-s_i_prev*(Residual-Residual_Prev)'*del_prev)*(1/denom));
           w = -(Residual-Residual_Prev)+alpha*Residual_Prev;
           Kinv = inv(Constan);
           l = length(Kinv);
           while(rnorm>ToL || bfgs_it<bfgs_max)
                Kinv = (eye(l)+v*w')*Kinv*(eye(l)+w*v') ;
                deltad = Kinv*Residual;
                d = d+s_i*deltad;
                rnorm = norm(deltad'*Residual)/norm(deltad0'*Residual0);
                bfgs_it = bfgs_it + 1;
               
           end
            
        end
        Residual_Prev = Residual; %just for BFGS vector calculations
        del_prev = deltad;
        s_i_prev = s_i;
        d1 = d(1);d2 = d(2);
    end
    %disp(d1);disp(d2);
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
        
     

function val = bisection(f,a,b,ToL,MAXIT)
    c=(a+b)/2;
    k=0;
    while (c-a>=ToL) && (k<=MAXIT)
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
