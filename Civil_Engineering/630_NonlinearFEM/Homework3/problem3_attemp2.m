%%%%%%%%%
% Homework 3 Problem 3: Solve the following problem using a newton raphson
% solver and incremental loading.
% BFGS method added to problem 2!
%
%%%%%%%%%
clear all
%Initialization
ToL = 1.e-12;
STOL = 0.001;
bfgs_max = 100;
d = zeros(2,1);
d1 = d(1);d2=d(2);

%Incremental Loading
loading_steps = 100;
d1_total = zeros(loading_steps,1);
d2_total = zeros(loading_steps,1);
t = linspace(1/loading_steps,1,loading_steps);
v = zeros(2,bfgs_max); % initialize v
w = zeros(2,bfgs_max); % initialize w

pt = 1;%pseudo time
%del_prev = [0;0];
%Residual_Prev = [0;0];
while (pt<loading_steps)
    fprintf('We are on load step %i \n',pt)
    Fext = [t(pt);t(pt)+2*t(pt)^2+2*t(pt)^3];
    %Newton Raphson Solver
    rnorm = 1;
        %disp(nriter)
        Fint = [d2-d1;d2+d2^2-d1^2];
        %Constan = [-1 1; -2*d1 1+2*d2];
        Residual = Fext - Fint;
        %disp(Residual)
        K = [-1 ,1; -2*d1, 1+2*d2];
        K0 = K;
        deltad = K\Residual;
        %s solve
        %s_func = @(s)deltad'*(Fext-Fint'*(d+s*deltad));
        %s_i = bisection(s_func,0,1,STOL,100);
        s_i = linesearch(d,deltad,Residual,Fext,ToL);
        d = d + s_i*deltad;
        d1 = d(1);d2 = d(2);
        %BFGS
        bfgs_it = 0;
        while(rnorm>ToL && bfgs_it<bfgs_max)
            fprintf('We are currently on BFGS iter %i \n',bfgs_it);
            fprintf('The norm is %f \n',rnorm);
            bfgs_it = bfgs_it+1;
            j = bfgs_it;
            if (bfgs_it == 1)
              deltad0=deltad;
              Residual0=Residual;
            end
            s0 = s_i;
            Residual_Last = Residual;
            Fint = [d2-d1;d2+d2^2-d1^2];
            Residual = Fext-Fint;


            v(:,j) = -deltad/(deltad'*(Residual-Residual_Last));
            alpha = (-s0*(Residual-Residual_Last)'*deltad)/(Residual_Last'*deltad);
            alpha = sqrt(max(alpha,0));
            w(:,j) = (Residual-Residual_Last)+alpha*Residual_Last;
            R_bar = Residual;
            for k=1:j
              R_bar = R_bar + (v(:,j-k+1)'*R_bar)*w(:,j-k+1);
            end
            deltad_bar = K0\R_bar;
            for k=1:j
              deltad_bar = deltad_bar + (w(:,k)'*deltad_bar)*v(:,k);
            end

            deltad=deltad_bar;
            rnorm = (deltad'*Residual_Last)/(deltad0'*Residual0);
            %s_func = @(s)deltad'*(Fext-Fint'*(d+s*deltad));
            %s_i = bisection(s_func,0,1,STOL,100);
            
            s_i = linesearch(d,deltad,Residual,Fext,ToL);
            d = d + s_i*deltad;
            d1 = d(1);d2 = d(2);


        end
  fprintf('We are currently on BFGS iter %i \n',bfgs_it);
  fprintf('The norm is %f \n',rnorm);
fprintf('\n');
%END BFGS
            
  
    %disp(d1);disp(d2);
    pt = pt + 1;
    d1_total(pt) = d1;
    d2_total(pt) = d2;
end  

h = figure;
plot(t,d1_total,t,d2_total)
title('displacement vs time')
xlabel('time')
ylabel('displacement')
xlim([0,1.])
legend('disp 1','disp 2')
saveas(h,'problem3.png')
     

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
