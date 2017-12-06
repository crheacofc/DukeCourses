%%%
% CEE 530 Homework 3 Problem 1:
% 
%4 nodes in element h=10 x0=0
%%%
function H3P1(~)
%Construct shape polynomial using lagrangian interpolation
    %x = 0:0.1:10;
    syms x
    num_nodes = 4; %%how many nodes in the element
    nodes = [0.0,3.33,6.66,10.0];
    SF = ones(1,4); %just creating the Shape Function Matrix
    SF = sym(SF);  %%making it symbolic
    for i = 1:num_nodes
            for k = 1:num_nodes;
                if k==i
                    SF(i) = SF(i) + 0 ;
                else
                    SF(i) = SF(i)*((x-nodes(k))/(nodes(i)-nodes(k))); 
                end
            end
    end
% Finished normal old lagrangian intepolation 
% Now for some plotting
x_val = 0:0.1:10;
SF1 = double(subs(SF(1),x,x_val));
SF2 = double(subs(SF(2),x,x_val));
SF3 = double(subs(SF(3),x,x_val));
SF4 = double(subs(SF(4),x,x_val));
plot(x_val,SF1,x_val,SF2,x_val,SF3,x_val,SF4)
title('Shape Functions');
legend('N1','N2','N3','N4')
savefig('shape_functions.fig')

%Great now on to part B
d = [1 0 1 4];
% we know that u^e = sum(n_i*d_i)
U = 0;
for i = 1:4
   U = U + SF(i)*d(i);
end
%Lets also get u' while we are here
U_p = diff(U);
U_vals = double(subs(U,x,x_val));
U_p_vals = double(subs(U_p,x,x_val));
plot(x_val,U_vals,x_val,U_p_vals)
title('Solution and its derivative');
legend('U','U primed')
savefig('U_Uprimed.fig')

%part c
d = [1 1 1 1];
% we know that u^e = sum(n_i*d_i)
U = 0;
for i = 1:4
   U = U + SF(i)*d(i);
end
%Lets also get u' while we are here
U_p = diff(U);
U_vals = double(subs(U,x,x_val));
U_p_vals = double(subs(U_p,x,x_val));
plot(x_val,U_vals,x_val,U_p_vals)
title('Solution and its derivative with d=(1,1,1,1)');
legend('U','U primed')
savefig('U_Uprimed round 2.fig')
% We can see that for part c this is just a straight bar...
end