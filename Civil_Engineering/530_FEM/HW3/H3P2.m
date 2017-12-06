% 
% We are taking the function f(x) = exp(-10*x). We are approximating the
% function using shape functions. Our element domain is [0,1]. 
%
% Part a) Construct a plot f(x) and F_n^e(x) for n = 2,3, and 5.

function H3P2(~)
    x = sym('x');
    function SF = node_approx(num_nodes,nodes)
        SF = ones(1,num_nodes); %just creating the Shape Function Matrix
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
    end
    function SFval = apply_vals(SF,f_i,num_nodes,x_val)
        SFval = zeros(num_nodes,length(x_val));
       for j = 1:num_nodes
           for xs = 1:length(x_val)
               SFval(j,xs) = double(subs(SF(j),x,x_val(xs))*f_i(j));
           end
       end
    end
    function expression = sf_exp(SF,f,num_nodes)
        expression = zeros(1,num_nodes);
        expression = sym(expression);
        for i = 1:num_nodes
           expression(1,i) = SF(i)*f(i); 
        end
        
    end
    function U = U_total(SF_val,num_nodes)
        U = zeros(1,length(SF_val(1,:)));
        for node = 1:num_nodes
            U = U + SF_val(node,:);
        end
    end
x_val = 0:0.01:1;
num_nodes = 2; %%how many nodes in the element
nodesa = [0.0,1.0];
SFa = node_approx(num_nodes,nodesa);
fa = [exp(-10*0),exp(-10*1)];
SFa = apply_vals(SFa,fa,num_nodes,x_val);
Ua = U_total(SFa,num_nodes);

num_nodes = 3;
nodesb = [0.0,0.5,1.0];
SFb = node_approx(num_nodes,nodesb);
fb = [exp(-10*0),exp(-10*0.5),exp(-10*1)];
SFb = apply_vals(SFb,fb,num_nodes,x_val);
Ub = U_total(SFb,num_nodes);

num_nodes = 5;
nodesc = [0.0,0.25,0.5,0.75,1.0];
SFc = node_approx(num_nodes,nodesc);
fc = [exp(-10*0),exp(-10*.25),exp(-10*0.5),exp(-10*.75),exp(-10*1)];
SFc = apply_vals(SFc,fc,num_nodes,x_val);
Uc = U_total(SFc,num_nodes);

% Finished normal old lagrangian intepolation and application f fi to the
% shape values
% Now for some plotting

plot(x_val,exp(-10*x_val),x_val,Ua,x_val,Ub,x_val,Uc)
title('exp(-10*x) vs approximations')
legend('Exact','2 node','3 node','5 node')
savefig('exp(-10*x)_and_approx.fig')


function E = L2_norm(f,fn)
    E = trapz(x_val,(f-fn).^2);
    E = E^(1/2);
end

%now for the norm business. First lets try to get all 14 version
    function f =f_var()
       f = exp(-10*x); 
    end

fileid = fopen('L2_norm.txt','w');

for n = 2:15
    x_val = 0:0.01:1;
    nodes = linspace(0.0,1.0,n);
    f = zeros(1,n);
    f_reg = zeros(1,n);
    for fval = 1:n
       f(1,fval) = exp(-10*nodes(fval)); 
    end
    for val = 1:length(x_val)
       f_reg(1,val) = exp(-10*val); 
    end
    SF = node_approx(n,nodes);
    %SF_expression = sf_exp(SF,f_var(),n);
    SF = apply_vals(SF,f,n,x_val);
    fn = U_total(SF,n);
    %h = figure;
    %plot(x_val,fn);
    %saveas(h,sprintf('FIG%d.png',k))
    fprintf(fileid,'For %i nodes we get the  L^2 norm value of %6f\n',n,L2_norm(f_reg,fn));

    %disp(L2_norm(f_reg,fn))
    
end
end

