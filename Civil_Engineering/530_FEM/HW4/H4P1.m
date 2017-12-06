%%Script for homework 4

function H4P1(~)
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

    function U = U_total(SF_val,num_nodes)
        U = zeros(1,length(SF_val(1,:)));
        for node = 1:num_nodes
            U = U + SF_val(node,:);
        end
    end

x_val = 0:0.01:1;
num_nodes = 3; %%how many nodes in the element
nodes = [0.0,0.5,1.0];
SF = node_approx(num_nodes,nodes);
SF1 = double(subs(SF(1),x,x_val));
SF2 = double(subs(SF(2),x,x_val));
SF3 = double(subs(SF(3),x,x_val));
%plot(x_val,SF1,x_val,SF2,x_val,SF3);
SF_p = diff(SF,x);
SF1_p = double(subs(SF_p(1),x,x_val));
SF2_p = double(subs(SF_p(2),x,x_val));
SF3_p = double(subs(SF_p(3),x,x_val));


%%K = int(B^t*E*B,x,0,1)
N = [SF(1),SF(2),SF(3)];
B = [SF_p(1),SF_p(2),SF_p(3)];
N = sym(N);
B = sym(B);
E = 1;
K = int(transpose(B)*E*B,x,0,1.0);

%load vector where b(x) = 10*sin(pi*x/l)
%b = [10.0*sin((pi*0)/1),10*sin((pi*0.5)/1),10.0*sin((pi*1)/1)];
F = int(10.0*sin(pi*x)*transpose(N),x,0,1.0); 
disp(F);
%ESB u(0)=0
% need to clean up k and f
K_new = zeros(2,2);
indices = {2,3}; %%indices without EBC
for r =1:2
    indr = indices{r};
    for c =1:2
        indc = indices{c}; 
        K_new(r,c) = K(indr,indc);
    end
end

f_new = zeros(2,1);
for fval =1:2
    indf = indices{fval};
    f_new(fval) = F(indf) - 0*K(indf,1); %multiply by zero because u(0) = 0  
end
%disp(K)
%disp(K_new)
%disp(b)
%disp(f_new)

U = K_new\f_new;
nodal_vals = [0,U(1),U(2)];
SF_val = apply_vals(SF,nodal_vals,3,x_val);
U_fin = U_total(SF_val,3);
U_r = ((10)/(pi^2))*(sin(pi*(x_val))+pi*(x_val));
plot(x_val,U_r,x_val,U_fin)
title('Approximate vs exact solutions');
legend('exact','approximate')
savefig('exact_vs_approximated.fig')

SF_corrected = dot(nodal_vals,SF);
sigma = diff(SF_corrected);
sigma_fin = double(subs(sigma,x,x_val));
sigma_r = (10/pi)*(cos(pi*x_val)+1);
plot(x_val,sigma_r,x_val,sigma_fin)
title('Approximate vs exact sigma');
legend('exact','approximate')
savefig('exact_vs_approximated_sigma.fig')

%%calcuate L^2 Norm
function E = L2_norm(f,fn)
    E = trapz(x_val,(f-fn).^2);
    E = E^(1/2);
end

fileid = fopen('L2_norm_p4.txt','w');
fprintf(fileid,'The  L^2 norm value is %6f',L2_norm(U_r,U_fin));


% HENCE WE CAN SEE THE L^2 NORM FOR THIS SET OF FUNCTIONS IS CONSIDERABLY BETTER!

% Now lets check out this kernel business...
disp(K*[1;1;1])

end