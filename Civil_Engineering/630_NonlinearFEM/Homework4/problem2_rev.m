%ASSUMPTIONS
%a simple nonlinear hardening constitutive relation
%sigma(epsilon) = E * epsilon + alpha * E * abs(epsilon) * epsilon
%left end fixed, right end subjected to load
%small deformation, no geometric nonlinearity
%piecewise linear displacement approximation

clear all;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User input
E=10; %modulus in constitutive relation
H = 4;
sig_y = 5;
alpha=100; %constant multiplier for nonlinear term
xrange=[0,1]; %domain
ndiv=50;%number of Finite Element divisions in the domain
finalLoad=10; %final load
NRtol=1e-8; %convergence threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End
%Edited by CLR to do incremental loading
len = xrange(2)-xrange(1);
nelem=ndiv;
nnode=ndiv+1;
nodalcoor=0:len/ndiv:len;
connectivity = zeros(nelem,2);
for i=1:nelem
    connectivity(i,1)=i;
    connectivity(i,2)=i+1;
end

dispVec=zeros(nnode,1);%%% initial disp=zeros
TangentK=zeros(nnode, nnode);
LoadVec=zeros(nnode-1, 1);%%% external load Vec

%incremental loading scheme
n_pseudo_max = 10; %1 pseudo time steps
n_pseudo = 1;
load_init = finalLoad/n_pseudo_max;
load_incr = finalLoad/n_pseudo_max;
        ResidualVec=zeros(nnode-1, 1); %%% initial residual vector
epsilon_history = zeros(n_pseudo_max,nelem);
sigma_history = zeros(n_pseudo_max,nelem);
    % Newton-Raphson at single iteration
        e_p = zeros(n_pseudo_max,nelem);
        Epsilon_p = 0;
        C = 0;
        sig_trial = 0;
        Pho_trial = 0;
        del_gamma = 0;
        while (n_pseudo<n_pseudo_max)

        LoadVec=zeros(nnode-1, 1);%%% external load Vec
        inc = load_init + n_pseudo*load_incr;
        ResidualVec(nnode-1)=inc;
        LoadVec=LoadVec+ResidualVec;
        NRiter=0;
        NRerror=inc;
            while(NRiter<=400) % maximum NR iteration

                TangentK=zeros(nnode, nnode); % initalize tangent stiffness matrix
                for i=1:nelem
                    nodes=connectivity(i, :);
                    node1=nodes(1); x1=nodalcoor(node1);
                    node2=nodes(2); x2=nodalcoor(node2);
                    Bmat=[1/(x1-x2), 1/(x2-x1)];
                    Epsilon = Bmat*dispVec([node1, node2]); %element strain)
                    %check for plasticity
                    e_p_current = e_p(n_pseudo,i);
                    sig_trial = E*(Epsilon-Epsilon_p);
                    Pho_trial = abs(sig_trial)-(sig_y+H*e_p_current);
                    
                    if (Pho_trial <= 0)
                       %elastic
                       sigma_history(n_pseudo+1,i) = sig_trial;
                       epsilon_history(n_pseudo+1,i) = Epsilon;
                       e_p_current(n_pseudo+1,i) = e_p_current;
                       C = E;
                    else
                        %plastic
                        
                        del_gamma = Pho_trial/(E+H);
                        sigma_history(n_pseudo+1,i) = sig_trial-E*del_gamma*sign(sig_trial);
                        e_p_current(n_pseudo+1,i) = e_p_current+del_gamma;
                        Epsilon_p = Epsilon_p + del_gamma*sign(sig_trial);
                        epsilon_history(n_pseudo+1,i) = Epsilon + Epsilon_p;
                        C = E*H/(E+H);
                    end
                    %K=E+alpha*E*2*abs(Epsilon); %update modulus
                    %K = C+alpha*C*2*abs(Epsilon);
                    K = C;
                    Kmat = K*(Bmat'*Bmat)*(x2-x1); %local stiffness matrix
                    TangentK(node1, node1)=TangentK(node1, node1)+Kmat(1,1);
                    TangentK(node1, node2)=TangentK(node1, node2)+Kmat(1,2);
                    TangentK(node2, node1)=TangentK(node2, node1)+Kmat(2,1);
                    TangentK(node2, node2)=TangentK(node2, node2)+Kmat(2,2);
                end

                if(NRiter>=2)
                     ResidualVec = LoadVec-TangentK(2:nnode, 2:nnode)*dispVec(2:nnode);
                end
                if(NRerror<NRtol)
                    break; %% break if error drops below tolerance
                end
                

                deltaDisp=TangentK(2:nnode, 2:nnode)\ResidualVec; 
                deltaDisp=[0; deltaDisp];
                dispVec=dispVec+deltaDisp;
                NRerror=norm(ResidualVec)/inc;
                NRiter=NRiter+1;
            end
            %end Newton-Raphson
        n_pseudo = n_pseudo + 1;
        end
%end incremental loading 
h = figure(1);
plot(nodalcoor,dispVec, '-*')
title('Displacement distribution in the bar')
saveas(h,'/Users/crhea/Dropbox/Duke_courses/Civil_Engineering/630_NonlinearFEM/Homework4/problem2_10.png')


h2 = figure(2);
plot(epsilon_history(:,2),sigma_history(:,2))
title('Stress Strain Relation')
xlabel('Strain')
ylabel('Stress')
saveas(h2,'/Users/crhea/Dropbox/Duke_courses/Civil_Engineering/630_NonlinearFEM/Homework4/problem2_stress-strain.png')





