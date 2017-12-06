%Nick Andriano
%CEE 630
%Homework 4
%Problem 2

%ASSUMPTIONS
%a simple nonlinear hardening constitutive relation
%sigma(epsilon) = E * epsilon + alpha * E * abs(epsilon) * epsilon
%left end fixed, right end subjected to load
%small deformation, no geometric nonlinearity
%piecewise linear displacement approximation

clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User input
E=10; %modulus in constitutive relation
alpha=100; %constant multiplier for nonlinear term
xrange=[0,1]; %domain
ndiv=50;%number of Finite Element divisions in the domain
finalLoad=10; %final load
NRtol=1e-8; %convergence threshold
g = 0;  %Essential Boundary Conditions
h = 10;  %Natural Boundary Tractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End


%Take in number of loading steps as an input
%n_loads = input('Enter number of loading steps: ');
n_loads = 10;

len = xrange(2)-xrange(1);
nelem=ndiv;
nnode=ndiv+1;
nodalcoor=0:len/ndiv:len;
connectivity = zeros(nelem,2);
for i=1:nelem
    connectivity(i,1)=i;
    connectivity(i,2)=i+1;
end

%initialize counter
count=0;

%--Plasticity parameters--
sigma_y = 5;
H = 4;


%History
Eps_old(nelem,1) = 0;
ep_old(nelem,1) = 0;
Stress_2 = zeros(1,n_loads+1);
Strain_2 = zeros(1,n_loads+1);

%First load for incremental loading
initial_load = finalLoad/n_loads;


ResidualVec=zeros(nnode-1, 1); %%% initial residual vector
dispVec=zeros(nnode,1);%%% initial disp=zeros

%LOADING
for t = initial_load:(finalLoad/n_loads):finalLoad
    %update iterate count
    count = count + 1;
    
    %Create matrices to store strain and stress
    strain = zeros(nelem,1);
    stress = zeros(nelem,1);
    %Initialize plastic strain
    ep_i = zeros(n_loads,1);
    Fint = zeros(nnode,1);
    
    %Reset values for the incremental loading step for t and h
    LoadVec=zeros(nnode-1, 1);%%% external load Vec
    ResidualVec(nnode-1)=t;
    
    %Update load vector with current load
    LoadVec = LoadVec + ResidualVec;
    
    %Begin NR iteration
    NRiter=0;
    NRerror=finalLoad;
    while(NRiter<=500) % maximum NR iteration
        
        TangentK=zeros(nnode, nnode); % initalize tangent stiffness matrix
        Fint = zeros(nnode,1);
        for i=1:nelem
            nodes=connectivity(i, :);
            node1=nodes(1); x1=nodalcoor(node1);
            node2=nodes(2); x2=nodalcoor(node2);
            Bmat=[1/(x1-x2), 1/(x2-x1)];
            Epsilon(i) = Bmat*dispVec([node1, node2]); %element strain
            
            %Check if elastic/plastic
            if count == 1
                stress_trial = E*(Epsilon(i));
                Phi = abs(stress_trial) - (sigma_y);
            else
                stress_trial = E*(Epsilon(i)-Eps_old(i,1));
                Phi = abs(stress_trial)-(sigma_y + H*ep_old(i,1));
            end
            %
            
            %for elastic
            if (Phi <= 0)
                stress(i,1) = stress_trial;
                strain(i,1) = Eps_old(i);
                ep_i(i,1) = ep_old(i,1);
                C_mod = E;
            %for plastic
            else
                %RETURN MAP
                deltaGamma = Phi/(E+H);
                ep_i(i,1) = ep_old(i,1)+deltaGamma;
                %Epsilon = Eps_old(i,1) + deltaGamma*sign(stress_trial);
                stress(i,1) = stress_trial -E*deltaGamma*sign(stress_trial);
                strain(i,1) = Eps_old(i,1)+deltaGamma*sign(stress_trial);
                C_mod = E*H/(E+H);
            end
            
            %ORIGINAL
            %K=E+alpha*E*2*abs(Epsilon); %update modulus
            
            K = C_mod;
            
            Kmat = K*(Bmat'*Bmat)*(x2-x1); %local stiffness matrix
            TangentK(node1, node1)=TangentK(node1, node1)+Kmat(1,1);
            TangentK(node1, node2)=TangentK(node1, node2)+Kmat(1,2);
            TangentK(node2, node1)=TangentK(node2, node1)+Kmat(2,1);
            TangentK(node2, node2)=TangentK(node2, node2)+Kmat(2,2);
            
            Fe_int = Bmat'*stress(i,1)*(x2-x1);
            Fint(node1,1) = Fint(node1,1)+Fe_int(1,1);
            Fint(node2,1) = Fint(node2,1)+Fe_int(2,1);
        end
        Strain_2(count+1) = Epsilon(2);
        
        if(NRiter>=2)
            ResidualVec = LoadVec - Fint(2:nnode,1);
        end
        if(NRerror<NRtol)
            break; %% break if error drops below tolerance
        end
        
        deltaDisp=TangentK(2:nnode, 2:nnode)\ResidualVec;
        deltaDisp=[0; deltaDisp];
        dispVec=dispVec+deltaDisp;
        NRerror=norm(ResidualVec)/t;
        NRiter=NRiter+1;
    end
    
    Stress_2(count+1) = stress(2);
    Eps_old = strain;
    ep_old = ep_i;
    
    
end


%UNLOADING
for t = finalLoad:-(finalLoad/n_loads):0
    %update iterate count
    count = count + 1;
    
    %Create matrices to store strain and stress
    strain = zeros(nelem,1);
    stress = zeros(nelem,1);
    %Initialize plastic strain
    ep_i = zeros(n_loads,1);
    Fint = zeros(nnode,1);
    
    %Reset values for the incremental loading step for t and h
    LoadVec=zeros(nnode-1, 1);%%% external load Vec
    ResidualVec(nnode-1)=t;
    
    %Update load vector with current load
    LoadVec = LoadVec + ResidualVec;
    
    %Begin NR iteration
    NRiter=0;
    NRerror=finalLoad;
    while(NRiter<=500) % maximum NR iteration
        
        TangentK=zeros(nnode, nnode); % initalize tangent stiffness matrix
        Fint = zeros(nnode,1);
        for i=1:nelem
            nodes=connectivity(i, :);
            node1=nodes(1); x1=nodalcoor(node1);
            node2=nodes(2); x2=nodalcoor(node2);
            Bmat=[1/(x1-x2), 1/(x2-x1)];
            Epsilon(i) = Bmat*dispVec([node1, node2]); %element strain
            
            %Check if elastic/plastic
            if count == 1
                stress_trial = E*(Epsilon(i));
                Phi = abs(stress_trial) - (sigma_y);
            else
                stress_trial = E*(Epsilon(i)-Eps_old(i,1));
                Phi = abs(stress_trial)-(sigma_y + H*ep_old(i,1));
            end
            %
            
            %for elastic
            if (Phi <= 0)
                stress(i,1) = stress_trial;
                strain(i,1) = Eps_old(i);
                ep_i(i,1) = ep_old(i,1);
                C_mod = E;
            %for plastic
            else
                %RETURN MAP
                deltaGamma = Phi/(E+H);
                ep_i(i,1) = ep_old(i,1)+deltaGamma;
                %Epsilon = Eps_old(i,1) + deltaGamma*sign(stress_trial);
                stress(i,1) = stress_trial -E*deltaGamma*sign(stress_trial);
                strain(i,1) = Eps_old(i,1)+deltaGamma*sign(stress_trial);
                C_mod = E*H/(E+H);
            end
            
            K = E;
            
            Kmat = K*(Bmat'*Bmat)*(x2-x1); %local stiffness matrix
            TangentK(node1, node1)=TangentK(node1, node1)+Kmat(1,1);
            TangentK(node1, node2)=TangentK(node1, node2)+Kmat(1,2);
            TangentK(node2, node1)=TangentK(node2, node1)+Kmat(2,1);
            TangentK(node2, node2)=TangentK(node2, node2)+Kmat(2,2);
            
            Fe_int = Bmat'*stress(i,1)*(x2-x1);
            Fint(node1,1) = Fint(node1,1)+Fe_int(1,1);
            Fint(node2,1) = Fint(node2,1)+Fe_int(2,1);
        end
        Strain_2(count+1) = Epsilon(2);
        
        if(NRiter>=2)
            ResidualVec = LoadVec - Fint(2:nnode,1);
        end
        if(NRerror<NRtol)
            break; %% break if error drops below tolerance
        end
        
        deltaDisp=TangentK(2:nnode, 2:nnode)\ResidualVec;
        deltaDisp=[0; deltaDisp];
        dispVec=dispVec+deltaDisp;
        NRerror=norm(ResidualVec)/t;
        NRiter=NRiter+1;
    end
    
    Stress_2(count+1) = stress(2);
    Eps_old = strain;
    ep_old = ep_i; 
end


%Plot displacement in bar
p1 = figure(1);
plot(nodalcoor,dispVec, '-*');
title('Displacement distribution in the bar');
saveas(p1,'/Users/crhea/Dropbox/Duke_courses/Civil_Engineering/630_NonlinearFEM/Homework4/problem2_10.png')

%Plot stress-strain for second element
p2 = figure(2);
plot(Strain_2,Stress_2,'-*');
title('Stress-Strain Plot');
xlabel('Strain');
ylabel('Stress');
ylim([0 12]);
saveas(p2,'/Users/crhea/Dropbox/Duke_courses/Civil_Engineering/630_NonlinearFEM/Homework4/problem2_stress-strain.png')







