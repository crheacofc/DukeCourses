%ASSUMPTIONS
%a simple nonlinear hardening constitutive relation
%sigma(epsilon) = E * epsilon + alpha * E * abs(epsilon) * epsilon
%left end fixed, right end subjected to load
%small deformation, no geometric nonlinearity
%piecewise linear displacement approximation

clear all;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User input
E=10; %modulus in constitutive relation
alpha=100; %constant multiplier for nonlinear term
xrange=[0,1]; %domain
ndiv=50;%number of Finite Element divisions in the domain
finalLoad=100; %final load
NRtol=1e-8; %convergence threshold
num_load_steps = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End

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

    ResidualVec=zeros(nnode-1, 1); %%% initial residual vector
    ResidualVec(nnode-1)=finalLoad;
    LoadVec=LoadVec+ResidualVec;
    
    NRiter=0;
    NRerror=finalLoad;
    while(NRiter<=500) % maximum NR iteration
        
        TangentK=zeros(nnode, nnode); % initalize tangent stiffness matrix
        for i=1:nelem
            nodes=connectivity(i, :);
            node1=nodes(1); x1=nodalcoor(node1);
            node2=nodes(2); x2=nodalcoor(node2);
            Bmat=[1/(x1-x2), 1/(x2-x1)];
            Epsilon = Bmat*dispVec([node1, node2]); %element strain
            K=E+alpha*E*2*abs(Epsilon); %update modulus
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
        NRerror=norm(ResidualVec)/finalLoad;
        NRiter=NRiter+1;
    end
    
figure(1)
plot(nodalcoor,dispVec, '-*')
title('Displacement distribution in the bar')