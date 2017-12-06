clc 
clear all
%%%Assignment: Complete the script so that it is functional.

% Mesh definition: Study what is inside this function to determine the
%element connectivity and global node numbering to be used in the problem.
[NodalCoord, Connectivity, EssentialBcs] = getMeshSimple();
%
close all;

% Inputs
E = 200E9;
nu = 0.3;
rho = 7800;
h = 0.1;
g = 9.8;
a = E/(1-nu^2);
D = a * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

% Assemble mesh: YOU NEED TO IMPLEMENT THIS FUNCTION
K = assembleStiffnessMatrix(NodalCoord, Connectivity, D, h);
[eigenvectors,Diagonal_mat] = eig(K);
eigenvalues = diag(Diagonal_mat);
%disp(eigenvalues)
figure;
hold on
plotMesh(NodalCoord, Connectivity,'-');
%for i = 1:length(eigenvalues)
 %plot(eigenvectors(:,i))
%end
%disp(eigenvectors(:,1))
% Get Force Vector: YOU NEED TO IMPLEMENT THIS FUNCTION
Fg = getForceFromGravity(NodalCoord, Connectivity, rho, g, h);

% add point forces
Fpoint= [0 0 0 0 10000 0 0 -10000 0 -10000 0 0]';
F = Fg+ Fpoint;
%Apply Essential Boundary conditions: YOU NEED TO IMPLEMENT THIS FUNCTION
[K, F] = applyEBC(K, F, EssentialBcs);

%Compute the displacements
d = K \ F;
%throw back EBC
count = 1;
total_disp = zeros(12,1);
for  i = 1:12
    if any(i == EssentialBcs)
       total_disp(i) = 0; 
    else
        total_disp(i) = d(count);
        count = count + 1;
    end
end

disp(total_disp)
disp(d)
scale=1e5;
plotdeformed(total_disp, NodalCoord, Connectivity, scale);

%ansys_def = [0 0 6.57E-6 0 8.78E-6 0 7.5E-6 -1.34E-6 8.7E-6 -1.11E-6 9.22E-6 -9.16E-6];
%plotANSYS(ansys_def,NodalCoord,Connectivity,scale)

title('Undeformed vs Deformed Mesh of 2 elements')
legend('- undeformed','-- deformed')

%now to plot the eigens...
%{
for i=1:12
    figure(i);
    eval = eigenvalues(i);
    disp(eval)
    plotMesh(NodalCoord, Connectivity,'-');
    plotdeformed(eigenvectors(:,i), NodalCoord, Connectivity,1);
    %title(strcat('Eigenvalue of ', num2str(eval)));
    %saveas(i,strcat('eigenvector_',num2str(i),'_fullint.png'))
    saveas(i,strcat('eigenvector_',num2str(i),'_redint.png'))

end
%}