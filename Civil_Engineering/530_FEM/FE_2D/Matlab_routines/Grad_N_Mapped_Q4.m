function [gradN, detJ] = Grad_N_Mapped_Q4(psi, eta, C)
%
% Computes the gradient of the shape functions and the determinant of 
% the Jacobian. This function returns a 2 x nodes matrix and a scalar.  
% The matrix C  contains the coordinates of the nodes arraged as [x1 y1; x2 y2; ...
% xn yn] 
% psi and eta are the coordinates in the parent element
%

% Calculate the Grad(N) matrix in parent Q4/Users/crhea/Downloads/Solution_hw_7/Problem4_MatlabCode/Grad_N_Mapped_Q4.m element
GN    = 0.25 * [eta-1  1-eta   1+eta   -eta-1;
    psi-1  -psi-1  1+psi    1-psi];

J     = GN*C;       % compute the Jacobian matrix
detJ  = det(J);     % determinant of Jacobian

gradN     = J\GN;   % compute the gradient of the shape functions: need to
                    % take J^{-1}*GradN_{parent} to get elemental SF

end
