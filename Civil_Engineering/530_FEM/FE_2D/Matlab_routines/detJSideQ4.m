function [detJside] = detJSideQ4(x, y)
%UNTITLED2 Summary of this function goes here
% psi is coordinate position along the side
% x and y are vectors with the coordinates of the nodes on the side

G = 0.25 * [-2  2]; % for Q4 element

dxdpsi = G*x;
dydpsi = G*y;

detJside = sqrt(dxdpsi.^2 + dydpsi.^2);

end

