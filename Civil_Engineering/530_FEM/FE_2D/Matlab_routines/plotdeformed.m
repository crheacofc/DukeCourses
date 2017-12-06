function [] = plotdeformed(d, NodalCoord, Connectivity, scale)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

dim=size(NodalCoord);

for i=1:dim(1)
    dof1 = 2*i-1;
    dof2 = 2*i;
    NodalCoord(i,1) = NodalCoord(i,1) + d(dof1)*scale;
    NodalCoord(i,2) = NodalCoord(i,2) + d(dof2)*scale;
end
hold on 

plotMesh(NodalCoord, Connectivity, '--');

end

