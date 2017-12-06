function [] = plotMesh(NodalCoord, Connectivity, symbol )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dim = size(Connectivity);

for i=1:dim(1)
    for j=1:dim(2)
        node = Connectivity(i,j);
        Xel(j) = NodalCoord(node,1);
        Yel(j) = NodalCoord(node,2);
    end
    Xel(j+1) =  Xel(1);
    Yel(j+1) = Yel(1);
    plot(Xel,Yel, symbol);
    hold on
end
end

