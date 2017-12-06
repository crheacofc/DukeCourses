function [ NodalCoord, Connectivity, essentialBcs] = getMeshSimple()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Xcoord = [0 10 20 2 5 17];
Ycoord = [0 2 0 6 8 6];
NodalCoord=[Xcoord' Ycoord'];
Connectivity=[1 2 5 4; 
              2 3 6 5];
essentialBcs=[1 2 4 6]; %These are dofs that will be fixed

end

