function [B] = get_B_plane_stress(gradN)
%This function arranges the entries in gradN into a matrix B for plane
%stress elements

dim = size(gradN);
B(1:3,2 * dim(2))=0;

for i=1:dim(2)
    B(1,2*(i-1)+1)=  gradN(1,i);
    B(2,2*i)=  gradN(2,i);
    B(3,2*(i-1)+1) = gradN(2,i); 
    B(3,2*i) = gradN(1,i);
   
end
