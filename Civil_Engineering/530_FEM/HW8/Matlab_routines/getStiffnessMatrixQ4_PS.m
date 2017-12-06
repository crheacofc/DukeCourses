function [ Ke ] = getStiffnessMatrixQ4_PS(C, D, thickness)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

weights = [1 ;1];
xIP = [-sqrt(3)/3; sqrt(3)/3];
%weights=[2];
%xIP = [0];


Ke = zeros(8,8);
for i=1:length(xIP)
    for j=1:length(xIP)
	[gradN, detJ] = Grad_N_Mapped_Q4(xIP(i), xIP(j), C);

	[B] = get_B_plane_stress(gradN); %grab different thank gradN since vector
    %disp(B)
        Ke = Ke + B' * D * B * detJ * weights(i) * weights(j) * thickness; 
    end 
end
   

end

