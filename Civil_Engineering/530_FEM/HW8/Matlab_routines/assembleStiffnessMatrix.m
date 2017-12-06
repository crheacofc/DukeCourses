function K = assembleStiffnessMatrix(NodalCoord, Connectivity, D, thickness)
%assembleStiffnessMatrix: This function assembles the global stiffness matrix using
%element contributions
% we will be using the getStiffnessMatrixQ4_PS function from the previous
% homework
    num_node = 4;
    num_dof = num_node*2; %per element
    num_dof_total = 12;
    A = zeros(length(Connectivity),num_dof); %create assembly matrix which is num_elements x num_dof
    
    for e = 1:2 %each element
         for i =1:num_node %all the columns of a row
             for j = 1:2 %num of dof per node
                % disp(j)
                  if (j==1)
                      A(e,(i-1)*2+j) = 2*Connectivity(e,i)-1 ; %x Dof
                    %disp(j)
                  else
                        A(e,(i-1)*2+j) = 2*Connectivity(e,i); %y Dof
                  end
             end %end j
         end %end i
         %disp(e)
    end%end e
%now to the actual assembly process
K = zeros(num_dof_total,num_dof_total);  % K matrix will be num_dof x num_dof
    %Loop over each element
    for e = 1:2
        Coords_el = getElementCoordinates(e, NodalCoord, Connectivity);
        Ke = getStiffnessMatrixQ4_PS(Coords_el, D, thickness);
        %print(K_e)

        %print(Coord_mat_el.T) %check that proper values are being returned (they are)

        for i= 1:num_dof
            dof_1 = A(e,i);  % get degrees of freedom
            for j = 1:num_dof
                dof_2 = A(e,j);
                K(dof_1, dof_2) = K(dof_1,dof_2) +  Ke(i,j);
            % end j loop
            end
        %end i loop
        end
%end elemental loop
    end
end
    



