function [Fg] = getForceFromGravity(NodalCoord, Connectivity, rho, g, thickness)
%This functions computes element force vectors from a gravity body force
%and assembles the contributions into a global vector.
%We want to apply a force of -pg on every y coordinate and nothing on the
%x.
num_node = 4;
    num_dof = num_node*2;
    A = zeros(length(NodalCoord),num_dof); %create assembly matrix which is num_elements x num_dof
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
    end%end e
    
 Fg = zeros(12,1);
 % now to actually apply the load!
    function f_grav = F_GRAV(trac)
        f_grav = zeros(8,1); %set initial
        %traction = zeros(8,1);
        %for i = 1:8
        %    current_coord = Assembly(element,i);%get nodal value
        %    if (mod(current_coord,2) == 0)
        %        traction(i) = trac(2);
        %    else
        %        traction(i) = 0;
        %    end
        %oend %filled traction
        %time to integrate
        X_ip =  [-1/sqrt(3),1/sqrt(3)];
        W_ip = [1,1];
        for k = 1:2
           s = X_ip(k);
           %disp(s)
           Ne = Shape_func(s,s);
           J_s = detJSideQ4(s,s);
           
           f_grav = f_grav + transpose(Ne)*trac*J_s(k)*W_ip(k);
        end
        
    end
 
 
for e = 1:2
   % Coords_el = getElementCoordinates(e, NodalCoord, Connectivity);
    fe = F_GRAV([0;-rho*g])*thickness;
    %print(K_e)

    %print(Coord_mat_el.T) %check that proper values are being returned (they are)

    for i= 1:8
        dof_1 = A(e,i);  % get degrees of freedom
        Fg(dof_1) = Fg(dof_1) +  fe(i);
    %end i loop
    end
%end elemental loop
end


end

