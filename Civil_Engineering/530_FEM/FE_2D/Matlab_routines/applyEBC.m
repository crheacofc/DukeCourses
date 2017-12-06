function [K_corrected, F_corrected] = applyEBC(K, F, EssentialBcs)
%This function enforces the essential boundary conditions on the global
%system of equations.

%first we need to see which nodes are not in the boundary!
    nodes = [];
    num_total_dof = 12;
    BC_vals = [0,0,0,0]; %zeros for the four DOFS
    for i = 1:num_total_dof
        nodes = [nodes i];
    Nodes_wout_IC = [];
        for i = 1:length(nodes)
            if any(nodes(i)==EssentialBcs)
                continue
            else
                Nodes_wout_IC = [Nodes_wout_IC i];
            end
        end
    end
    %disp(Nodes_wout_IC)
    %disp(len(K.T))
    %now to fix F and K matrices
    K_corrected = zeros(length(Nodes_wout_IC), length(Nodes_wout_IC));
    F_corrected = zeros(length(Nodes_wout_IC),1);
    %disp(F_corrected)
    for i = 1:length(Nodes_wout_IC)
        current_node_without = Nodes_wout_IC(i);
        for j = 1:length(EssentialBcs)
            %disp(len(BC_nodes))
            current_node_with = EssentialBcs(j);
            %disp(current_node_with,current_node_without)
            %now take away BC_val at node with bc from affected nodes
            if j==1 %if the first one affected
                F_corrected(i) = F(current_node_without) - BC_vals(j)*K(current_node_without,current_node_with);
            else %all others based on what we already have for F_corrected just calculated
                F_corrected(i) = F_corrected(i) - BC_vals(j)*K(current_node_without,current_node_with);
            end
        end 
    end 





    for k = 1:length(Nodes_wout_IC)
        row = Nodes_wout_IC(k);
        for c = 1:length(Nodes_wout_IC)
            col = Nodes_wout_IC(c);
            K_corrected(k,c) = K(row,col);
        end
    end
    
   

end

