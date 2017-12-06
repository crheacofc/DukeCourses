%Calculating the determinant of the jacobain of two elements
function HW_7_4(~)
    %x = sym('x');
    %y = sym('y');
    %first element

    a = 1;
    b = 1;
    function N_1 = N1(x,y)
       N_1 = (1/4.)*(b-x)*(a-y);
    end
    function N_2 = N2(x,y)
       N_2 = (1/4.)*(b+x)*(a-y);
    end
    function N_3 = N3(x,y)
       N_3 = (1/4.)*(b+x)*(a+y);
    end
    function N_4 = N4(x,y)
       N_4 = (1/4.)*(b-x)*(a+y);
    end

    %The following will create a vector of shape functions given a point
    function [N] = getN_Q4(x,y)
        N = zeros(1,4);
        N(1,1) = N1(x,y);
        N(1,2) = N2(x,y);
        N(1,3) = N3(x,y);
        N(1,4) = N4(x,y);
    end
    function [N] = getN_Q4_vec(x,y)
        N = zeros(2,8);
        N(1,1) = N1(x,y);
        N(1,3) = N2(x,y);
        N(1,5) = N3(x,y);
        N(1,7) = N4(x,y);
        N(2,2) = N1(x,y);
        N(2,4) = N2(x,y);
        N(2,6) = N3(x,y);
        N(2,8) = N4(x,y);
    end


    function [gradN, detJ] = Grad_N_Mapped_Q4(x, y, C)
        gradN = zeros(2,4);
        gradN(1,1) = (-1/4.)*(1-y);%N1x
        gradN(1,2) = (1/4.)*(1-y);%N2x
        gradN(1,3) = (1/4.)*(1+y);%N3x
        gradN(1,4) = (-1/4.)*(1+y);%N4x
        gradN(2,1) = (-1/4.)*(1-x);%N1y
        gradN(2,2) = (-1/4.)*(1+x);%N2y
        gradN(2,3) = (1/4.)*(1+x);%N3y
        gradN(2,4) = (1/4.)*(1-x);%N4y
        J = gradN*C;
        detJ = det(J);
    end
    C = [[-10,-10];[5,-10];[7,10];[-10,10]];
    %length(transpose(C))
    [gradN1, detJ1] = Grad_N_Mapped_Q4(0.5,0.5,C);
    disp(gradN1)
    disp(detJ1)
    
    %part b plane stress
    function [B] = get_B_plane_stress(gradN)
        B = zeros(3,8);
        for i=1:length(gradN)
            B(1,2*i-1) = gradN(1,i);
            B(2,2*i) = gradN(2,i);
            B(3,2*i) = gradN(1,i);
            B(3,2*i-1) = gradN(2,i);
        end
    end

    get_B_plane_stress(gradN1)
    %part c - stiffness matrix
    function [Ke] = getStiffnessMatrixQ4_PS (C, D, thickness)
       Ke = zeros(2*length(transpose(C)),2*length(transpose(C)));
       X_ip = [-sqrt(3)/3,sqrt(3)/3];
       W_ip = [1,1];
       for i=1:length(X_ip)
           for j=1:length(X_ip)
               xi_1 = X_ip(i);xi_2=X_ip(j);
               [gradN, detJ] = Grad_N_Mapped_Q4(xi_1, xi_2, C);
               B = get_B_plane_stress(gradN);
               
               Ke = Ke + (transpose(B)*D*B*detJ*W_ip(i)*W_ip(j));
           end %end j
       end%end i
       Ke = thickness*Ke;
        
    end
    E = 1;
    v = 0.3;
    D = (E/(1-v^2))*[[1,v,0];[v,1,0];[0,0,(1-v)/2]];

    thickness = 0.1;
    
    Ke = getStiffnessMatrixQ4_PS (C, D, thickness);
    disp(Ke)
    
    trac = [1;1];
    function [fe] = Traction(trac)
        fe = zeros(8,1);
        X_ip = [-sqrt(3)/3,sqrt(3)/3];
        W_ip = [1,1];
        for i=1:length(X_ip)
            s = X_ip(i);
            N = getN_Q4_vec(s,s);
            for j=1:4
               N(1,j) = 0;
               N(2,j) = 0;
            end
            
            [~, detJ] = Grad_N_Mapped_Q4(s, s, C);
            fe = fe + N'*trac*detJ*W_ip(i); 
        end %end i
        
    end
    
    Traction(trac)
end