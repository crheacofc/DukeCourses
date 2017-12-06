%
% Carter Rhea
% Homework 6
% Start Date: October 19,2016
% This program is for question two
%
function H6P2(~)
    %Define the shape functions:
    function N_1 = N1(x,y)
       N_1 = (1/4.)*(x^2-x)*(y^2-y);
    end
    function N_2 = N2(x,y)
       N_2 = (1/4.)*(x^2+x)*(y^2-y);
    end
    function N_3 = N3(x,y)
       N_3 = (1/4.)*(x^2+x)*(y^2+y);
    end
    function N_4 = N4(x,y)
       N_4 = (1/4.)*(x^2-x)*(y^2+y);
    end
    function N_5 = N5(x,y)
        N_5 = (1/2.)*(1-x^2)*(1+y);
    end 
    function N_6 = N6(x,y)
        N_6 = (1/2.)*(1-x)*(1-y^2);
    end
    function N_7 = N7(x,y)
        N_7 = (1/2.)*(1-x^2)*(1-y);
    end
    function N_8 = N8(x,y)
        N_8 = (1/2.)*(1+x)*(1-y^2);
    end
    %The following will create a vector of shape functions given a point
    function [N] = getN_Q8(x,y)
        N = zeros(1,8);
        N(1) = N1(x,y);
        N(2) = N2(x,y);
        N(3) = N3(x,y);
        N(4) = N4(x,y);
        N(5) = N5(x,y);
        N(6) = N6(x,y);
        N(7) = N7(x,y);
        N(8) = N8(x,y);
    end

     function [G] = getGradN_Q8(x,y)
        G = zeros(2,8);
        G = sym(G);
        G(1,1) = (1/4.)*(2*x-1)*(y.^2-y);
        G(1,2) = (1/4.)*(2*x+1)*(y.^2-y);
        G(1,3) = (1/4.)*(2*x+1)*(y.^2+y);
        G(1,4) = (1/4.)*(2*x-1)*(y.^2+y);
        G(1,5) = -x*(1+y);
        G(1,6) = -(1/2.)*(1-y^2);
        G(1,7) = -x*(1-y);
        G(1,8) = (1/2.)*(1-y^2);
        G(2,1) = (1/4.)*(x.^2-x)*(2*y-1);
        G(2,2) = (1/4.)*(x.^2+x)*(2*y-1);
        G(2,3) = (1/4.)*(x.^2+x)*(2*y+1);
        G(2,4) = (1/4.)*(x.^2-x)*(2*y+1);
        G(2,5) = (1/2.)*(1-x^2);
        G(2,6) = -y*(1-x);
        G(2,7) = -(1/2.)*(1-x.^2);
        G(2,8) = -y*(1+x);
     end
    function K = Stiff(k)
       x = sym('x');
       y = sym('y');
       G = getGradN_Q8(x,y);
       GT = transpose(G);
       K = zeros(8,8);
       Int = (GT*k)*G;
       %disp(size(Int))
       % now lets integrate each element
       for i=1:8
           for j=1:8
               int = matlabFunction(Int(i,j),'Vars',[x y]);
               K(i,j) = integral2(int,-1,1,-1,1);
               %disp(GT(j,i)*k*G(i,j))
           end %end j
       end %end i
    end
    
    disp(Stiff(100))
end

