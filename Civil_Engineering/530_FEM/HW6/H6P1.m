%
% Carter Rhea
% Homework 6
% Start Date: October 19,2016
% This program is for question one
%
function H6P1(~)
    %Define the shape functions:
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
        N = zeros(2,8);
        N(1,1) = N1(x,y);
        N(1,3) = N2(x,y);
        N(1,4) = N3(x,y);
        N(1,5) = N4(x,y);
        N(2,2) = N1(x,y);
        N(2,4) = N2(x,y);
        N(2,6) = N3(x,y);
        N(2,8) = N4(x,y);
    end
    
    %disp(getN_Q4(0,1));
    
    %now to get the gradient...
    
    function [G] = getGradN_Q4(x,y)
        G = zeros(2,4);
        G(1,1) = (-1/4.)*(a-y);%N1x
        G(1,2) = (1/4.)*(a-y);%N2x
        G(1,3) = (1/4.)*(a+y);%N3x
        G(1,4) = (-1/4.)*(a+y);%N4x
        G(2,1) = (-1/4.)*(b-x);%N1y
        G(2,2) = (-1/4.)*(b+x);%N2y
        G(2,3) = (1/4.)*(b+x);%N3y
        G(2,4) = (1/4.)*(b-x);%N4y

    end
    %disp(getGradN_Q4(0,1));
    
   
    %Lets now take a look at getting an interpolated value given the nodal
    %values(d)
    function [u,v] = getDisp(d,x,y)
        Vals = getN_Q4(x,y)*d; %multiple N by the displacement vector
        u = Vals(1); %U is the x-component
        v = Vals(2); %V is the y-component
    end
    %and for the desired calculation for part d...
    d = transpose([0 0 0 0 1 1 1 -1]);
    [u,v] = getDisp(d,0.5,-0.5);
    disp(strcat('The x value at (0.5,-0.5) is ' , num2str(u)))
    disp(strcat('The y value at (0.5,-0.5) is ' , num2str(v)))
    
    
    %And now the strains (part e) This is nearly the exact same thing as
    %the previous problem...
    %Note that we have to change up the function from part c a bit to solve
    %this... not an issue though because we still use the values from c!
    
    function G = getGradN_Q4_ep(Grad)
        G = zeros(3,8);
        G(1,1) = Grad(1,1);
        G(1,3) = Grad(1,2);
        G(1,4) = Grad(1,3);
        G(1,5) = Grad(1,4);
        G(2,2) = Grad(2,1);
        G(2,4) = Grad(2,2);
        G(2,6) = Grad(2,3);
        G(2,8) = Grad(2,4);
        G(3,1) = Grad(1,1);
        G(3,2) = Grad(1,2);
        G(3,3) = Grad(1,3);
        G(3,4) = Grad(1,4);
        G(3,5) = Grad(2,1);
        G(3,6) = Grad(2,2);
        G(3,7) = Grad(2,3);
        G(3,8) = Grad(2,4);
    end
    function [epsilon] = getStrain(d,x,y)
        GN4 = getGradN_Q4(x,y);
        epsilon = getGradN_Q4_ep(GN4)*d; %multiple B by the displacement vector
        
    end
    
    epsilon = getStrain(d,0.5,-0.5);
    disp(strcat('The xx-strain value at (0.5,-0.5) is ' , num2str(epsilon(1))))
    disp(strcat('The yy-strain value at (0.5,-0.5) is ' , num2str(epsilon(2))))
    disp(strcat('The xy-strain value at (0.5,-0.5) is ' , num2str(epsilon(3))))
end

