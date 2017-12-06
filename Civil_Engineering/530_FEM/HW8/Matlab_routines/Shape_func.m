%
% Carter Rhea
% Homework 6
% Start Date: October 19,2016
% This program is for question two
%
function [N] = Shape_func(x,y)    %Define the shape functions:
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

