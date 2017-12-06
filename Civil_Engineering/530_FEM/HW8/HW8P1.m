%%Homework 8 Problem 1
function HW8P1(~)
    x = sym('x'); y = sym('y'); c = sym('c');
    %d = c*[1;-1/2;-1;-1/2;1;-1/2;-1;-1/2];
    %d = c*[1;-1/2;0;0;-1;-1/2;0;0;1;-1/2;0;0;-1;-1/2;0;0];    
    d = c*[1;-1/2;-1;-1/2;1;-1/2;-1;-1/2;0;0;0;-1/2;0;0;0;-1/2];

    N_1 = -(1/4)*(1-x)*(1-y)*(1+x+y);

    N_2 = -(1/4.)*(1+x)*(1-y)*(1-x+y);

    N_3 = -(1/4.)*(1+x)*(1+y)*(1-x-y);

    N_4 = -(1/4.)*(1-x)*(1+y)*(1+x-y);

    N_5 = (1/2.)*(1-x^2)*(1-y);

    N_6 = (1/2.)*(1+x)*(1-y^2);

    N_7 = (1/2.)*(1-x^2)*(1+y);

    N_8 = (1/2.)*(1-x)*(1-y^2);

%The following will create a vector of shape functions given a point
    N = zeros(2,16);
    N = sym('N');
    N(1,1) = N_1;
    N(1,3) = N_2;
    N(1,5) = N_3;
    N(1,7) = N_4;
    N(1,9) = N_5;
    N(1,11) = N_6;
    N(1,13) = N_7;
    N(1,15) = N_8;
    N(2,2) = N_1;
    N(2,4) = N_2;
    N(2,6) = N_3;
    N(2,8) = N_4;
    N(2,10) = N_5;
    N(2,12) = N_6;
    N(2,14) = N_7;
    N(2,16) = N_8;
    
    
    U = N*d;
    U = simplify(U);
    disp(U)

    E_xx = diff(U(1),x);
    E_yy = diff(U(2),y);
    E_xy = 1./2*(diff(U(2),x)+diff(U(1),y));

    disp(E_xx)
    disp(E_yy)
    disp(E_xy)













end