//
//  main.cpp
//  FE_solver_2D
//
// Based on my readings and notes on the FEM during late July 2016. This program is going to (as of now) simply create the K - global matrix given some element connectivity matrix and K.

//  We want to use triangular elements so K can be easily computed since we dont have to worry about an integral!

//  For now we are using a rectangular mesh!
//  Created by Carter Rhea on 7/30/16.
//  Copyright © 2016 Carter Rhea. All rights reserved.
//

#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;



/*
 The following functions are based on the FEM mesh using tirangles starting 0,0 at the bottom left corner.
 dx0/dy0 --> N_k
 dx1/dy1 --> N_i
 dx2/dy2 --> N_j
 
 //WOULD LIKE TO USE THESE FUNCTIONS AT SOME POINT AS A REFINEMENT!
 
 double dx0(vector<double> P1,vector<double> P2){
 return -P1[1];
 }
 double dx1(vector<double> P1,vector<double> P2){
 return -P2[1]+P1[1];
 }
 double dx2(vector<double> P1,vector<double> P2){
 return P2[1];
 }
 double (*pdx0)(vector<double>,vector<double>) = &dx0;
 double (*pdx1)(vector<double>,vector<double>) = &dx1;
 double (*pdx2)(vector<double>,vector<double>) = &dx2;
 
 
 double dy0(vector<double> P1,vector<double> P2){
 return P1[0];
 }
 double dy1(vector<double> P1,vector<double> P2){
 return P2[0]-P1[0];
 }
 double dy2(vector<double> P1,vector<double> P2){
 return -P2[0];
 }
 double (*pdy0)(vector<double>,vector<double>) = &dy0;
 double (*pdy1)(vector<double>,vector<double>) = &dy1;
 double (*pdy2)(vector<double>,vector<double>) = &dy2;
 */
//This function is going to take in our three points P0,P1,P2 --> I,J,K and spit out the local values. This will be used in the main function to build the K_global.
//Create L and U from LU factorization
tuple<vector<vector<double>> ,vector<vector<double>>> LU(vector<vector<double>> A_matrix,int p,int n){
    //for a 10 by 10 Matrix A
    vector<vector<double>> L(n,vector<double>(n));
    vector<vector<double>> U(n,vector<double>(n));
    for (int i=0;i<n;i++){
        L[i][i] = 1.0; //set diagonal to 1.0
        U[0][i] = A_matrix[0][i];//set the first row of U to the first row of A
    }
    for (int k=1;k<n;k++){
        for (int j=k;j<n;j++){
            double add=0.0;
            for (int s=k-p;s<k;s++){
                add += L[k][s]*U[s][j];
            }
            U[k][j] = A_matrix[k][j] - add;
        }
        for (int i=k+1;i<n;i++){
            double add = 0.0;
            for (int s=k-p;s<k;s++){
                add += L[i][s]*U[s][k];
            }
            L[i][k] = (A_matrix[i][k]-add)/U[k][k];
        }
    }
    
    return make_tuple(L,U);
}


vector<double> solveUnitLowerTriangularBanded(vector<vector<double>> L, vector<double> b, int p, int n){
    vector<double> c(n);
    for (int i=0;i<n;i++){
        double add = 0;
        for (int j=0;j<i;j++){
            add += L[i][j]*c[j];
        }
        c[i] = b[i] - add;
    }
    return c;
}


vector<double> solveUnitUpperTriangularbanded(vector<vector<double>> U,vector<double> c,int p, int n){
    vector<double> x(n);
    for (int i=n-1;i>-1;i--){
        double add=0.0;
        for (int j=i-(p+1);j<n;j++){
            add += U[i][j]*x[j];
        }
        x[i] = (c[i] - add)/U[i][i];
    }
    return x;
}
//shape functions


vector<vector<double>> K_local(vector<double> P0,vector<double> P1,vector<double> P2){
    vector<vector<double>> K_loc(3,vector<double>(3));
    double D = (P1[0]*P2[1]-P2[0]*P1[1]);
    double A = 1./2.;//since we just have triangles of identicel area and base = 1 and height = 1.
    K_loc[0][0] = A*((-P1[1])*(-P1[1]) + (P1[0])*(P1[0]))/(D*D); //(dx0*dx0+dy0*dy0)/(2*A)
    K_loc[0][1] = A*((-P1[1])*(-P2[1] + P1[1])+(P1[0])*(P2[0]-P1[0]))/(D*D);//(dx0*dx1_dy0*dy1)/(2*A)
    K_loc[0][2] = A*((-P1[1])*(P2[1]) + (P1[0])*(-P2[0]))/(D*D);//(dx0*dx2+dy0*dy2)/(2*A)
    K_loc[1][0] = A*((-P2[1]+P1[1])*(-P1[1]) + (P2[0]-P1[0])*(P1[0]))/(D*D);//(dx1*dx0+dy1*dy0)/(2*A)
    K_loc[1][1] = A*((-P2[1]+P1[1])*(-P2[1]+P1[1]) + (P2[0]-P1[0])*(P2[0]-P1[0]))/(D*D);//(dx1*dx1+dy1*dy1)/(2*A)
    K_loc[1][2] = A*((-P2[1]+P1[1])*(P2[1]) + (P2[0]-P1[0])*(-P2[0]))/(D*D);//(dx1*dx2+dy1*dy2)/(2*A)
    K_loc[2][0] = A*((P2[1])*(-P1[1]) + (-P2[0])*(P1[0]))/(D*D);//(dx2*dx0+dy2*dy0)/(2*A)
    K_loc[2][1] = A*((P2[1])*(-P2[1]+P1[1]) + (-P2[0])*(P2[0]-P1[0]))/(D*D);//(dx2*dx1+dy2*dy1)/(2*A)
    K_loc[2][2] = A*((P2[1])*(P2[1]) + (-P2[0])*(-P2[0]))/(D*D);//(dx2*dx2+dy2*dy2)/(2*A)
    
    return K_loc;
}

vector<double> F_local(vector<double> P0,vector<double> P1,vector<double> P2){
    vector<double> F_loc(3);
    F_loc[0] ;
    
    return F_loc;
}

int main(int argc, const char * argv[]) {
    
    
    
    //define our nodal points
    map<int,vector<double>> Node;
    double x_start = 0.0; double y_start = 0.0;
    double x_end = 4.0; double y_end = 2.0;
    int x_steps = 5; int y_steps = 3;
    double h_x = (x_end-x_start)/(x_steps-1);double h_y=(y_end-y_start)/(y_steps-1);
    double x = 0.0; double y = 0.0;
    vector<double> vec;
    
    int count = 0;
    for (int i=0;i<y_steps;i++){
        y = y_start+i*h_y;
        for (int j=0;j<x_steps;j++){
            vec.clear();
            x = x_start+j*h_x;
            vec.push_back(x);vec.push_back(y);
            Node.insert(pair<int,vector<double>>(count,vec));
            count++;
        }
    }
    
    
    //Create global connectivity matrix (cheating the same way again!)
    //for the triangular
    vector<vector<double>> GCM(16,vector<double>(3));
    GCM[0] = {0,1,6};
    GCM[1] = {1,2,7};
    GCM[2] = {2,3,8};
    GCM[3] = {3,4,9};
    GCM[4] = {0,6,5};
    GCM[5] = {1,7,6};
    GCM[6] = {2,8,7};
    GCM[7] = {3,9,8};
    GCM[8] = {5,6,11};
    GCM[9] = {6,7,12};
    GCM[10] = {7,8,13};
    GCM[11] = {8,9,14};
    GCM[12] = {5,11,10};
    GCM[13] = {6,12,11};
    GCM[14] = {7,13,12};
    GCM[15] = {8,14,13};
    
    
    
    
    
    
    
    //Creation of K global matrix given a K matrix and a Cell connectivity matrix.
    
    vector<vector<double>> K_global(15,vector<double>(15));
    for (int j=0;j<16;j++){
        for (int i=0;i<3;i++){
            for(int k=0;k<3;k++){
                K_global[GCM[j][i]][GCM[j][k]] += K_local(Node[GCM[j][0]],Node[GCM[j][1]],Node[GCM[j][2]])[i][k];
            }
        }
    }
    
    //create vector for the Solutions (from boundary conditions)
    vector<double> BC(15);
    BC = {75.0,50.0,50.0,50.0,150.0,100.0,0.0,0.0,0.0,250.,150.0,200.0,200.0,200.0,200.0};
    //create vector for F
    //vector<double> F(15);
    //F = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    
    
    
    //Now to solve!
    //int n_elements = 15;
    //Before I can go about solving this I need to clean up the [K]{X}={F} equation
    //the general idea is to use the known solutions (u) and the k values and esentially just move them over to the right hand side of the equation and make it a smaller system to solve.
    
    //for now we need a list specifying which nodes have solutions and which dont...
    vector<int>w_IC = {0,1,2,3,4,5,6,10,11,12,13,14,15};
    //little code to get the nodal ids without
    vector<int>wo_IC ={};
    for (int i=0;i<Node.size();i++){
        if(std::find(w_IC.begin(), w_IC.end(), i) != w_IC.end()) {
            /* v contains x */
        }
        else {
            wo_IC.push_back(i);
        }
    }
    
    vector<double> F_new(wo_IC.size()); //no F yet going to have to make more robust for that
    vector<vector<double>> K_new(wo_IC.size(),vector<double>(wo_IC.size()));
    //now actually fixing the equations...
    for (int j=0;j<wo_IC.size();j++){
        int current_node_without = wo_IC[j];
        //int countj=0; //cycle through without
        //int counti=0; //cyclle through within nodal ids
        for(int i=0;i<w_IC.size();i++){
            int current_node_with = w_IC[i];
            F_new[j] -= BC[current_node_with]*K_global[current_node_without][current_node_with];
            
        }
    }
    for (int kr=0;kr<wo_IC.size();kr++){
        int row = wo_IC[kr];
        for(int kc=0;kc<wo_IC.size();kc++){
            int col = wo_IC[kc];
            K_new[kr][kc] = K_global[row][col];
        }
    }
    
    
    
    vector<vector<double>> L,U;
    tuple<vector<vector<double>> ,vector<vector<double>>> LUs = LU(K_new, 1, 3);
    L = get<0>(LUs);
    U = get<1>(LUs);
    vector<double> c = solveUnitLowerTriangularBanded(L, F_new, 1, 3); //just BC not F+BC since there is no F --i.e. no additional function
    vector<double> x_sol = solveUnitUpperTriangularbanded(U, c, 1, 3);
    
    for (int i=0;i<3;i++){
        cout<<x_sol[i]<<endl;
    }
    for(int i=0;i<wo_IC.size();i++){
        int sol = wo_IC[i];
        BC[sol] = x_sol[sol];
    }
    
    //now lets put everything together (i.e. reput this into BC.
    
    
    ofstream mini;
    mini.open("/Users/crhea/Desktop/FEM2D.txt");
    
    
    for (int i=0;i<15;i++){
        for (int j=0;j<15;j++){
            mini<<K_global[i][j]<<" ";
        }
        mini<<endl;
    }
    mini.close();
    
    
    ofstream output;
    output.open("/Users/crhea/Desktop/FEM2D.vtk");
    output<< "# vtk DataFile Version 2.0" <<'\n';
    output<<"Simple output for verification 1 "<<'\n';
    output<< "ASCII"<<'\n';
    output<<"DATASET UNSTRUCTURED_GRID"<<'\n';
    output << '\n';
    output<<"POINTS "<< 15 <<" double"<<'\n';
    for (int count=0;count<15;count++){
        output<<Node[count][0]<<" "<<Node[count][1]<<" "<<0.0<<'\n';
    }
    output<<'\n';
    int number_cell_elements = 16;
    output<<"CELLS "<<number_cell_elements<<" "<<(number_cell_elements)*4<<'\n';
    
    for (int i=0;i<number_cell_elements;i++){
        output<<"3"<<" "<<GCM[i][0]<<" "<<GCM[i][1]<<" "<<GCM[i][2]<<'\n';
        
    }
    
    output<<'\n';
    output<<"CELL_TYPES "<<number_cell_elements<<'\n';
    for (int i=0;i<number_cell_elements;i++){
        output<<"5"<<'\n';
    }
    output<<'\n';
    
    output<<"POINT_DATA "<< 15<<'\n';
    output<<"SCALARS Value double"<<'\n';
    output<<"LOOKUP_TABLE default"<<'\n';
    for (int count=0;count<15;count++){
        output<<BC[count]<<'\n';
    }
    output<<'\n';
    
    output.close();
    
    return 0;
}
