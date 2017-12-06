//
//  main.cpp
//  FE_solver_1D
//  This is created during the summer before my first year as a graduate student. I am going to model this solver off the one I created in python for Math 445 at COFC in order to review the FEM for a one dimensional problem and allow me to practice my c++ programming.
//  Created by Carter Rhea on 7/28/16.
//  Copyright © 2016 Carter Rhea. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <vector>
#include <tuple>
#include <fstream>
#include <sstream>

using namespace std;



double s_i(double x){
    return (1./2.)*(1.-x);
}

double s_j(double x){
    return (1./2.)*(1.+x);
}

double ds_i(){
    return -1./2.;
}

double ds_j(){
    return 1./2.;
}


double K1(double J, double x){
    double val=0.0;
    int q=0;
    val = ((ds_i()*(1./J))*(ds_i()*(1./J))  +  q*s_i(x)*s_i(x))*J;
    return val;
}


double K2(double J, double x){
    double val;
    int q=0;
    
    val = ((ds_i()*(1./J))*(ds_j()*(1./J))+q*s_i(x)*s_j(x))*J;
    return val;
}

double  F1(double x_i,double x_j,double J,double x){
    double val = sin(((x_j-x_i)/2.)*x+(x_i+x_j)/2.)*s_i(x)*J;
    return val;
}

double  F2(double x_i,double x_j,double J,double x){
    double val = sin(((x_j-x_i)/2.)*x+(x_i+x_j)/2.)*s_j(x)*J;
    return val;
}

//should really change this!
double trap_k(function<double(double,double)> func,int n, double J){
    double h = (2.)/(100*n);
    double Trap = 0.0;
    double x,Trap_n = 0.0;
    for (int i=0;i<(100*n);i++){
        x = -1.+i*h;
        Trap = Trap + func(J,x);
        Trap_n = h*(func(J,-1)/2+Trap+func(J,1)/2);
    }
    return Trap_n;
    
}

double trap_f(function<double(double,double,double,double)> func,int n, double x_i, double x_j, double J){
    double h = (2.)/(100*n);
    double Trap = 0.0;
    double x,Trap_n = 0.0;
    for (int i=0;i<(100*n);i++){
        x = -1.+i*h;
        Trap = Trap + func(x_i,x_j,J,x);
        Trap_n = h*(func(x_i,x_j,J,-1)/2+Trap+func(x_i,x_j,J,1)/2);
    }
    return Trap_n;
    
}







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






/*
vector<vector<double>> Stiff_matrix(int n_elements, int s_elements, double start, double end){
    vector<vector<double>> Stiff(n_elements,vector<double>(n_elements));
    double h = (end-start)/n_elements;
    double J = h/2.;
    for (int k=0;k<n_elements-1;k++){ //this walks us through all of the stiffnesss matrix
        for (int i=k;i<k+s_elements;i++){
            for (int j=k;j<k+s_elements;j++){
                if (i==j){
                    Stiff[i][j] += trap_k(K1, n_elements, J);
                }
                else{
                    Stiff[i][j] += trap_k(K2, n_elements, J);
                }
            }
        }
    }
    
    //apply BC

    
    
    //Stiff[0][0] *= 2.;
    //Stiff[n_elements-1][n_elements-1] *= 2.;
    return Stiff;
}

vector<double> F_matrix(int n_elements, int s_elements, double start, double end){
    vector<double> F(n_elements);
    double h = (end-start)/n_elements;
    double J = h/2.;
    for (int k=0;k<n_elements-1;k++){
        double x_i = start+k*h;
        double x_j = x_i+h;
        F[k] += trap_f(F1, n_elements, x_i, x_j, J);
        F[k+1] += trap_f(F2, n_elements, x_i, x_j, J);
    }
    F[0] *= 2.;
    F[n_elements-1] *= 2.;
    
    return F;
}
*/


//Need to create a K local and F local matrix and then tack on the assembly matrix...



int main(int argc, const char * argv[]) {
    
    
    double start = 0.0;
    double end = 10;
    int n_elements = 100;
    int s_elements = 2;
    vector<vector<double>> K = Stiff_matrix(n_elements,s_elements,start,end);
    //cout<<K[1][2]<<endl;
    vector<double> F = F_matrix(n_elements, s_elements, start, end);
    //cout<<F[10]<<endl;
    //now to apply EBC
    vector<int> node_set = {0};
    vector<double> BC = {0.0};
    vector<int>wo_IC ={};
    for (int i=0;i<node_set.size();i++){
        if(std::find(node_set.begin(), node_set.end(), i) != node_set.end()) {
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
        for(int i=0;i<node_set.size();i++){
            int current_node_with = node_set[i];
            F_new[j] -= BC[current_node_with]*K[current_node_without][current_node_with];
            
        }
    }
    for (int kr=0;kr<wo_IC.size();kr++){
        int row = wo_IC[kr];
        for(int kc=0;kc<wo_IC.size();kc++){
            int col = wo_IC[kc];
            K_new[kr][kc] = K[row][col];
        }
    }

    
    
    vector<vector<double>> L,U;
    tuple<vector<vector<double>> ,vector<vector<double>>> LUs = LU(K, 1, n_elements);
    L = get<0>(LUs);
    U = get<1>(LUs);
    //cout<<U[10][99]<<endl;
    vector<double> c = solveUnitLowerTriangularBanded(L, F, 1, n_elements);
    vector<double> x_sol = solveUnitUpperTriangularbanded(U, c, 1, n_elements);
    
    vector<double> x_axis;
    double h = (end-start)/n_elements;
    
    vector<double> fin_sol;
    fin_sol.push_back(0.0);
    for (int i=0;i<x_sol.size();i++){
        fin_sol.push_back(x_sol[i]);
    }
    
    ofstream output;
    output.open("/Users/crhea/Desktop/FEM.txt");
    output<<"X_axis "<<"Value "<<endl;
    for(int k=0;k<n_elements;k++){
        output<<start+k*h<<" "<<fin_sol[k]<<endl;
    }
    output.close();
    
    
    
    return 0;
}
