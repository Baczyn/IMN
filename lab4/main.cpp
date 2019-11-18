#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

#define epsilon 1
#define delta 0.1
#define n_x 150
#define n_y 100
#define V1 10
#define TOL pow(10,-5)

double X_max = delta*n_x;
double Y_max = delta*n_y;
double sigma_X = 0.1*X_max;
double sigma_Y= 0.1*Y_max;







double gestosc_q(double x,double y){
    return exp(-(pow(x-0.35*X_max,2)/pow(sigma_X,2)) - (pow(y-0.5*Y_max,2)/pow(sigma_Y,2) )) -
           exp(-(pow(x-0.65*X_max,2)/pow(sigma_X,2)) - (pow(y-0.5*Y_max,2)/pow(sigma_Y,2) ) );


}
double S(double tab[][n_y + 1]){
double s = 0;
    for(int i = 0; i < n_x; i++) 
		for(int j = 0; j < n_y; j++) 
			s += pow(delta,2) * ( 0.5 * pow( ( tab[i+1][j] - tab[i][j])/delta, 2)
              +0.5 * pow( (tab[i][j+1] - tab[i][j])/delta, 2) - gestosc_q(i*delta,j*delta)*tab[i][j]);

    return s;
}

void global(double omega,ofstream &file,ofstream &file2){
    double Sn=0;//nowe
    double Ss=0;//stare
    double Vn[n_x+1][n_y+1] = {0};
    double Vs[n_x+1][n_y+1] = {0};

    for(int i =0;i<n_x+1;i++){
        Vn[i][0]=V1;
        Vs[i][0]=V1;
    }

    int it=0;
    do{
        it++;
        for(int i=1;i<n_x;i++){
            for(int j=1;j<n_y;j++){
                    Vn[i][j] = 0.25 * ( Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1]
                    +  delta*delta/epsilon * gestosc_q(i*delta,j*delta) );
                }
            }
        for(int j = 1; j < n_y; j++){
            Vn[0][j] = Vn[1][j];
            Vn[n_x][j] = Vn[n_x - 1][j];
            }

        for(int i = 0 ; i < n_x+1 ; i++){
            for(int j = 0; j < n_y+1; j++){
                Vs[i][j] = (1-omega)*Vs[i][j] + omega*Vn[i][j];
            }
        }

        Ss = Sn;
        Sn = S(Vs);
        file<<it<<" "<<Ss<<"\n";
    }while(fabs( (Sn - Ss)/Ss ) > TOL);
    
    for(int i=0;i<n_x+1;i++){
    for(int j=0;j<n_y+1;j++){
        file2<<i*delta<<" "<<j*delta<<" "<<Vn[i][j]<<" "<< ( Vs[i+1][j] -2*Vs[i][j] + Vs[i-1][j])/(delta*delta) +
                    (Vs[i][j+1] -2*Vs[i][j] + Vs[i][j-1])/(delta*delta) + gestosc_q(i*delta,j*delta)/epsilon <<"\n";
        file2<<"\n";
        }
    }
}

void local(double omega,  ofstream &file,ofstream &file2){

    double Sn=0;//nowe
    double Ss=0;//stare
    double V[n_x+1][n_y+1] = {0};

    for(int i =0;i<n_x+1;i++){
        V[i][0]=V1;
    }

    int it=0;
    do{
        it++;

        for(int i = 1; i < n_x; i++){
            for(int j = 1; j < n_y; j++){
                V[i][j] = (1 - omega)* V[i][j] +(omega/4) * ( V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1]
                    +  delta*delta/epsilon * gestosc_q(i*delta,j*delta) );
                }
        }
        for(int j = 1; j < n_y; j++){
            V[0][j] = V[1][j];
            V[n_x][j] = V[n_x - 1][j];
            }

        Ss = Sn;
        Sn = S(V);
        file<<it<<" "<<Ss<<"\n";
    }while(fabs( (Sn - Ss)/Ss ) > TOL);

    for(int i=0;i<n_x+1;i++){
        for(int j=0;j<n_y+1;j++){
        file2<<i*delta<<" "<<j*delta<<" "<<V[i][j]<<"\n";
            file2<<"\n";
        }
}

}



int main(){
    ofstream file1, file2,file3,file4,file5,file6,file7,file8,file9,file10,file11,file12;
    file1.open("calka_global_06.txt");
    file2.open("Vglobal06.txt");
    file3.open("calka_global_1.txt");
    file4.open("Vglobal1.txt");
    global(0.6,file1,file2);
    global(1,file3,file4);
    

     file5.open("calka_local1.txt");
     file6.open("Vlocal1.txt");
     local(1,file5,file6);

     file7.open("calka_local1_4.txt");
     file8.open("Vlocal1_4.txt");
     local(1.4,file7,file8);

     file9.open("calka_local1_8.txt");
     file10.open("Vlocal1_8.txt");
     local(1.8,file9,file10);

     file11.open("calka_local1_9.txt");
     file12.open("Vlocal1_9.txt");
     local(1.9,file11,file12);




}
