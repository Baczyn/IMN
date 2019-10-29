#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void jawny_euler(double step, int t_max, int t_min,ofstream &file){
    int lb=-1;
    double y0=1;
    double y;
    double y_dokladne;
    for(double i=t_min;i<t_max; i+=step){
        y_dokladne=exp(lb*i);
        y = y0 + step*lb*y0;
        y0=y;
        file<< i << " " << y << " "<<y-y_dokladne <<"\n";
    }   
}

void jawny_RK2(double step, int t_max, int t_min,ofstream &file){
    
    int lb=-1;
    double y0=1;
    double y=1;
    double k1;
    double k2;
    double y_dokladne;
    for(double i=t_min;i<t_max; i+=step){
        y_dokladne=exp(lb*(i+step));

        k1=lb*y;
        k2=lb*(y+step*k1);

        y = y0 + (step/2.0)*(k1+k2);
        y0=y;

        file<< i +step<< " " << y << " "<<y-y_dokladne <<"\n";
        
    }

}

void jawny_RK4(double step, int t_max, int t_min,ofstream &file){
        
    int lb=-1;
    double y0=1;
    double y=1;
    double k1;
    double k2;
    double k3;
    double k4;


    double y_dokladne;
    for(double i=t_min;i<t_max; i+=step){
        y_dokladne=exp(lb*(i+step));
        k1=lb*y;
        k2=lb*(y+(step/2.0)*k1);
        k3=lb*(y+(step/2.0)*k2);
        k4=lb*(y+step*k3);

        y=y0+(step/6.0)*(k1+2*k2+2*k3+k4);
        y0=y;

        file<< i +step<< " " << y << " "<<y-y_dokladne <<"\n";
       


}}

double V(double W,double t){return 10*sin(W*t);}

double f(double t,double Q,double I){return I;}

double g(double t,double Q,double I,double W,double L,double R,double C){
    return (V(W,t)/L) - R/(L*I) -Q/(L*C);
}

void RK4_RZII(double step, fstream &file,double R,double C,double L,double f){
double k1_q,k1_i,k2_q,k2_i,k3_q,k3_i,k4_q,k4_i;

double q0=0;
double i0=0;
double Q=q0;

double I=i0;

double W0=1/sqrt(L*C);
double WV=f*W0;












}


int main(){
    ///////zad1//////////////
    ofstream file1, file2, file3;
    file1.open("wynik_euler1");
    file2.open("wynik_euler2");
    file3.open("wynik_euler3");

    jawny_euler(1,5,0,file1);
    jawny_euler(0.1,5,0,file2);
    jawny_euler(0.01,5,0,file3);
    ///////////zad2//////////////////
    ofstream file11, file22, file33;
    file11.open("RK2_1");
    file22.open("RK2_2");
    file33.open("RK2_3");

    jawny_RK2(1,5,0,file11);
    jawny_RK2(0.1,5,0,file22);
    jawny_RK2(0.01,5,0,file33);
        ///////////zad3//////////////////
    ofstream file111, file222, file333;
    file111.open("RK4_1");
    file222.open("RK4_2");
    file333.open("RK4_3");

    jawny_RK4(1,5,0,file111);
    jawny_RK4(0.1,5,0,file222);
    jawny_RK4(0.01,5,0,file333);
    //////////zad4///////////


}