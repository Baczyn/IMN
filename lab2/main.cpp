#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define Beta 0.001 //czestosc kontaktu chorych ze zdrowymi
#define N 500 //ilosc ludzi
#define Gamma 0.1 //czas trwania choroby
#define TOL 0.000001
#define dt 0.1 
#define M 20   //max iteracji
#define u_0 1 //poczatkowa liczba chorych
#define t_max 100      //

//wspolczynik do zadania 3
double a11 = 0.25;
double a22 = a11;
double a12 = 0.25 - sqrt(3)/6.0;
double a21 = 0.25 + sqrt(3)/6.0;


double alfa(){return Beta*N - Gamma;}
//zad1//
double f(double u){
    return alfa()*u - Beta*pow(u,2);
}
void picard(ofstream &file){
    double u_n=u_0;
    double u_next=u_0;
    
    for(double i=0;i<t_max;i+=dt){
        u_n=u_next;
        double u_mi=0;
        int iter=0;
        while( (abs(u_next-u_mi)>TOL) && (iter<M)){
            u_mi=u_next;
            u_next=u_n+(dt/2.0)*(f(u_n)+f(u_mi));
            iter++;
        }  
        file<<i<<" "<<u_next<<" "<<N-u_next<<"\n"; 
        cout<<i<<" "<<iter<<"\n";

    }
    cout<<"\n\n";
}
//zad2//
double df(double u){
    return alfa() - Beta*2.0*u;
}
void newton(ofstream &file){
    double u_n=u_0;
    double u_next=u_0;
    
    for(double i=0;i<t_max;i+=dt){
        u_n=u_next;
        double u_mi=0;
        int iter=0;
        while( (abs(u_next-u_mi)>TOL) && (iter<M)){
            u_mi=u_next;
            u_next =u_mi-(u_mi-u_n-(dt/2.0)*(f(u_n)+f(u_mi)))/(1.0-(dt/2.0)*df(u_mi));
            iter++;
        }  
        file<<i<<" "<<u_next<<" "<<N-u_next<<"\n"; 
        cout<<i<<" "<<iter<<"\n";
    }
    cout<<"\n\n";
}
//zad3//
double F1(double U1, double U2, double u_n){
    return (U1-u_n-dt*(a11*(alfa()*U1-Beta*pow(U1,2)) + a12*(alfa()*U2 - Beta*pow(U2,2))));
}
double F2(double U1, double U2,double u_n){
    return (U2-u_n-dt*(a21*(alfa()*U1-Beta*pow(U1,2)) + a22*(alfa()*U2 - Beta*pow(U2,2))));
}
double m11(double U1){
    return (1-dt*a11*(alfa()-2.0*Beta*U1));
}

double m12(double U2){
    return ((-1)*dt*a12*(alfa() - 2.0*Beta*U2));
}

double m21(double U1){
    return ((-1)*dt*a21*(alfa() - 2.0*Beta*U1));
}

double m22(double U2){
    return (1-dt*a22*(alfa() - 2.0*Beta*U2));
}

void RK2(ofstream &file){
    double U1,U2,deltaU1,deltaU2;
    double u_n=u_0;
    double u_next=u_0;

    for(double i=0;i<t_max;i+=dt){
        u_n=u_next;
        double u_mi=0;
        int iter=0;
        U1=U2=u_n;
        while( (abs(u_next-u_mi)>TOL) && (iter<M)){
            u_mi=u_next;
            deltaU1=(F2(U1,U2,u_n)*m12(U2)-F1(U1,U2,u_n)*m22(U2))/
                    (m11(U1)*m22(U2)-m12(U2)*m21(U1));
            deltaU2=(F1(U1,U2,u_n)*m21(U1)-F2(U1,U2,u_n)*m11(U1))/
                    (m11(U1)*m22(U2)-m12(U2)*m21(U1));
            
            U1=deltaU1+U1;
            U2=deltaU2+U2;

            u_next = u_n + dt*(0.5*f(U1) + 0.5*f(U2));
            iter++;
        }
    file<<i<<" "<<u_next<<" "<<N-u_next<<"\n"; 
    cout<<i<<" "<<iter<<"\n";
   }
   cout<<"\n\n";
}


int main(){
    ofstream file1, file2, file3;
    file1.open("picard.txt");
    file2.open("newton.txt");
    file3.open("RK2.txt");
    picard(file1);
    newton(file2);
    RK2(file3);

}
