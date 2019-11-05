#include <iostream>
#include <cmath>
#include <utility> //dla stdpair
#include <functional>
#include <fstream>

#define x0 0.01
#define v0 0.0
#define dt_0 1.0
#define S 0.75
#define p 2.0
#define alfa 5
#define t_max 40


using namespace std;

using pair_d=pair<double,double>;

double g(pair_d xv){
    return alfa*(1-pow(xv.first,2))*xv.second - xv.first;
}

double a11(){
    return 1;
}

double a12(double dt){
    return (-1)*(dt/2.0);
}

double a21(double dt,pair_d xv_k){
    return (-1)*(dt/2.0)*((-2.0)*alfa*xv_k.first*xv_k.second-1.0);
}

double a22(double dt,pair_d xv_k){
    return (1 - (dt/2.0)*alfa*(1-pow(xv_k.first,2)) );
}

double F(pair_d xv,pair_d xv1, double dt){

    return xv1.first - xv.first - (dt/2.0)*(xv.second+xv1.second);
}

double G(pair_d xv,pair_d xv1, double dt){

    return xv1.second - xv.second - (dt/2.0)*(g(xv)+g(xv1));
}


pair_d trapez(pair_d xv,double dt){
    double sigma= pow(10,-10);
    double dx=100,dv=110;
    pair_d xv_n(xv);

    while( fabs(dx) > sigma || fabs(dv) > sigma){
       

         dx = ((-1)*F(xv,xv_n, dt)*a22(dt,xv_n)-(-1)*G(xv,xv_n,dt)*a12(dt))/
              (a11()*a22(dt,xv_n) - a12(dt)*a21(dt,xv_n));

         dv = (a11()*(-1)*G(xv,xv_n,dt) - a21(dt,xv_n)*(-1)*F(xv,xv_n,dt))/
              (a11()*a22(dt,xv_n) - a12(dt)*a21(dt,xv_n));

    xv_n=pair_d(xv_n.first+dx,xv_n.second+dv);

    }

    return xv_n;

}



//wzor 25-31
pair_d RK2(pair_d xv,double dt){
    pair_d k1(xv.second,    g(xv)); 
    pair_d k2(xv.second+dt*k1.second,  g(pair_d(xv.first+ dt*k1.first,xv.second +dt*k1.second)) );
    return pair_d(xv.first+ (dt/2.0)*(k1.first+k2.first),xv.second+ (dt/2.0)*(k1.second+k2.second));

}

void kontrola(function<pair_d(pair_d,double)> schemat_num,double TOL,ofstream &file){
    double t=0; double dt=dt_0; pair_d xv(x0,v0);
    pair_d xv_1;
    pair_d xv_2;
    double Ex,Ev;

    do{
        xv_2 = schemat_num(xv,dt);
        xv_2 = schemat_num(xv_2,dt);

        xv_1 = schemat_num(xv,2.0*dt);

        Ex = ( xv_2.first - xv_1.first )/(pow(2,p)-1.0);
        Ev = ( xv_2.second - xv_1.second)/(pow(2,p)-1.0);

        if( max(fabs(Ex),fabs(Ev)) < TOL ){
            t=t+2*dt;
            xv = pair_d(xv_2);
            file<<t<<" "<<dt<<" "<<xv.first<<" "<<xv.second<<"\n";
        }
        dt= pow(((S*TOL)/(max(fabs(Ex),fabs(Ev)) )), (1.0/(p+1.0)))*dt;
        }while(t<t_max);

    
       
}

int main(){
    ofstream file1, file2;
    file1.open("RK2.txt");
    file2.open("trapez.txt");

    kontrola(RK2,pow(10,-2),file1);
    kontrola(RK2,pow(10,-5),file1);

    kontrola(trapez,pow(10,-2),file2);
    kontrola(trapez,pow(10,-5),file2);

    file1.close();
    file2.close();

}
