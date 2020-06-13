#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <fstream>

using namespace std;



vector<double> gauss(vector< vector<double> > A) {
    int n = A.size();

    for (int i=0; i<n; i++) {
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }


        for (int k=i; k<n+1;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        for (int k=i+1; k<n; k++) {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }


    vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}

int main() {
    ofstream Output;
    Output.open("output.csv");
    int n,i,j,ii,tt;double z,a,b,c,d,dt,dx,alpha,u,Q,k,h,Ta;
    printf("Enter Value of n\n");
    scanf("%d",&n);
    vector<double> line(n+1,0);
    vector< vector<double> > mat(n,line);
    vector<double> T(n+1,0);
    printf("Enter Value of ambient temperature\n");
    scanf("%lf",&Ta);
    printf("Enter Value of dt\n");
    scanf("%lf",&dt);
    printf("Enter Value of total time\n");
    scanf("%d",&tt);
    printf("Enter Value of dx\n");
    scanf("%lf",&dx);
    printf("Enter Value of alpha\n");
    scanf("%lf",&alpha);
    printf("Enter Value of u\n");
    scanf("%lf",&u);
    printf("Enter Value of Q\n");
    scanf("%lf",&Q);
    printf("Enter Value of  k\n");
    scanf("%lf",&k);
    printf("Enter Value of h\n");
    scanf("%lf",&h);

    Output<<"n="<<n<<","<<"Ta="<<Ta<<","<<"dt="<<dt<<","<<"tt="<<tt<<","<<"k="<<k<<",\n";
   Output<<"dx="<<dx<<","<<"alpha="<<alpha<<","<<"u="<<u<<","<<"Q="<<Q<<","<<"h="<<h<<",\n";
    for(i=0;i<n;i++)
    {
            T[i]=Ta;
    }
    for(ii=dt;ii<tt;ii++)
    {

    for(i=0;i<n;i++)
    {
        if(i==0)
        {
            z=(1+(alpha*dt/(dx*dx))+0*u*dt/dx);
            mat[i][n]=T[i]+(1*Q*u*dt/k)+(alpha*dt*Q/(dx*k));
            mat[i][0]=z;
            mat[i][1]=-1*(alpha*dt/(dx*dx));
            for(j=2;j<n;j++)
        {
            mat[i][j]=0;
        }
        }

        else if(i==n-1)
        {

            mat[i][n]=T[i]+h*Ta*alpha*dt/(k*dx*dx);
            mat[i][n-1]=(1+alpha*dt*h/(k*dx)+1.0*u*dt/dx+1*alpha*dt/(dx*dx));
            mat[i][n-2]=(-1*u*dt/dx-alpha*dt/(dx*dx));
            for(j=0;j<n-2;j++)
        {
            mat[i][j]=0;
        }
        }
        else
        {
            z=(1+2*alpha*dt/(dx*dx)+1.0*u*dt/dx);
             for(j=0;j<n;j++)
            {
                if(j==i-1)
                {

                    mat[i][j]=-1*(u*dt/dx+alpha*dt/(dx*dx));
                    mat[i][j+1]=1*z;
                    mat[i][j+2]=-1*(alpha*dt/(dx*dx));
                    j+=2;
                }
                else    {mat[i][j]=0;}
            }
            mat[i][n]=T[i];

        }

    }


    T = gauss(mat);
    for (int i=0; i<n; i++) {
        Output << T[i] << ",";
    }
   Output << endl;
    }
    Output.close();
    return 0;
}

