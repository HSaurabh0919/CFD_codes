#include <stdio.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>


using namespace std;

double maxx(double a,double b)
{if(a>=b) return a;
else return b;
}

double abs(double a)
{
    if(a<0) a=a*(-1);
    return a;
}


vector<double> gauss(vector< vector<double> > A,double p,int m,int nn) {
    printf("5\n");
    int n = nn*m;double tmp,tmp2,tmp3;vector<double> result(n+1,0);int i,j;
    for(i=0;i<n;i++){result[i]=0;}
    double error=1.0;
    //int var=1;

	while(error>=0.01)
	{
	    error=0.001;
		for(i=0;i<n;i++)
		{
		 tmp=0;
		 for(j=0;j<n;j++)
			{
			 if(i!=j) {tmp+=A[i][j]*result[j];}
			}
			tmp2=result[i];
		    result[i]=result[i]*(1-p)+p*(A[i][n]-tmp)/A[i][i];
		    tmp3=abs(result[i]-tmp2);
            error=maxx(error,tmp3);

		 }

	}
	printf("6\n");
	return result;

}

int main()
{
    ofstream Output;
    Output.open("output.csv");
int i,j,ll,kkk;
int m=5;
int n=5;
double af=1.0;double w=1.0;double l=1.0;
double dx=l/m;double dy=w/n;
double dt=1.0;
int tamb=30;
double tc=500.0;
double afa=1.0;
double A=tc+afa;
double omg=1.0;
int h=1;
int k=1;
int t=5;   //time
int c=0;
printf("1\n");
int iii;
double ap=(1+(4*af*dt/(dx*dx)));
printf("A %lf %lf %lf %lf \n",af,dt,dx,(af*dt/(dx*dx)));
double an=-(af*dt/(dx*dx));
double aw=an;
double as=an;
double ae=an;
double Q[n][m][7];
double ccc=(1.0+(dx*h/(2*k)));
double hk=(1.0-(dx*h/(2*k)))/ccc;
double cont=2*af*dt*h*tamb/((2*k-dx*h)*(dx*dx));
printf("%lf \n",omg);
//%initialization
vector<double> T(m*n+1,0);
printf("%lf %lf %lf %lf \n",dx,ccc,hk,cont);
vector<double> line2(m*n+1,0);
vector< vector<double> > B(m*n,line2);
for(i=0;i<n;i++)
    {
        for(j=0;j<m;j++)
        {
            for(iii=0;iii<=6;iii++)
           {

            Q[i][j][iii]=0.0;
           }

        }

    }
printf("2 %lf %lf %lf %lf\n",aw,A,omg,an);

printf("aw*2*A*omg %lf \n",aw*2*A*omg);

for(kkk=0;kkk<=m*n;kkk++)
{
    T[kkk]=30.0;
}
printf("%lf \n",T[0*m+m-1]+tamb*h*dx/(k*ccc));
for(i=0;i<n;i++)
{


    for(j=0;j<m;j++)
    {


         if (i==0 && j==0)
            {Q[i][j][1]=ap;
             Q[i][j][3]=ae;
             Q[i][j][4]=an;
             Q[i][j][6]=T[i*m+j]- aw*2*A*omg;}

         else if (i==n-1 && j==0)
             {Q[i][j][1]=ap;
             Q[i][j][3]=ae;
             Q[i][j][5]=as;
             Q[i][j][6]=T[i*m+j]- aw*2*A*omg-an*dy*c;}

         else if(i==0 && j==m-1)
             {Q[i][j][1]=ap+ae*hk+as;
             Q[i][j][2]=aw;
             Q[i][j][4]=an;
             Q[i][j][6]=T[i*m+j]+tamb*h*dx/(k*ccc);
             }

         else if(i==n-1 && j==m-1)
             {Q[i][j][1]=ap+ae*hk+as+an;
             Q[i][j][2]=aw;
             Q[i][j][5]=as;
             Q[i][j][6]=T[i*m+j] + aw*dy*c+tamb*h*dx/(k*ccc);
             }

         else if( j==0)
             {Q[i][j][1]=ap-aw;
             Q[i][j][3]=ae;
             Q[i][j][5]=as;
             Q[i][j][4]=an;
             Q[i][j][6]=T[i*m+j]- aw*2*A*omg;
             }
         else if(i==n-1)
            {

             Q[i][j][1]=ap+aw;
             Q[i][j][3]=ae;
             Q[i][j][2]=aw;
             Q[i][j][5]=as;
             Q[i][j][6]=T[i*m+j]+aw*dy*c;
            }

         else if(j==m-1)
             {Q[i][j][1]=ap+ae*hk;
             Q[i][j][2]=aw;
             Q[i][j][5]=as;
             Q[i][j][4]=an;
             Q[i][j][6]=T[i*m+j]+tamb*h*dx/(k*ccc);
             }

         else if (i==0)
             {Q[i][j][1]=ap+as;
             Q[i][j][3]=ae;
             Q[i][j][2]=aw;
             Q[i][j][4]=an;
             Q[i][j][6]=T[i*m+j];
             }
          else
          {


             Q[i][j][1]=ap;
             Q[i][j][2]=aw;
             Q[i][j][3]=ae;
             Q[i][j][4]=an;
             Q[i][j][5]=as;
             Q[i][j][6]=T[i*m+j];
          }
         }
    }
printf("3\n");
    // second part

    //double C[m*n];
    for(i=0;i<n*m;i++)
    {
        for(j=0;j<=m*n;j++)
        {
            B[i][j]=0.0;
        }
    }

  Output<<"B5\n";

for(i=0;i<n;i++)
    {
        for(j=0;j<m;j++)
        {
            for(iii=1;iii<=6;iii++)
           {

            Output<<Q[i][j][iii]<<",";
           }
           Output<<"\n";
        }
        Output<<"\n";
    }

    i=0;
    j=0;
    int kk=0;
    printf("10\n");
    for(kk=0;kk<m*n;kk++)
    {
       if(i==0&&j==0)
       {
           B[kk][i*m+j]=Q[i][j][1];
        //B[kk][i*m+j-1]=Q[i][j][2];
        B[kk][i*m+j+1]=Q[i][j][3];
        B[kk][(i+1)*m+j]=Q[i][j][4];
        //B[kk][(i-1)*m+j]=Q[i][j][5];
        B[kk][m*n]=Q[i][j][6];
       }
       else if(i==n-1&&j==m-1)
       {
           B[kk][i*m+j]=Q[i][j][1];
        B[kk][i*m+j-1]=Q[i][j][2];
        //B[kk][i*m+j+1]=Q[i][j][3];
        //B[kk][(i+1)*m+j]=Q[i][j][4];
        B[kk][(i-1)*m+j]=Q[i][j][5];
        B[kk][m*n]=Q[i][j][6];
       }
       else if(i==n-1&&j==0)
       {
        B[kk][i*m+j]=Q[i][j][1];
       // B[kk][i*m+j-1]=Q[i][j][2];
        B[kk][i*m+j+1]=Q[i][j][3];
       // B[kk][(i+1)*m+j]=Q[i][j][4];
        B[kk][(i-1)*m+j]=Q[i][j][5];
        B[kk][m*n]=Q[i][j][6];
       }
       else if(i==0&&j==m-1)
       {
        B[kk][i*m+j]=Q[i][j][1];
        B[kk][i*m+j-1]=Q[i][j][2];
        //B[kk][i*m+j+1]=Q[i][j][3];
        B[kk][(i+1)*m+j]=Q[i][j][4];
       // B[kk][(i-1)*m+j]=Q[i][j][5];
        B[kk][m*n]=Q[i][j][6];
       }
       else if(i==0)
       {
        B[kk][i*m+j]=Q[i][j][1];
        B[kk][i*m+j-1]=Q[i][j][2];
        B[kk][i*m+j+1]=Q[i][j][3];
        B[kk][(i+1)*m+j]=Q[i][j][4];
        //B[kk][(i-1)*m+j]=Q[i][j][5];
        B[kk][m*n]=Q[i][j][6];
       }
       else if(i==n-1)
       {
           B[kk][i*m+j]=Q[i][j][1];
        B[kk][i*m+j-1]=Q[i][j][2];
        B[kk][i*m+j+1]=Q[i][j][3];
       // B[kk][(i+1)*m+j]=Q[i][j][4];
        B[kk][(i-1)*m+j]=Q[i][j][5];
        B[kk][m*n]=Q[i][j][6];
       }
       else if(j==0)
       {
           B[kk][i*m+j]=Q[i][j][1];
        //B[kk][i*m+j-1]=Q[i][j][2];
        B[kk][i*m+j+1]=Q[i][j][3];
        B[kk][(i+1)*m+j]=Q[i][j][4];
        B[kk][(i-1)*m+j]=Q[i][j][5];
        B[kk][m*n]=Q[i][j][6];
       }
       else if(j==m-1)
       {

           B[kk][i*m+j]=Q[i][j][1];
        B[kk][i*m+j-1]=Q[i][j][2];
        //B[kk][i*m+j+1]=Q[i][j][3];
        B[kk][(i+1)*m+j]=Q[i][j][4];
        B[kk][(i-1)*m+j]=Q[i][j][5];
        B[kk][m*n]=Q[i][j][6];
       }


        else{

        B[kk][i*m+j]=Q[i][j][1];
        B[kk][i*m+j-1]=Q[i][j][2];
        B[kk][i*m+j+1]=Q[i][j][3];
        B[kk][(i+1)*m+j]=Q[i][j][4];
        B[kk][(i-1)*m+j]=Q[i][j][5];
        B[kk][m*n]=Q[i][j][6];
        }


        if (j==m-1)
           {

            i=i+1;
            j=0;
           }
        else {j=j+1;}
       // printf("A %d %d \n",i,j);
    }
printf("B11 \n");
for(i=0;i<=m*n;i++)
{
    T[i]=0.0;
}
printf("12 \n");
/*

 for(i=0;i<m*n;i++)
    {
        for(j=0;j<=m*n;j++)
        {
            printf("%lf ",B[i][j]);
        }
        printf("\n");
    }
printf("B ENDS\n");
*/
    Output<<"\n";
    T=gauss(B,0.3,m,n);
    i=0;
    j=0;
    //int kk;
    for(kk=0;kk<m*n;kk++)
    {
       Output<<T[kk]<<",";
        if (j==m-1)
           {

            i=i+1;
            j=0;
             Output<<"\n";
           }
        else {j=j+1;}
    }
    printf("4\n");

    return 0;
    }










