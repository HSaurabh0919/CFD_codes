#include <iostream>

#include<bits/stdc++.h>
using namespace std;

#define N 100        // Number of unknowns


// function to reduce matrix to r.e.f.  Returns a value to
// indicate whether matrix is singular or not
int forwardElim(double mat[N][N+1],int n);

// function to calculate the values of the unknowns
void backSub(double mat[N][N+1],int n);

// function to get matrix content
void gaussianElimination(double mat[N][N+1],int n)
{
    /* reduction into r.e.f. */
    int singular_flag = forwardElim(mat,n);

    /* if matrix is singular */
    if (singular_flag != -1)
    {
        printf("Singular Matrix.\n");

        /* if the RHS of equation corresponding to
           zero row  is 0, * system has infinitely
           many solutions, else inconsistent*/
        if (mat[singular_flag][n])
            printf("Inconsistent System.");
        else
            printf("May have infinitely many "
                   "solutions.");

        return;
    }

    /* get solution to system and print it using
       backward substitution */
    backSub(mat,n);
}

// function for elemntary operation of swapping two rows
void swap_row(double mat[N][N+1], int i, int j,int n)
{
    //printf("Swapped rows %d and %d\n", i, j);

    for (int k=0; k<=n; k++)
    {
        double temp = mat[i][k];
        mat[i][k] = mat[j][k];
        mat[j][k] = temp;
    }
}

// function to print matrix content at any stage
void print(double mat[N][N+1],int n)
{
    for (int i=0; i<n; i++, printf("\n"))
        for (int j=0; j<=n; j++)
            printf("%lf ", mat[i][j]);

    printf("\n");
}

// function to reduce matrix to r.e.f.
int forwardElim(double mat[N][N+1],int n)
{
    for (int k=0; k<n; k++)
    {
        // Initialize maximum value and index for pivot
        int i_max = k;
        int v_max = mat[i_max][k];

        /* find greater amplitude for pivot if any */
        for (int i = k+1; i < n; i++)
            if (abs(mat[i][k]) > v_max)
                v_max = mat[i][k], i_max = i;

        /* if a prinicipal diagonal element  is zero,
         * it denotes that matrix is singular, and
         * will lead to a division-by-zero later. */
        if (!mat[k][i_max])
            return k; // Matrix is singular

        /* Swap the greatest value row with current row */
        if (i_max != k)
            swap_row(mat, k, i_max,n);


        for (int i=k+1; i<n; i++)
        {
            /* factor f to set current row kth elemnt to 0,
             * and subsequently remaining kth column to 0 */
            double f = mat[i][k]/mat[k][k];

            /* subtract fth multiple of corresponding kth
               row element*/
            for (int j=k+1; j<=n; j++)
                mat[i][j] -= mat[k][j]*f;

            /* filling lower triangular matrix with zeros*/
            mat[i][k] = 0;
        }

        //print(mat);        //for matrix state
    }
    //print(mat);            //for matrix state
    return -1;
}

// function to calculate the values of the unknowns
void backSub(double mat[N][N+1],int n)
{
    double x[n];  // An array to store solution

    /* Start calculating from last equation up to the
       first */
    for (int i = n-1; i >= 0; i--)
    {
        /* start with the RHS of the equation */
        x[i] = mat[i][n];

        /* Initialize j to i+1 since matrix is upper
           triangular*/
        for (int j=i+1; j<n; j++)
        {
            /* subtract all the lhs values
             * except the coefficient of the variable
             * whose value is being calculated */
            x[i] -= mat[i][j]*x[j];
        }

        /* divide the RHS by the coefficient of the
           unknown being calculated */
        x[i] = x[i]/mat[i][i];
    }

    printf("\nSolution for the system:\n");
    for (int i=0; i<n; i++)
        printf("%lf\n", x[i]);
}

// Driver program
int main()
{
    /* input matrix */
    int n,i,j;double z,a,b,c,d,dt,dx,alpha,u,Q,k,h,Ta;
    printf("Enter Value of n\n");
    scanf("%d",&n);
    double mat[N][N+1];double T[n];
    printf("Enter Value of ambient temperature\n");
    scanf("%lf",&Ta);
    printf("Enter Value of dt\n");
    scanf("%lf",&dt);
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

    for(i=0;i<n;i++)
    {
            T[i]=Ta;
    }
    for(i=0;i<n;i++)
    {
        if(i==0)
        {
            z=(1+alpha*dt/(dx*dx));
            mat[i][n]=T[i]+(2*Q*u*dt/k)+(alpha*dt*Q/(dx*k));
            mat[i][0]=z;
            mat[i][1]=-1*(alpha*dt/(dx*dx));
            for(j=2;j<n;j++)
        {
            mat[i][j]=0;
        }
        }
        else if(i==1)
        {
            z=(1+2*alpha*dt/(dx*dx)+1.5*u*dt/dx);
            mat[i][0]=-(alpha*dt/(dx*dx)-2*u*dt/dx);
            mat[i][1]=z;
            mat[i][2]=-1*alpha*dt/(dx*dx);
            mat[i][n]=T[i]-Q*u*dt/k;
            for(j=3;j<n;j++)
        {
            mat[i][j]=0;
        }
        }
        else if(i==n-1)
        {
            mat[i][n]=T[i]+alpha*dt*h*Ta/(k*dx);
            mat[i][n-1]=(1+alpha*dt*h/(k*dx)+1.5*u*dt/dx+alpha*dt/(dx*dx));
            mat[i][n-2]=(-2*u*dt/dx-alpha*dt/(dx*dx));
            mat[i][n-3]=(0.5*u*dt/dx);
            for(j=0;j<n-3;j++)
        {
            mat[i][j]=0;
        }
        }
        else
        {
            z=(1+2*alpha*dt/(dx*dx)+3*u*dt/dx);
             for(j=0;j<n;j++)
            {
                if(j==i-2)
                {
                    mat[i][j]=0.5*u*dt/dx;
                    mat[i][j+1]=-1*(2*u*dt/dx+alpha*dt/(dx*dx));
                    mat[i][j+2]=1*z;
                    mat[i][j+3]=-1*(alpha*dt/(dx*dx));
                    j+=3;
                }
                else    {mat[i][j]=0;}
            }
            mat[i][n]=T[i];

        }

    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<=n;j++)
        {
            printf("%lf ",mat[i][j]);
        }
        printf("\n");
    }



    gaussianElimination(mat,n);

    return 0;
}
