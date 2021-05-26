
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "datadef.h"
#include "init.h"
#include <mpi.h>
#include <omp.h>



#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

extern int nprocs, proc;
extern int iLeft, iRight, eachNode;

MPI_Status status;

/* Computation of tentative velocity field (f, g) */
void compute_tentative_velocity(float **u, float **v, float **f, float **g,
    char **flag, int imax, int jmax, float del_t, float delx, float dely,
    float gamma, float Re)
{
    int  i, j;
    float du2dx, duvdy, duvdx, dv2dy, laplu, laplv;


    int iRightTemp = iRight;

    if (proc == nprocs - 1){
        iRightTemp -= 1;
    }


    #pragma omp parallel for schedule(static) private(du2dx,duvdy,laplu)
    for (i = iLeft; i <= iRightTemp; i++) {
        for (j=1; j<=jmax; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                du2dx = ((u[i][j]+u[i+1][j])*(u[i][j]+u[i+1][j])+
                    gamma*fabs(u[i][j]+u[i+1][j])*(u[i][j]-u[i+1][j])-
                    (u[i-1][j]+u[i][j])*(u[i-1][j]+u[i][j])-
                    gamma*fabs(u[i-1][j]+u[i][j])*(u[i-1][j]-u[i][j]))
                    /(4.0*delx);
                duvdy = ((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1])+
                    gamma*fabs(v[i][j]+v[i+1][j])*(u[i][j]-u[i][j+1])-
                    (v[i][j-1]+v[i+1][j-1])*(u[i][j-1]+u[i][j])-
                    gamma*fabs(v[i][j-1]+v[i+1][j-1])*(u[i][j-1]-u[i][j]))
                    /(4.0*dely);
                laplu = (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/delx/delx+
                    (u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dely/dely;

                f[i][j] = u[i][j]+del_t*(laplu/Re-du2dx-duvdy);
            } else {
                f[i][j] = u[i][j];
            }
        }
    }

    #pragma omp parallel for schedule(static) private(duvdx,dv2dy,laplv)
    for (i = iLeft; i <= iRight; i++) {
        for (j=1; j<=jmax-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                duvdx = ((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])+
                    gamma*fabs(u[i][j]+u[i][j+1])*(v[i][j]-v[i+1][j])-
                    (u[i-1][j]+u[i-1][j+1])*(v[i-1][j]+v[i][j])-
                    gamma*fabs(u[i-1][j]+u[i-1][j+1])*(v[i-1][j]-v[i][j]))
                    /(4.0*delx);
                dv2dy = ((v[i][j]+v[i][j+1])*(v[i][j]+v[i][j+1])+
                    gamma*fabs(v[i][j]+v[i][j+1])*(v[i][j]-v[i][j+1])-
                    (v[i][j-1]+v[i][j])*(v[i][j-1]+v[i][j])-
                    gamma*fabs(v[i][j-1]+v[i][j])*(v[i][j-1]-v[i][j]))
                    /(4.0*dely);

                laplv = (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/delx/delx+
                    (v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dely/dely;

                g[i][j] = v[i][j]+del_t*(laplv/Re-duvdx-dv2dy);
            } else {
                g[i][j] = v[i][j];
            }
        }
    }

    /* f & g at external boundaries */
    #pragma omp parallel for schedule(static) private(j)
    for (j=1; j<=jmax; j++) {
        f[0][j]    = u[0][j];
        f[imax][j] = u[imax][j];
    }
    
    #pragma omp parallel for schedule(static) private(i)
    for (i=iLeft; i<=iRight; i++) {
        g[i][0]    = v[i][0];
        g[i][jmax] = v[i][jmax];
    }

    // MPI_Barrier(MPI_COMM_WORLD);

    if (proc != 0) {
        MPI_Recv(&f[iLeft-1][0], jmax+2, MPI_FLOAT, proc-1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&g[iLeft-1][0], jmax+2, MPI_FLOAT, proc-1, 0, MPI_COMM_WORLD, &status);

        MPI_Send(&f[iLeft][0], jmax+2, MPI_FLOAT, proc-1, 0, MPI_COMM_WORLD);
        MPI_Send(&g[iLeft][0], jmax+2, MPI_FLOAT, proc-1, 0, MPI_COMM_WORLD);
    }
    if (proc != nprocs - 1) {
        MPI_Send(&f[iRight][0], jmax+2, MPI_FLOAT, proc+1, 0, MPI_COMM_WORLD);
        MPI_Send(&g[iRight][0], jmax+2, MPI_FLOAT, proc+1, 0, MPI_COMM_WORLD);

        MPI_Recv(&f[iRight+1][0], jmax+2, MPI_FLOAT, proc+1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&g[iRight+1][0], jmax+2, MPI_FLOAT, proc+1, 0, MPI_COMM_WORLD, &status);
    }
}


/* Calculate the right hand side of the pressure equation */
void compute_rhs(float **f, float **g, float **rhs, char **flag, int imax,
    int jmax, float del_t, float delx, float dely)
{
    int i, j;
    #pragma omp parallel for schedule(static) private(i,j)
    for (i=iLeft;i<=iRight;i++) {
        for (j=1;j<=jmax;j++) {
            if (flag[i][j] & C_F) {
                /* only for fluid and non-surface cells */
                rhs[i][j] = (
                             (f[i][j]-f[i-1][j])/delx +
                             (g[i][j]-g[i][j-1])/dely
                            ) / del_t;
            }
        }
    }
}
/* Red/Black SOR to solve the poisson equation */
int poisson(float **p, float **rhs, char **flag, int imax, int jmax,
    float delx, float dely, float eps, int itermax, float omega,
    float *res, int ifull)
{
    int i, j, iter;
    float add, beta_2, beta_mod;
    float p0 = 0.0;

    int rb; /* Red-black value. */

    float rdx2 = 1.0/(delx*delx);
    float rdy2 = 1.0/(dely*dely);
    beta_2 = -omega/(2.0*(rdx2+rdy2));

    float p0_Sum = 0.0;   // used to reduction later

    /* Calculate sum of squares */
    #pragma omp parallel for schedule(static) private(i,j) reduction(+:p0)
    for (i = iLeft; i <= iRight; i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[i][j] & C_F) { p0 += p[i][j]*p[i][j]; }
        }
    }

    // Reduce p0_Sum to p0 (summation)
    MPI_Allreduce(&p0, &p0_Sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    p0 = sqrt(p0_Sum/ifull);

    if (p0 < 0.0001) { p0 = 1.0; }

    /* Red/Black SOR-iteration */
    for (iter = 0; iter < itermax; iter++) {
        for (rb = 0; rb <= 1; rb++) {
          	#pragma omp parallel for schedule(static) private(beta_mod)
            for (i = iLeft; i <= iRight; i++) {
                for (j = 1; j <= jmax; j++) {
                    if ((i+j) % 2 != rb) { continue; }
                    if (flag[i][j] == (C_F | B_NSEW)) {
                        /* five point star for interior fluid cells */
                        p[i][j] = (1.-omega)*p[i][j] -
                              beta_2*(
                                    (p[i+1][j]+p[i-1][j])*rdx2
                                  + (p[i][j+1]+p[i][j-1])*rdy2
                                  -  rhs[i][j]
                              );
                    } else if (flag[i][j] & C_F) {
                        /* modified star near boundary */
                        beta_mod = -omega/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2);
                        p[i][j] = (1.-omega)*p[i][j] -
                            beta_mod*(
                                  (eps_E*p[i+1][j]+eps_W*p[i-1][j])*rdx2
                                + (eps_N*p[i][j+1]+eps_S*p[i][j-1])*rdy2
                                - rhs[i][j]
                            );
                    }
                } /* end of j */
            } /* end of i */

        if (proc != 0) {
      // Left
          MPI_Recv(&p[iLeft - 1][0], jmax + 2, MPI_FLOAT, proc - 1, 0, MPI_COMM_WORLD, &status);
          MPI_Send(&p[iLeft][0], jmax + 2, MPI_FLOAT, proc - 1, 0, MPI_COMM_WORLD);
      }
        if (proc != nprocs - 1) {
      // Right
          MPI_Send(&p[iRight][0], jmax + 2, MPI_FLOAT, proc + 1, 0, MPI_COMM_WORLD);
          MPI_Recv(&p[iRight + 1][0], jmax + 2, MPI_FLOAT, proc + 1, 0, MPI_COMM_WORLD, &status);
        }
        } /* end of rb */



        /* Partial computation of residual */
        *res = 0.0;
        float addSquared = 0.0;
        float addSquared_Sum = 0.0;
        for (i = iLeft; i <= iRight;i++) {
            for (j = 1; j <= jmax; j++) {
                if (flag[i][j] & C_F) {
                    /* only fluid cells */
                    add = (eps_E*(p[i+1][j]-p[i][j]) -
                        eps_W*(p[i][j]-p[i-1][j])) * rdx2  +
                        (eps_N*(p[i][j+1]-p[i][j]) -
                        eps_S*(p[i][j]-p[i][j-1])) * rdy2  -  rhs[i][j];

                    addSquared += add * add;
             }
            }
        }
        MPI_Allreduce(&addSquared, &addSquared_Sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        *res = sqrt((addSquared_Sum)/ifull)/p0;

        /* convergence? */
        if (*res<eps) break;
    } /* end of iter */
    return iter;
}


/* Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
 void update_velocity(float **u, float **v, float **f, float **g, float **p,
     char **flag,  int imax, int jmax, float del_t, float delx, float dely)
 {
     int i, j;



     int iRightTemp = iRight;
     
     // Account for the last region
     if (iRightTemp == imax) {
         iRightTemp -= 1;
     }
     #pragma omp parallel for schedule(static) private(i, j)
     for (i=iLeft; i<=iRightTemp; i++) {
         for (j=1; j<=jmax; j++) {
             /* only if both adjacent cells are fluid cells */
             if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                 u[i][j] = f[i][j]-(p[i+1][j]-p[i][j])*del_t/delx;
             }
         }
     }
     #pragma omp parallel for schedule(static) private(i, j)
     for (i=iLeft; i<=iRight; i++) {
         for (j=1; j<=jmax-1; j++) {
             /* only if both adjacent cells are fluid cells */
             if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                 v[i][j] = g[i][j]-(p[i][j+1]-p[i][j])*del_t/dely;
             }
         }
     }

     // MPI_Barrier(MPI_COMM_WORLD);

     if (proc != nprocs - 1) {
     // Pass to right
         MPI_Send(&u[iRight][0], jmax + 2, MPI_FLOAT, proc + 1, 0, MPI_COMM_WORLD);
         MPI_Recv(&u[iRight + 1][0], jmax + 2, MPI_FLOAT, proc + 1, 0, MPI_COMM_WORLD, &status);
       } if (proc != 0) {
     // Pass to left
         MPI_Recv(&u[iLeft - 1][0], jmax + 2, MPI_FLOAT, proc - 1, 0, MPI_COMM_WORLD, &status);
         MPI_Send(&u[iLeft][0], jmax + 2, MPI_FLOAT, proc - 1, 0, MPI_COMM_WORLD);
     }

     // MPI_Barrier(MPI_COMM_WORLD);

     MPI_Allgather(MPI_IN_PLACE, eachNode * (jmax + 2), MPI_FLOAT, &u[1][0], eachNode * (jmax + 2), MPI_FLOAT, MPI_COMM_WORLD);
     MPI_Allgather(MPI_IN_PLACE, eachNode * (jmax + 2), MPI_FLOAT, &v[1][0], eachNode * (jmax + 2), MPI_FLOAT, MPI_COMM_WORLD);
     MPI_Allgather(MPI_IN_PLACE, eachNode * (jmax + 2), MPI_FLOAT, &p[1][0], eachNode * (jmax + 2), MPI_FLOAT, MPI_COMM_WORLD);


   }

/* Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions (ie no particle moves more than one cell width in one
 * timestep). Otherwise the simulation becomes unstable.
 */
 void set_timestep_interval(float *del_t, int imax, int jmax, float delx,
     float dely, float **u, float **v, float Re, float tau)
 {
     int i, j;
     float umax, vmax, deltu, deltv, deltRe;

     float umaxSum, vmaxSum;

     /* del_t satisfying CFL conditions */
     if (tau >= 1.0e-10) { /* else no time stepsize control */
         umax = 1.0e-10;
         vmax = 1.0e-10;

         int iLeftTemp, iRightTemp;

         iLeftTemp = iLeft;
         if (proc == 0){
           iLeftTemp = 0;
         }


         iRightTemp = iRight;

         if (proc == nprocs -1){
           iRightTemp = imax + 1;
         }
         #pragma omp parallel for schedule(static) private(i,j) reduction(max:umax)
         for (i=iLeftTemp; i<=iRightTemp; i++) {
             for (j=1; j<=jmax+1; j++) {
                 umax = max(fabs(u[i][j]), umax);
             }
         }
           #pragma omp parallel for schedule(static) private(i,j) reduction(max:vmax)
           for (i=iLeftTemp; i<=iRightTemp; i++) {
             for (j=0; j<=jmax+1; j++) {
                 vmax = max(fabs(v[i][j]), vmax);
             }
         }

         MPI_Allreduce(&umax, &umaxSum, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
         MPI_Allreduce(&vmax, &vmaxSum, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

         deltu = delx/umaxSum;

         deltv = dely/vmaxSum;

         deltRe = 1/(1/(delx*delx)+1/(dely*dely))*Re/2.0;

         if (deltu<deltv) {
             *del_t = min(deltu, deltRe);
         } else {
             *del_t = min(deltv, deltRe);
         }
         *del_t = tau * (*del_t); /* multiply by safety factor */
     }
 }