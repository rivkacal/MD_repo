//********************************************************************************************/
// electrostatics energy term
//
//********************************************************************************************/
#include <stdio.h>
#include <math.h>
void debyehuckel_(int IC[],
		  int JC[],
		  float Q1Q2[],
		  float X[],
		  float Y[],
		  float Z[],
		  float Fx[],
		  float Fy[],
		  float Fz[] ,
		  int *NE,
		  float *Eelec,
		  float *DebyeHuckelPotentials,
		  float *DebyeHuckelForces) 
{
  //float K = 332.0; //?? not sure this is the correct number.
  //float epsilon = *deConstant;
  //float kappa = *screeningFactor;
  //float B_kappa = *saltCoefficient;
  //int i;
  int i,intDistSquared;
  float helper,F_over_r;

  for(i =0;i < *NE;i++)
  {
      int C1 = IC[i]-1;
      int C2 = JC[i]-1;
      float dx = X[C1] - X[C2];
      float dy = Y[C1] - Y[C2];
      float dz = Z[C1] - Z[C2];

      float r2 = pow(dx,2) + pow(dy,2) + pow(dz,2);
      if (r2<10)
      {
        printf ("*****************\n");
        printf ("R2: %8.3f; INDEX1: %d; INDEX2: %d\n", r2,C1,C2);
        //printf ("INDEX1: %d\n", C1); 
        //printf ("INDEX2: %d\n", C2);
      }
      if (r2 < 80000)
      {
	helper=100*r2;   //make sure to multiple at the inverse of interval in dhEnergyTable.c
	intDistSquared=(int) helper;
	*Eelec+= DebyeHuckelPotentials[intDistSquared]*Q1Q2[i];
        // force in the direction C1 to C2 ,devided by r (is used to compute forces in x,y,z directions)
        F_over_r = DebyeHuckelForces[intDistSquared]*Q1Q2[i];
      } else { 
        *Eelec = 0;
        F_over_r = 0;
      }

		
       Fx[C1]+=  F_over_r*dx;
       Fx[C2]-=  F_over_r*dx;
       Fy[C1]+=  F_over_r*dy;
       Fy[C2]-=  F_over_r*dy;
       Fz[C1]+=  F_over_r*dz;
       Fz[C2]-=  F_over_r*dz;
  }
}

void debyehuckelfactor_(float *sigma,
			float *deConstant,
			float *screeningFactor,
			float *saltCoefficient,
			float *esEnergy)
{
  float K = 332.0; //?? not sure this is the correct number.
  *esEnergy = K*(*saltCoefficient)*exp(-(*screeningFactor)*(sqrt(*sigma)))/((*deConstant)*(sqrt(*sigma)));
}

