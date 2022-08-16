//********************************************************************************************/
// this script generates array of debye-huckel potential in an intervals of r^2
//
//********************************************************************************************/
#include <stdio.h>
#include <math.h>

void dhenergytable_(float *screeningFactor,
                     float *saltCoefficient,
                     float *deConstant,
                     //float *esCutoffDistance,
                     float DebyeHuckelPotentials[],
                     float DebyeHuckelForces_over_r[])

{
    float DistanceSquared=0.01,Dist,interval=0.01,cutofflength,K=332.0;
    int i,intcutoff;
    
    //cutofflength=*esCutoffDistance/interval//delete
    //intcutoff=(int) cutofflength;    //DELETE this is the number of segments (the length of the array)
    for (i=0 ; i<8000000 ; i++)
    {
        Dist=sqrt(DistanceSquared);
        DebyeHuckelPotentials[i]= K*(*saltCoefficient)*exp(-(*screeningFactor)*(Dist))/((*deConstant)*(Dist));
        //printf ("Potentials %8.3f\n", DebyeHuckelPotentials[i]);
        DebyeHuckelForces_over_r[i]=DebyeHuckelPotentials[i]*(1/DistanceSquared + *screeningFactor/Dist);
        //printf ("DebyeHuckelForces_over_r %8.3f\n", DebyeHuckelForces_over_r[i]);
        DistanceSquared+= interval;
    }
    printf ("finished creating table 1\n");
    fflush(stdout);
}
