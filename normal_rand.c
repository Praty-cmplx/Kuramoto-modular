#include<math.h>
double std_norm()
{
    double const pi = 3.1415926;
    double x1 = genrand64_real3();
    double x2 = genrand64_real3();
    double z = sqrt(-2*log(x1))*cos(2*pi*x2);
    return (z);
}

/*
void main()
{
    int seed = (unsigned) time(NULL);
	init_genrand64(seed);

    char fname[100];
    sprintf(fname, "normal_dist.dat");
    FILE *fp;
    fp = fopen(fname, "w");
    for (int i=0; i<200; i++)
        fprintf(fp, "%.5f \n", std_norm());
    fclose(fp);
}
*/
