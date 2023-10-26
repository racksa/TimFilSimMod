// Online C++ compiler to run C++ program online
#include <iostream>
#include <cmath>
// typedef double Real;
typedef float Real;

double func1(const double a){
    return a*2.0;
}

int main() {
    // Write C++ code here
    Real scalar_part = -0.9999999689378366;
    scalar_part = -1.0;
    printf("test=%.16f\n", sqrt(0.5*(1.0 + scalar_part)));
    // double scalar_part = scalar_part_r;
    
    Real vector_part[3] = {0.0001603911234386, -0.0001193993535068, 0.0001488045963561};
    // double vector_part[3] ={vector_part_r[0], vector_part_r[1], vector_part_r[2] };
    
    printf("norm = %.16f\n",sqrt(scalar_part*scalar_part+ vector_part[0]*vector_part[0]+ vector_part[1]*vector_part[1]+ vector_part[2]*vector_part[2]) );
    printf("(%.16f %.16f %.16f %.16f)\n", scalar_part, vector_part[0], vector_part[1], vector_part[2]);
    
    if (scalar_part < -0.9999999) {
      const Real dir_norm = sqrt(vector_part[0]*vector_part[0] + vector_part[1]*vector_part[1] + vector_part[2]*vector_part[2]);
        printf("norm=%.16f\n", dir_norm*dir_norm);
      if(dir_norm*dir_norm < 0.00000001){
          printf("direct\n");
        scalar_part = 0.0;
        vector_part[0] /= dir_norm;
        vector_part[1] /= dir_norm;
        vector_part[2] /= dir_norm;
      }
      else{
         if(scalar_part==-1.0){
          scalar_part=-0.9999998;
        }
        printf("new\n");
        printf("1+scalar part=%.16f\n", 1.0 + scalar_part);
        scalar_part = sqrt(0.5*(1.0 + scalar_part));
        printf("scalar=%.16f\n", scalar_part);
        const Real temp = 2.0*scalar_part;
        vector_part[0] /= temp;
        vector_part[1] /= temp;
        vector_part[2] /= temp;
        
        Real q_norm = sqrt(scalar_part*scalar_part +vector_part[0]*vector_part[0] + vector_part[1]*vector_part[1] + vector_part[2]*vector_part[2]);
        scalar_part *= q_norm;
        vector_part[0] /= q_norm;
        vector_part[1] /= q_norm;
        vector_part[2] /= q_norm;
        
        printf("qnorm=%.16f\n", q_norm);
      }

    } else {
        printf("compute\n");
      scalar_part = sqrt(0.5*(1.0 + scalar_part));
      const Real temp = 2.0*scalar_part;
      vector_part[0] /= temp;
      vector_part[1] /= temp;
      vector_part[2] /= temp;

    }
    
    float b = func1(scalar_part);
    
     printf("(%.16f %.16f %.16f %.16f)\n", scalar_part, vector_part[0], vector_part[1], vector_part[2]);
    printf("norm = %.16f\n",sqrt(scalar_part*scalar_part+ vector_part[0]*vector_part[0]+ vector_part[1]*vector_part[1]+ vector_part[2]*vector_part[2]) );
    printf("\nnorm = 1.0000000000000002\n(-0.9999999689378366 0.0001603911234386 -0.0001193993535068 0.0001488045963561)\n norm=0.0000000621243260\nnew\n1+scalar part=0.0000000310621634\n(0.0001246237605695 0.6435013784915440 -0.4790392817594920 0.5970153511499336)\nnorm = 1.0000000026289275");
    return 0;
}