#include "smol.h"
 
int main()
{
          
  SMOL <float> matOp;
  float **mat1 = matOp.criarMatriz(3,2),
      **mat2 = matOp.criarMatriz(2,3),
      **mat3 = matOp.criarMatriz(3,3),
      **resultado1,
      **resultado2,
      **resultado3;
  float **resultado2Inversa = matOp.criarMatriz(3,3);
  
  int valor = 0;
  int i,j;
  
  for (i=0; i<3; i++){
    mat3[i][i] = 1;
    for (j=0; j<2; j++){      
      mat1[i][j] = ++valor;
      mat2[j][i] = valor;      
    }
  }
  mat3[0][0] = 2;
  mat3[0][1] = -4;
  mat3[0][2] = 4;
  mat3[1][0] = 5;
  mat3[1][1] = 4;
  mat3[1][2] = 6;
  mat3[2][0] = -3;
  mat3[2][1] = 0;
  mat3[2][2] = 2;
  
  printf("mat1 \n");  
  matOp.printMatriz(mat1,3,2,"\t%f");
  printf("mat2 \n");  
  matOp.printMatriz(mat2,2,3,"\t%f");
  printf("mat3 \n");  
  matOp.printMatriz(mat3,3,3,"\t%f");
  
  printf("mat1 * mat2\n");  
  //multiplicacao
  resultado1 = matOp.matrizMult(mat1,3,2,mat2,2,3);
  matOp.printMatriz(resultado1,3,3,"\t%f");
  //printf("\n");
  //multiplicacao
  //resultado2 = matOp.matrizMult(resultado1,3,3,mat3,3,3);
  //matOp.printMatriz(resultado2,3,3,"|%d|");  
  
  printf("\n");  
  printf("The mat3's determinant is %f \n",matOp.calculaDeterminante(mat3,3));
  matOp.inversaoMatriz(mat3,3,resultado2Inversa);
  matOp.printMatriz(resultado2Inversa,3,3,"\t%f");
  
  printf("\n");
  //transposta
  //resultado3 = matOp.matrizTrans(mat1,3,2); 
  //matOp.printMatriz(resultado3,2,3,"|%d|");
     
  //desalocando a matriz  
  matOp.desalocaMatriz(mat1,2);
  matOp.desalocaMatriz(mat2,3);
  matOp.desalocaMatriz(mat3,3);
  //matOp.desalocaMatriz(resultado1,3);
  //matOp.desalocaMatriz(resultado2,3);
  //matOp.desalocaMatriz(resultado3,3);
  matOp.desalocaMatriz(resultado2Inversa,3);
  
  
  return 0;
}
