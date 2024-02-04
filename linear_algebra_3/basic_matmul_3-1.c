/*
This work is licensed under the Creative Commons
Attribution-Noncommercial-Share Alike 3.0 Unported License. To view a
copy of this license, visit
http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to
Creative Commons, 171 Second Street, Suite 300, San Francisco,
California, 94105, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#define DIM 32
#define REPEAT 100000

void input_mat(int mat[DIM][DIM])
{
  int l,c;
  for (l=0;l<DIM;l++) {
    for (c=0;c<DIM;c++) {
      scanf("%d",&mat[l][c]);
    }
  }
}

void print_mat(int mat[DIM][DIM])
{
  int l,c;
  for (l=0;l<DIM;l++) {
    for (c=0;c<DIM;c++) {
      printf("%d ",mat[l][c]);
    }
    printf("\n");
  }
}

void Mul(int res[DIM][DIM], int mat1[DIM][DIM], int mat2[DIM][DIM])
{
  int l,c,k;
  for (l=0;l<DIM;l++) {
    for (c=0;c<DIM;c++) { res[l][c]=0;
      for (k=0;k<DIM;k++) {
        res[l][c]+=mat1[l][k]*mat2[k][c];
      }
      res[l][c]%=2;
    }
  }
}

main()
{
  int mat1[DIM][DIM]; int mat2[DIM][DIM];
  int mat3[DIM][DIM]; int count;
  printf("Input Mat1\n");  input_mat(mat1);
  printf("Input Mat2\n");  input_mat(mat2);
  for (count=0;count<REPEAT;count++) Mul(mat3,mat1,mat2);
  printf("Product :\n");  print_mat(mat3);
}
