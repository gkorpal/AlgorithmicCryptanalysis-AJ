#include <stdio.h>
#include <stdlib.h>
#define DIM 32
#define WORD unsigned int
#define REPEAT 100000

#define bit(M,l,c) ((M[l]>>c)&1)
#define flipbit(M,l,c) if (1) {M[l]^=(1UL<<c);} else

void input_mat(WORD mat[DIM])
{
  int l,c,val;
  for (l=0;l<DIM;l++) {
    mat[l]=0;
    for (c=0;c<DIM;c++) {
      scanf("%d",&val);
      if (val) flipbit(mat,l,c);
    }
  }
}

void print_mat(WORD mat[DIM])
{
  int l,c;
  for (l=0;l<DIM;l++) {
    for (c=0;c<DIM;c++) {
      printf("%d ",bit(mat,l,c));
    }
    printf("\n");
  }
} 

void Mul(WORD res[DIM], WORD mat1[DIM], WORD mat2[DIM])
{
  int l,c,k;
  WORD val,mask;
  for (l=0;l<DIM;l++) res[l]=0;
  for (c=0;c<DIM;c++) { mask=0;
    for (k=0;k<DIM;k++) {
      mask^=(bit(mat2,k,c)<<k);
    }
    for (l=0;l<DIM;l++) {
      val=mat1[l]&mask;
      val^=(val>>16); val^=(val>>8);
      val^=(val>>4);  val^=(val>>2);
      val^=(val>>1);
      if (val&1) flipbit(res,l,c); 
    }
  }
}

main()
{
  WORD mat1[DIM]; WORD mat2[DIM]; WORD mat3[DIM]; int count;
  printf("Input Mat1\n");  input_mat(mat1);
  printf("Input Mat2\n");  input_mat(mat2);
  for(count=0; count<REPEAT; count++) Mul(mat3,mat1,mat2);
  printf("Product :\n");  print_mat(mat3);
}
