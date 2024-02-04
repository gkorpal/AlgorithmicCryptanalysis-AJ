#include <stdio.h>
#include <stdlib.h>
#define DIM 32
#define WORD unsigned int
#define bit(M,l,c) ((M[l]>>c)&1)
#define flipbit(M,l,c) if (1) {M[l]^=(1UL<<c);} else
#define FALSE 0
#define TRUE 1

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

void Mul(WORD res[DIM],WORD mat1[DIM],WORD mat2[DIM])
{
  int l,c,k,val;
  for (l=0;l<DIM;l++) {
    res[l]=0;
    for (c=0;c<DIM;c++) {
      val=0;
      for (k=0;k<DIM;k++) {
        val^=bit(mat1,l,k)&bit(mat2,k,c);
      }
      if (val) flipbit(res,l,c);
    }
  }
}

/* Warning: mat is transformed during MatInv */
int MatInv(WORD mat[DIM], WORD inv[DIM])
{
  int piv,l,c,k;
  WORD val,vali,mask;
  for(piv=0,mask=1;piv<DIM;piv++,mask<<=1)
    inv[piv]=mask;
  for(piv=0,mask=1;piv<DIM;piv++,mask<<=1) {
    for (c=piv;c<DIM;c++) if (mask&mat[c]) break;
    if (c>=DIM) return(FALSE);
    val=mat[c];mat[c]=mat[piv];mat[piv]=val;
    vali=inv[c];inv[c]=inv[piv];inv[piv]=vali;
    for(c=0;c<DIM;c++) if ((c!=piv)&&(mask&mat[c])) {
      mat[c]^=val;inv[c]^=vali; }}
  return(TRUE);
}

main()
{
  WORD mat1[DIM]; WORD mat2[DIM]; WORD matinv[DIM];
  int i;

  printf("Input Matrix\n");  input_mat(mat1);
  for(i=0;i<DIM;i++) mat2[i]=mat1[i];
  if (MatInv(mat2,matinv)==FALSE)
    printf("Non invertible matrix\n");
  else {
    printf("Inverse :\n");
    print_mat(matinv);
    Mul(mat2,mat1,matinv);
    printf("Correctness Test (should give Identity matrix) :\n");
    print_mat(mat2);
  }
}
