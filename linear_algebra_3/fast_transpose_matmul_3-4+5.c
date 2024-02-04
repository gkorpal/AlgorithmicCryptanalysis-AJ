#include <stdio.h>
#include <stdlib.h>

#define DIM 32
#define WORD unsigned int
#define REPEAT 100000

#define bit(M,l,c) ((M[l]>>c)&1)
#define xorbit(M,l,c) if (1) {M[l]^=(1UL<<c);} else

#ifndef NO_MAIN
void input_mat(WORD mat[DIM])
{
  int l,c,val;
  
  for (l=0;l<DIM;l++) {
    mat[l]=0;
    for (c=0;c<DIM;c++) {
      scanf("%d",&val);
      if (val) xorbit(mat,l,c);
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
#endif

void Transpose(WORD transp[DIM], WORD mat[DIM])
{
  int l0,l,c;
  WORD val1,val2;

  for (l=0;l<DIM/2;l++) {
    transp[l]=(mat[l]&0xffff)|((mat[l+(DIM/2)]&0xffff)<<16);
    transp[l+(DIM/2)]=((mat[l]&0xffff0000)>>16)|(mat[l+(DIM/2)]&0xffff0000);
  }
  for (l0=0;l0<DIM;l0+=DIM/2) {
    for (l=l0;l<l0+(DIM/4);l++) {
      val1=(transp[l]&0xff00ff)|((transp[l+(DIM/4)]&0xff00ff)<<8);
      val2=((transp[l]&0xff00ff00)>>8)|(transp[l+(DIM/4)]&0xff00ff00);
      transp[l]=val1;
      transp[l+(DIM/4)]=val2;
    }
  }
  for (l0=0;l0<DIM;l0+=DIM/4) {
    for (l=l0;l<l0+(DIM/8);l++) {
      val1=(transp[l]&0xf0f0f0f)|((transp[l+(DIM/8)]&0xf0f0f0f)<<4);
      val2=((transp[l]&0xf0f0f0f0)>>4)|(transp[l+(DIM/8)]&0xf0f0f0f0);
      transp[l]=val1;
      transp[l+(DIM/8)]=val2;
    }
  }
  for (l0=0;l0<DIM;l0+=DIM/8) {
    for (l=l0;l<l0+(DIM/16);l++) {
      val1=(transp[l]&0x33333333)|((transp[l+(DIM/16)]&0x33333333)<<2);
      val2=((transp[l]&0xcccccccc)>>2)|(transp[l+(DIM/16)]&0xcccccccc);
      transp[l]=val1;
      transp[l+(DIM/16)]=val2;
    }
  }

  for (l=0;l<DIM;l+=2) {
    val1=(transp[l]&0x55555555)|((transp[l+1]&0x55555555)<<1);
    val2=((transp[l]&0xaaaaaaaa)>>1)|(transp[l+1]&0xaaaaaaaa);
    transp[l]=val1;
    transp[l+1]=val2;
  }
}

void Mul(WORD res[DIM], WORD mat1[DIM], WORD mat2[DIM])
{
  int l,c,k;
  WORD transp[DIM];
  WORD tmp[DIM];
  WORD val;

  Transpose(transp,mat2);

  for (l=0;l<DIM;l++) {
    for (c=0;c<DIM;c++) {
      val=mat1[l]&transp[c];
      val^=(val>>16); val&=0xffff;
      tmp[c]=val;
    }
    for (c=0;c<DIM/2;c++) {
      val=tmp[c]|(tmp[c+(DIM/2)]<<(DIM/2));
      tmp[c]=(val&0xff00ff)^((val>>(DIM/4))&0xff00ff);
    }
    for (c=0;c<DIM/4;c++) {
      val=tmp[c]|(tmp[c+(DIM/4)]<<(DIM/4));
      tmp[c]=(val&0xf0f0f0f)^((val>>(DIM/8))&0xf0f0f0f);
    }
    for (c=0;c<DIM/8;c++) {
      val=tmp[c]|(tmp[c+(DIM/8)]<<(DIM/8));
      tmp[c]=(val&0x33333333)^((val>>(DIM/16))&0x33333333);
    }
    val=tmp[0]|(tmp[2]<<2);
    tmp[0]=(val&0x55555555)^((val>>1)&0x55555555);
    val=tmp[1]|(tmp[3]<<2);
    tmp[1]=(val&0x55555555)^((val>>1)&0x55555555);
    val=tmp[0]|(tmp[1]<<1);
    res[l]=val;
  }
}

#ifndef NO_MAIN
main()
{
  WORD mat1[DIM]; WORD mat2[DIM]; WORD mat3[DIM]; int count;

  printf("Input Mat1\n");  input_mat(mat1);
  printf("Input Mat2\n");  input_mat(mat2);
  for(count=0; count<REPEAT; count++) Mul(mat3,mat1,mat2);
  printf("Product :\n");  print_mat(mat3);
}
#endif
