#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>

#define WORD unsigned int

#define DIM 32
#define Cst16 0xffff
#define Cst8 0xff00ff
#define Cst4 0xf0f0f0f
#define Cst2 0x33333333
#define Cst1 0x55555555


#define CUTOFF VarCutOff
#define MODUP VarModUp

#define access(M,i,j,bsize) (&M[((i)*(bsize)*DIM)+((j)*DIM)])

int VarCutOff=1;
int VarModUp=4;

/* External procedure for 32x32 boolean matrix multiplication */
extern void Mul(WORD res[DIM], WORD mat1[DIM], WORD mat2[DIM]);

double runtime ()
{
  struct tms t;
  long clk_tck;
  double ret;

  clk_tck=sysconf(_SC_CLK_TCK);
  times(&t);
  ret=t.tms_utime;
  ret/=clk_tck;
  return(ret);
}


void matmul(WORD * A, WORD * B, WORD *Res, int bsize) {
  int i,j,k,l;
  WORD tmp[DIM];
  
  for(i=0;i<DIM*bsize*bsize;i++) 
    Res[i]=0;

  for(i=0;i<bsize;i++) 
    for(j=0;j<bsize;j++) {
      for(k=0;k<bsize;k++) {
	Mul(tmp, access(A,i,k,bsize), access(B,k,j,bsize));
	for(l=0;l<DIM;l++) {
	  access(Res,i,j,bsize)[l]^=tmp[l];
	}
      }
    }
}

int mult_bloc(WORD * A, WORD * B, WORD *Res, int block_size)
{
  int i,j,l;
  int hsize, bsize;
  WORD *tmp1, *tmp2, *tmp3;
  WORD val;


  if (block_size<=CUTOFF){
    matmul(A,B,Res,block_size);
  }
  else {
    for(i=0;i<block_size*block_size*DIM;i++) Res[i]=0;
    if (block_size==3) {
      hsize=1;
      bsize=2;
    }
    else if (block_size%(MODUP)==(MODUP-1)) {
      hsize=(block_size+1)/2;
      bsize=block_size;
    }
    else {
      hsize=block_size/2;
      bsize=2*hsize;
    }
    tmp1=(WORD*)malloc((hsize*hsize)*DIM*sizeof(*tmp1));
    tmp2=(WORD*)malloc((hsize*hsize)*DIM*sizeof(*tmp2));
    tmp3=(WORD*)malloc((hsize*hsize)*DIM*sizeof(*tmp3));
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  val=access(A,i,j,block_size)[l];
	  if ((i+hsize)<bsize) val^=access(A,i+hsize,j,block_size)[l];
	  access(tmp1,i,j,hsize)[l]=val;
	
	  val=access(B,i,j,block_size)[l];
	  if ((j+hsize)<bsize) val^=access(B,i,j+hsize,block_size)[l];
	  access(tmp2,i,j,hsize)[l]=val;
	}
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  val=access(Res,i,j,block_size)[l]^access(tmp3,i,j,hsize)[l];
	  access(Res,i,j,block_size)[l]=val;
	}
      }
    }
    
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  if ((j+hsize)<bsize) {
	    val=access(A,i,j+hsize,block_size)[l];
	    if ((i+hsize)<bsize) val^=access(A,i+hsize,j+hsize,block_size)[l];
	  }
	  else val=0;
	  access(tmp1,i,j,hsize)[l]=val;
	
	  if ((i+hsize)<bsize) {
	    val=access(B,i+hsize,j,block_size)[l];
	    if ((j+hsize)<bsize) val^=access(B,i+hsize,j+hsize,block_size)[l];
	  }
	  else val=0;
	  access(tmp2,i,j,hsize)[l]=val;
	}
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<bsize-hsize;i++){
      for(j=0;j<bsize-hsize;j++){
	for(l=0;l<DIM;l++) {
	  val=access(Res,i+hsize,j+hsize,block_size)[l];
	  val^=access(tmp3,i,j,hsize)[l];
	  access(Res,i+hsize,j+hsize,block_size)[l]=val;
	}
      }
    }
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  if ((j+hsize)<bsize) 
	    val=access(A,i,j+hsize,block_size)[l];
	  else
	    val=0;
	  if ((i+hsize)<bsize)
	    val^=access(A,i+hsize,j,block_size)[l];
	  access(tmp1,i,j,hsize)[l]=val;
	  
	  if ((j+hsize)<bsize) 
	    val=access(B,i,j+hsize,block_size)[l];
	  else
	    val=0;
	  if ((i+hsize)<bsize) 
	    val^=access(B,i+hsize,j,block_size)[l];
	  access(tmp2,i,j,hsize)[l]=val;
	}
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  val=access(Res,i,j,block_size)[l];
	  val^=access(tmp3,i,j,hsize)[l];
	  access(Res,i,j,block_size)[l]=val;
	}
      }
    }

    for(i=0;i<bsize-hsize;i++){
      for(j=0;j<bsize-hsize;j++){
	for(l=0;l<DIM;l++) {
	  val=access(Res,i+hsize,j+hsize,block_size)[l];
	  val^=access(tmp3,i,j,hsize)[l];
	  access(Res,i+hsize,j+hsize,block_size)[l]=val;
	}
      }
    }
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  val=access(B,i,j,block_size)[l];
	  if ((i+hsize)<bsize) {
	    access(tmp1,i,j,hsize)[l]=access(A,i+hsize,j,block_size)[l];
	    val^=access(B,i+hsize,j,block_size)[l];
	  }
	  else
	    access(tmp1,i,j,hsize)[l]=0;
	  access(tmp2,i,j,hsize)[l]=val;
	}
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  val=access(Res,i,j,block_size)[l];
	  val^=access(tmp3,i,j,hsize)[l];
	  access(Res,i,j,block_size)[l]=val;
	  if ((i+hsize)<bsize) {
	    val=access(Res,i+hsize,j,block_size)[l];
	    val^=access(tmp3,i,j,hsize)[l];
	    access(Res,i+hsize,j,block_size)[l]=val;
	  }
	}
      }
    }
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  val=access(A,i,j,block_size)[l];
	  if ((j+hsize)<bsize)
	    val^=access(A,i,j+hsize,block_size)[l];
	  access(tmp1,i,j,hsize)[l]=val;
	  if ((j+hsize)<bsize)
	    access(tmp2,i,j,hsize)[l]=access(B,i,j+hsize,block_size)[l];
	  else
	    access(tmp2,i,j,hsize)[l]=0;
	}
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  val=access(Res,i,j,block_size)[l];
	  val^=access(tmp3,i,j,hsize)[l];
	  access(Res,i,j,block_size)[l]=val;
	  if ((j+hsize)<bsize){
	    val=access(Res,i,j+hsize,block_size)[l];
	    val^=access(tmp3,i,j,hsize)[l];
	    access(Res,i,j+hsize,block_size)[l]=val;
	  }
	}
      }
    }
    
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  if ((j+hsize)<bsize) {
	    access(tmp1,i,j,hsize)[l]=access(A,i,j+hsize,block_size)[l];
	    val=access(B,i,j+hsize,block_size)[l];
	    if ((i+hsize)<bsize) val^=access(B,i+hsize,j+hsize,block_size)[l];
	    access(tmp2,i,j,hsize)[l]=val;	
	  }
	  else {
	    access(tmp1,i,j,hsize)[l]=0;
	    access(tmp2,i,j,hsize)[l]=0;
	  }
	}
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<hsize;i++){
      for(j=0;j<bsize-hsize;j++){
	for(l=0;l<DIM;l++) {
	  if ((i+hsize)<bsize) {
	    val=access(Res,i+hsize,j+hsize,block_size)[l];
	    val^=access(tmp3,i,j,hsize)[l];
	    access(Res,i+hsize,j+hsize,block_size)[l]=val;
	  }
	  val=access(Res,i,j+hsize,block_size)[l];
	  val^=access(tmp3,i,j,hsize)[l];
	  access(Res,i,j+hsize,block_size)[l]=val;
	}
      }
    }
    
    for(i=0;i<hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  if ((i+hsize)<bsize) {
	    val=access(A,i+hsize,j,block_size)[l];
	    if ((j+hsize)<bsize)
	      val^=access(A,i+hsize,j+hsize,block_size)[l];
	    access(tmp1,i,j,hsize)[l]=val;
	    access(tmp2,i,j,hsize)[l]=access(B,i+hsize,j,block_size)[l];
	  }
	  else {
	    access(tmp1,i,j,hsize)[l]=0;
	    access(tmp2,i,j,hsize)[l]=0;
	  }
	}
      }
    }
    mult_bloc(tmp1,tmp2,tmp3,hsize);
    for(i=0;i<bsize-hsize;i++){
      for(j=0;j<hsize;j++){
	for(l=0;l<DIM;l++) {
	  if ((j+hsize)<bsize) {
	    val=access(Res,i+hsize,j+hsize,block_size)[l];
	    val^=access(tmp3,i,j,hsize)[l];
	    access(Res,i+hsize,j+hsize,block_size)[l]=val;
	  }
	  val=access(Res,i+hsize,j,block_size)[l];
	  val^=access(tmp3,i,j,hsize)[l];
	  access(Res,i+hsize,j,block_size)[l]=val;
	}
      }
    }
    
    free(tmp1); free(tmp2); free(tmp3);

    if (bsize==(block_size-1)) {
      for (i=0;i<bsize;i++) {
	for (j=0;j<bsize;j++) {
	  WORD tmpmat[DIM];
	  
	  Mul(tmpmat, access(A,i,block_size-1,block_size),
	      access(B,block_size-1,j,block_size));
	  for(l=0;l<DIM;l++)
	    access(Res,i,j,block_size)[l]^=tmpmat[l];
	}
      }

      for (i=0;i<block_size-1;i++) {
	for (j=0;j<block_size;j++) {
	  WORD tmpmat[DIM];
	  
	  Mul(tmpmat, access(A,i,j,block_size),
	      access(B,j,block_size-1,block_size));
	  for(l=0;l<DIM;l++)
	    access(Res,i,block_size-1,block_size)[l]^=tmpmat[l];
	  Mul(tmpmat, access(A,block_size-1,j,block_size),
	      access(B,j,i,block_size));
	  for(l=0;l<DIM;l++)
	    access(Res,block_size-1,i,block_size)[l]^=tmpmat[l];
	}
      }
      for (i=0;i<block_size;i++) {
	WORD tmpmat[DIM];
	
	Mul(tmpmat, access(A,block_size-1,i,block_size),
	    access(B,i,block_size-1,block_size));
	for(l=0;l<DIM;l++)
	  access(Res,block_size-1,block_size-1,block_size)[l]^=tmpmat[l];
      }
    }
  }
}

#define MAXSIZE 512
main()
{
  int i,j,k,sum;
  double timer;

  int BSIZE;
  WORD *m1,*m2,*m3,*m4;

  for(BSIZE=2;BSIZE<MAXSIZE;BSIZE++) {
    m1=malloc(BSIZE*BSIZE*DIM*sizeof(*m1));
    m2=malloc(BSIZE*BSIZE*DIM*sizeof(*m2));
    m3=malloc(BSIZE*BSIZE*DIM*sizeof(*m3));
    m4=malloc(BSIZE*BSIZE*DIM*sizeof(*m4));
    
    for(i=0;i<BSIZE;i++)
      for(j=0;j<BSIZE;j++) {
	for(k=0;k<DIM;k++) {
	  access(m1,i,j,BSIZE)[k]=random()^(random()<<16);
	  access(m2,i,j,BSIZE)[k]=random()^(random()<<16);
	}
      }

    timer=runtime();
    matmul(m1,m2,m4,BSIZE);
    timer=runtime()-timer;
    printf("%d %0.2f 0\n",BSIZE,timer);
    fflush(stdout);
    
    for (VarCutOff=1;VarCutOff<=8;VarCutOff++) {
      for (VarModUp=4;VarModUp<=64;VarModUp<<=1) {
	timer=runtime();
	mult_bloc(m1,m2,m3,BSIZE);
	timer=runtime()-timer;
	printf("%d %0.2f %d\n",BSIZE,timer,VarCutOff+2*VarModUp);
	fflush(stdout);

#if 1
	/* test correctness */
	for(i=0;i<BSIZE;i++)
	  for(j=0;j<BSIZE;j++) {
	    for(k=0;k<DIM;k++) {
	      if (access(m3,i,j,BSIZE)[k]!=access(m4,i,j,BSIZE)[k]) {
		printf("ARGGGG (%d, %d, %d)!!\n",i,j,k);
	      }
	    }
	  }
#endif

      }
    }
    
    free(m1);
    free(m2);
    free(m3);
    free(m4);
  }
}
