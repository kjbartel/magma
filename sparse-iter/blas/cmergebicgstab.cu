/*
    -- MAGMA (version 1.5.0-beta3) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date July 2014

       @generated from zmergebicgstab.cu normal z -> c, Fri Jul 18 17:34:28 2014
       @author Hartwig Anzt

*/
#include "common_magma.h"

#define BLOCK_SIZE 512

#define PRECISION_c


// These routines merge multiple kernels from cmergebicgstab into one
// The difference to cmergedbicgstab2 is that the SpMV is not merged into the
// kernes. This results in higher flexibility at the price of lower performance.

/* -------------------------------------------------------------------------- */

__global__ void 
magma_cbicgmerge1_kernel(  
                    int n, 
                    magmaFloatComplex *skp,
                    magmaFloatComplex *v, 
                    magmaFloatComplex *r, 
                    magmaFloatComplex *p ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    magmaFloatComplex beta=skp[1];
    magmaFloatComplex omega=skp[2];
    if( i<n ){
        p[i] =  r[i] + beta * ( p[i] - omega * v[i] );

    }

}

/**
    Purpose
    -------

    Mergels multiple operations into one kernel:

    p = beta*p
    p = p-omega*beta*v
    p = p+r
    
    -> p = r + beta * ( p - omega * v ) 

    Arguments
    ---------

    @param
    n           int
                dimension n

    @param
    skp         magmaFloatComplex*
                set of scalar parameters

    @param
    v           magmaFloatComplex*
                input v

    @param
    r           magmaFloatComplex*
                input r

    @param
    p           magmaFloatComplex*
                input/output p


    @ingroup magmasparse_cgegpuk
    ********************************************************************/

extern "C" int
magma_cbicgmerge1(  int n, 
                    magmaFloatComplex *skp,
                    magmaFloatComplex *v, 
                    magmaFloatComplex *r, 
                    magmaFloatComplex *p ){

    
    dim3 Bs( BLOCK_SIZE );
    dim3 Gs( (n+BLOCK_SIZE-1)/BLOCK_SIZE );
    magma_cbicgmerge1_kernel<<<Gs, Bs, 0>>>( n, skp, v, r, p );

   return MAGMA_SUCCESS;
}

/* -------------------------------------------------------------------------- */

__global__ void 
magma_cbicgmerge2_kernel(  
                    int n, 
                    magmaFloatComplex *skp, 
                    magmaFloatComplex *r,
                    magmaFloatComplex *v, 
                    magmaFloatComplex *s ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    magmaFloatComplex alpha=skp[0];
    if( i<n ){
        s[i] =  r[i] - alpha * v[i] ;
    }

}

/**
    Purpose
    -------

    Mergels multiple operations into one kernel:

    s=r
    s=s-alpha*v
        
    -> s = r - alpha * v

    Arguments
    ---------

    @param
    n           int
                dimension n

    @param
    skp         magmaFloatComplex*
                set of scalar parameters

    @param
    r           magmaFloatComplex*
                input r

    @param
    v           magmaFloatComplex*
                input v

    @param
    s           magmaFloatComplex*
                input/output s


    @ingroup magmasparse_cgegpuk
    ********************************************************************/

extern "C" int
magma_cbicgmerge2(  int n, 
                    magmaFloatComplex *skp, 
                    magmaFloatComplex *r,
                    magmaFloatComplex *v, 
                    magmaFloatComplex *s ){

    
    dim3 Bs( BLOCK_SIZE );
    dim3 Gs( (n+BLOCK_SIZE-1)/BLOCK_SIZE );

    magma_cbicgmerge2_kernel<<<Gs, Bs, 0>>>( n, skp, r, v, s );

   return MAGMA_SUCCESS;
}

/* -------------------------------------------------------------------------- */

__global__ void 
magma_cbicgmerge3_kernel(  
                    int n, 
                    magmaFloatComplex *skp, 
                    magmaFloatComplex *p,
                    magmaFloatComplex *se,
                    magmaFloatComplex *t,
                    magmaFloatComplex *x, 
                    magmaFloatComplex *r ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    magmaFloatComplex alpha=skp[0];
    magmaFloatComplex omega=skp[2];
    if( i<n ){
        magmaFloatComplex s;
        s = se[i];
        x[i] = x[i] + alpha * p[i] + omega * s;
        r[i] = s - omega * t[i];
    }

}

/**
    Purpose
    -------

    Mergels multiple operations into one kernel:

    x=x+alpha*p
    x=x+omega*s
    r=s
    r=r-omega*t
        
    -> x = x + alpha * p + omega * s
    -> r = s - omega * t

    Arguments
    ---------

    @param
    n           int
                dimension n

    @param
    skp         magmaFloatComplex*
                set of scalar parameters

    @param
    p           magmaFloatComplex*
                input p

    @param
    s           magmaFloatComplex*
                input s

    @param
    t           magmaFloatComplex*
                input t

    @param
    x           magmaFloatComplex*
                input/output x

    @param
    r           magmaFloatComplex*
                input/output r


    @ingroup magmasparse_cgegpuk
    ********************************************************************/

extern "C" int
magma_cbicgmerge3(  int n, 
                    magmaFloatComplex *skp,
                    magmaFloatComplex *p,
                    magmaFloatComplex *s,
                    magmaFloatComplex *t,
                    magmaFloatComplex *x, 
                    magmaFloatComplex *r ){

    
    dim3 Bs( BLOCK_SIZE );
    dim3 Gs( (n+BLOCK_SIZE-1)/BLOCK_SIZE );
    magma_cbicgmerge3_kernel<<<Gs, Bs, 0>>>( n, skp, p, s, t, x, r );

   return MAGMA_SUCCESS;
}

/* -------------------------------------------------------------------------- */

__global__ void 
magma_cbicgmerge4_kernel_1(  
                    magmaFloatComplex *skp ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        magmaFloatComplex tmp = skp[0];
        skp[0] = skp[4]/tmp;
    }
}

__global__ void 
magma_cbicgmerge4_kernel_2(  
                    magmaFloatComplex *skp ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        skp[2] = skp[6]/skp[7];
        skp[3] = skp[4];
    }
}

__global__ void 
magma_cbicgmerge4_kernel_3(  
                    magmaFloatComplex *skp ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        magmaFloatComplex tmp1 = skp[4]/skp[3];
        magmaFloatComplex tmp2 = skp[0] / skp[2];
        skp[1] =  tmp1*tmp2;
        //skp[1] =  skp[4]/skp[3] * skp[0] / skp[2];

    }
}

/**
    Purpose
    -------

    Performs some parameter operations for the BiCGSTAB with scalars on GPU.

    Arguments
    ---------

    @param
    type        int
                kernel type

    @param
    skp         magmaFloatComplex*
                vector with parameters


    @ingroup magmasparse_cgegpuk
    ********************************************************************/

extern "C" int
magma_cbicgmerge4(  int type, 
                    magmaFloatComplex *skp ){

    dim3 Bs( 2 );
    dim3 Gs( 1 );
    if( type == 1 )
        magma_cbicgmerge4_kernel_1<<<Gs, Bs, 0>>>( skp );
    else if( type == 2 )
        magma_cbicgmerge4_kernel_2<<<Gs, Bs, 0>>>( skp );
    else if( type == 3 )
        magma_cbicgmerge4_kernel_3<<<Gs, Bs, 0>>>( skp );
    else
        printf("error: no kernel called\n");

   return MAGMA_SUCCESS;
}

