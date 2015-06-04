#include "testings.h"

int main( int argc, char** argv )
{
    magma_opts opts;
    parse_opts( argc, argv, &opts );
    
    for( int i = 0; i < opts.ntest; ++i ) {
        printf( "m %5d, n %5d, k %5d\n", opts.msize[i], opts.nsize[i], opts.ksize[i] );
    }
    printf( "\n" );
    
    printf( "ntest    %d\n", opts.ntest );
    printf( "mmax     %d\n", opts.mmax  );
    printf( "nmax     %d\n", opts.nmax  );
    printf( "kmax     %d\n", opts.kmax  );
    printf( "\n" );
    
    printf( "nb       %d\n", opts.nb       ); 
    printf( "nrhs     %d\n", opts.nrhs     );
    printf( "nstream  %d\n", opts.nstream  );
    printf( "ngpu     %d\n", opts.ngpu     );
    printf( "niter    %d\n", opts.niter    );
    printf( "nthread  %d\n", opts.nthread  );
    printf( "itype    %d\n", opts.itype    );
    printf( "svd_work %d\n", opts.svd_work );
    printf( "\n" );
    
    printf( "check    %s\n", (opts.check  ? "true" : "false") );
    printf( "lapack   %s\n", (opts.lapack ? "true" : "false") );
    printf( "warmup   %s\n", (opts.warmup ? "true" : "false") );
    printf( "all      %s\n", (opts.all    ? "true" : "false") );
    printf( "\n" );
    
    printf( "uplo     %c\n", opts.uplo   );
    printf( "transA   %c\n", opts.transA );
    printf( "transB   %c\n", opts.transB );
    printf( "side     %c\n", opts.side   );
    printf( "diag     %c\n", opts.diag   );
    printf( "jobu     %c\n", opts.jobu   );
    printf( "jobvt    %c\n", opts.jobvt  );
    printf( "jobz     %c\n", opts.jobz   );
    printf( "jobvr    %c\n", opts.jobvr  );
    printf( "jobvl    %c\n", opts.jobvl  );
    
    return 0;
}
