/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MAGMA_zpotri_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( N ) + FMULS_POTRI( N ))
#define _FADDS (FADDS_POTRF( N ) + FADDS_POTRI( N ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    MAGMA_enum uplo = PlasmaLower;

    LDA = max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );

    /* 
     * Initialize Data 
     * It's done in static to avoid having the same sequence than one 
     * the function we want to trace
     */
    MAGMA_zplghe_Tile( (double)N, descA, 51 );
    
    /* MAGMA ZPOTRF / ZTRTRI / ZLAUUM  */
    /*
     * Example of the different way to combine several asynchonous calls
     */
#if defined(TRACE_BY_SEQUENCE)
    {
        magma_sequence_t *sequence[3];
        magma_request_t request[3] = { PLASMA_REQUEST_INITIALIZER, 
                                      MAGMA_REQUEST_INITIALIZER, 
                                      MAGMA_REQUEST_INITIALIZER };
        
        MAGMA_Sequence_Create(&sequence[0]);
        MAGMA_Sequence_Create(&sequence[1]);
        MAGMA_Sequence_Create(&sequence[2]);
        
        if ( ! iparam[IPARAM_ASYNC] ) {
            START_TIMING();

            MAGMA_zpotrf_Tile_Async(uplo, descA,                sequence[0], &request[0]);
            MAGMA_Sequence_Wait(sequence[0]);
            
            MAGMA_ztrtri_Tile_Async(uplo, PlasmaNonUnit, descA, sequence[1], &request[1]);
            MAGMA_Sequence_Wait(sequence[1]);
            
            MAGMA_zlauum_Tile_Async(uplo, descA,                sequence[2], &request[2]);
            MAGMA_Sequence_Wait(sequence[2]);
            STOP_TIMING();

        } else {

            START_TIMING();
            MAGMA_zpotrf_Tile_Async(uplo, descA,                sequence[0], &request[0]);
            MAGMA_ztrtri_Tile_Async(uplo, PlasmaNonUnit, descA, sequence[1], &request[1]);
            MAGMA_zlauum_Tile_Async(uplo, descA,                sequence[2], &request[2]);

            MAGMA_Sequence_Wait(sequence[0]);
            MAGMA_Sequence_Wait(sequence[1]);
            MAGMA_Sequence_Wait(sequence[2]);
            STOP_TIMING();
        }
        
        MAGMA_Sequence_Destroy(sequence[0]);
        MAGMA_Sequence_Destroy(sequence[1]);
        MAGMA_Sequence_Destroy(sequence[2]);
    }       
#else
    {
        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            MAGMA_zpotrf_Tile(uplo, descA);
            MAGMA_ztrtri_Tile(uplo, PlasmaNonUnit, descA);
            MAGMA_zlauum_Tile(uplo, descA);
            STOP_TIMING();

        } else {

            /* Default: we use Asynchonous call with only one sequence */
            magma_sequence_t *sequence;
            magma_request_t request[2] = { PLASMA_REQUEST_INITIALIZER, 
                                          MAGMA_REQUEST_INITIALIZER };
        
            START_TIMING();
            MAGMA_Sequence_Create(&sequence);
            MAGMA_zpotrf_Tile_Async(uplo, descA, sequence, &request[0]);
            MAGMA_zpotri_Tile_Async(uplo, descA, sequence, &request[1]);
            MAGMA_Sequence_Wait(sequence);
            STOP_TIMING();
        
            MAGMA_Sequence_Destroy(sequence);       
        }
    }
#endif
    
    /* Check the solution */
    if ( check )
    {
        dparam[IPARAM_ANORM] = 0.0;
        dparam[IPARAM_XNORM] = 0.0;
        dparam[IPARAM_BNORM] = 0.0;
        dparam[IPARAM_RES]   = 0.0;
    }

    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
