/*
    -- MAGMA (version 1.5.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date September 2014

       @generated from magmasparse_z.h normal z -> s, Wed Sep 17 15:08:42 2014
       @author Hartwig Anzt
*/

#ifndef MAGMASPARSE_S_H
#define MAGMASPARSE_S_H

#include "magma_types.h"
#include "magmasparse_types.h"

#define PRECISION_s


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_SPARSE Matrix Descriptors
*/
/* CSR Matrix descriptor */
typedef struct {
    int type;

    magma_int_t   m;
    magma_int_t   n;
    magma_int_t nnz;

    float *d_val;
    magma_int_t *d_rowptr;
    magma_int_t *d_colind;

} magma_smatrix_t;


/* BCSR Matrix descriptor */
typedef struct {
    int type;

    magma_int_t   rows_block;
    magma_int_t   cols_block;

    magma_int_t nrow_blocks;
    magma_int_t  nnz_blocks;

    float *d_val;
    magma_int_t *d_rowptr;
    magma_int_t *d_colind;

} magma_sbcsr_t;


#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_SPARSE Auxiliary functions
*/


magma_int_t
magma_sparse_opts(      int argc, 
                        char** argv, 
                        magma_sopts *opts, 
                        int *matrices );

magma_int_t 
read_s_csr_from_binary( magma_int_t* n_row, 
                        magma_int_t* n_col, 
                        magma_int_t* nnz, 
                        float **val, 
                        magma_index_t **row, 
                        magma_index_t **col,
                        const char * filename);

magma_int_t 
read_s_csr_from_mtx(    magma_storage_t *type, 
                        magma_location_t *location,
                        magma_int_t* n_row, 
                        magma_int_t* n_col, 
                        magma_int_t* nnz, 
                        float **val, 
                        magma_index_t **row, 
                        magma_index_t **col, 
                        const char *filename);

magma_int_t 
magma_s_csr_mtx(        magma_s_sparse_matrix *A, 
                        const char *filename );

magma_int_t 
magma_s_csr_mtxsymm(    magma_s_sparse_matrix *A, 
                        const char *filename );

magma_int_t 
magma_s_csr_compressor( float ** val, 
                        magma_index_t ** row, 
                        magma_index_t ** col, 
                        float ** valn, 
                        magma_index_t ** rown, 
                        magma_index_t ** coln, 
                        magma_int_t *n );

magma_int_t 
magma_s_csrtranspose(   magma_s_sparse_matrix A, 
                        magma_s_sparse_matrix *B );

magma_int_t 
magma_s_cucsrtranspose( magma_s_sparse_matrix A, 
                        magma_s_sparse_matrix *B );

magma_int_t 
s_transpose_csr(        magma_int_t n_rows, 
                        magma_int_t n_cols, 
                        magma_int_t nnz,
                        float *val, 
                        magma_index_t *row, 
                        magma_index_t *col, 
                        magma_int_t *new_n_rows, 
                        magma_int_t *new_n_cols, 
                        magma_int_t *new_nnz, 
                        float **new_val, 
                        magma_index_t **new_row, 
                        magma_index_t **new_col );

magma_int_t
magma_scsrsplit(        magma_int_t bsize,
                        magma_s_sparse_matrix A,
                        magma_s_sparse_matrix *D,
                        magma_s_sparse_matrix *R );

magma_int_t
magma_smscale(          magma_s_sparse_matrix *A, 
                        magma_scale_t scaling );

magma_int_t 
magma_smdiff(           magma_s_sparse_matrix A, 
                        magma_s_sparse_matrix B, 
                        real_Double_t *res );

magma_int_t
magma_smdiagadd(        magma_s_sparse_matrix *A, 
                        float add );

magma_int_t 
magma_smsort(           magma_s_sparse_matrix *A );


magma_int_t
magma_ssymbilu(         magma_s_sparse_matrix *A, 
                        magma_int_t levels,
                        magma_s_sparse_matrix *L,
                        magma_s_sparse_matrix *U );


magma_int_t 
write_s_csr_mtx(        magma_int_t n_row, 
                        magma_int_t n_col, 
                        magma_int_t nnz, 
                        float **val, 
                        magma_index_t **row, 
                        magma_index_t **col, 
                        magma_order_t MajorType,
                        const char *filename );

magma_int_t 
write_s_csrtomtx(        magma_s_sparse_matrix A,
                        const char *filename );

magma_int_t 
print_s_csr(            magma_int_t n_row, 
                        magma_int_t n_col, 
                        magma_int_t nnz, 
                        float **val, 
                        magma_index_t **row, 
                        magma_index_t **col );

magma_int_t 
print_s_csr_mtx(        magma_int_t n_row, 
                        magma_int_t n_col, 
                        magma_int_t nnz, 
                        float **val, 
                        magma_index_t **row, 
                        magma_index_t **col, 
                        magma_order_t MajorType );


magma_int_t 
magma_s_mtranspose(     magma_s_sparse_matrix A, 
                        magma_s_sparse_matrix *B );


magma_int_t 
magma_s_mtransfer(      magma_s_sparse_matrix A, 
                        magma_s_sparse_matrix *B, 
                        magma_location_t src, 
                        magma_location_t dst );

magma_int_t 
magma_s_vtransfer(      magma_s_vector x, 
                        magma_s_vector *y, 
                        magma_location_t src, 
                        magma_location_t dst );

magma_int_t 
magma_s_mconvert(       magma_s_sparse_matrix A, 
                        magma_s_sparse_matrix *B, 
                        magma_storage_t old_format, 
                        magma_storage_t new_format );


magma_int_t
magma_s_vinit(          magma_s_vector *x, 
                        magma_location_t memory_location,
                        magma_int_t num_rows, 
                        float values );

magma_int_t
magma_s_vvisu(          magma_s_vector x, 
                        magma_int_t offset, 
                        magma_int_t displaylength );

magma_int_t
magma_s_vread(          magma_s_vector *x, 
                        magma_int_t length,
                        char * filename );
magma_int_t
magma_s_mvisu(          magma_s_sparse_matrix A );

magma_int_t 
magma_sdiameter(        magma_s_sparse_matrix *A );

magma_int_t 
magma_srowentries(      magma_s_sparse_matrix *A );

magma_int_t
magma_s_mfree(          magma_s_sparse_matrix *A );

magma_int_t
magma_s_vfree(          magma_s_vector *x );

magma_int_t
magma_sresidual(        magma_s_sparse_matrix A, 
                        magma_s_vector b, 
                        magma_s_vector x, 
                        float *res );
magma_int_t
magma_smgenerator(  magma_int_t n,
                    magma_int_t offdiags,
                    magma_index_t *diag_offset,
                    float *diag_vals,
                    magma_s_sparse_matrix *A );

magma_int_t
magma_sm_27stencil(  magma_int_t n,
                     magma_s_sparse_matrix *A );

magma_int_t
magma_sm_5stencil(  magma_int_t n,
                     magma_s_sparse_matrix *A );

magma_int_t
magma_ssolverinfo(  magma_s_solver_par *solver_par, 
                    magma_s_preconditioner *precond_par );

magma_int_t
magma_ssolverinfo_init( magma_s_solver_par *solver_par, 
                        magma_s_preconditioner *precond );

magma_int_t
magma_ssolverinfo_free( magma_s_solver_par *solver_par, 
                        magma_s_preconditioner *precond );


/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_SPARSE function definitions / Data on CPU
*/




/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_SPARSE function definitions / Data on CPU / Multi-GPU
*/

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_SPARSE solvers (Data on GPU)
*/

magma_int_t 
magma_scg(             magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par );

magma_int_t 
magma_scg_res(         magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par );

magma_int_t 
magma_scg_merge(       magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par );

magma_int_t 
magma_sgmres(          magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par );

magma_int_t
magma_sbicgstab(       magma_s_sparse_matrix A, magma_s_vector b, magma_s_vector *x,  
                       magma_s_solver_par *solver_par );

magma_int_t
magma_sbicgstab_merge( magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par );

magma_int_t
magma_sbicgstab_merge2( magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par );

magma_int_t
magma_spcg(            magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par, 
                       magma_s_preconditioner *precond_par );

magma_int_t
magma_spbicgstab(      magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par, 
                       magma_s_preconditioner *precond_par );

magma_int_t
magma_spgmres(         magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par, 
                       magma_s_preconditioner *precond_par );
magma_int_t
magma_sjacobi(         magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par );

magma_int_t
magma_sbaiter(         magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par );

magma_int_t
magma_siterref(        magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par, 
                       magma_s_preconditioner *precond_par );

magma_int_t
magma_silu(            magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par, 
                       magma_int_t *ipiv );

magma_int_t
magma_sbcsrlu(         magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par );


magma_int_t
magma_sbcsrlutrf(      magma_s_sparse_matrix A, 
                       magma_s_sparse_matrix *M,
                       magma_int_t *ipiv, 
                       magma_int_t version );

magma_int_t
magma_sbcsrlusv(       magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par, 
                       magma_int_t *ipiv );



magma_int_t
magma_silucg(          magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par );

magma_int_t
magma_silugmres(       magma_s_sparse_matrix A, magma_s_vector b, 
                       magma_s_vector *x, magma_s_solver_par *solver_par ); 


magma_int_t
magma_slobpcg_shift(    magma_int_t num_rows,
                        magma_int_t num_vecs, 
                        magma_int_t shift,
                        float *x );
magma_int_t
magma_slobpcg_res(      magma_int_t num_rows,
                        magma_int_t num_vecs, 
                        float *evalues, 
                        float *X,
                        float *R, 
                        float *res );

magma_int_t
magma_slobpcg_maxpy(    magma_int_t num_rows,
                        magma_int_t num_vecs, 
                        float *X,
                        float *Y);


/*/////////////////////////////////////////////////////////////////////////////
    -- MAGMA_SPARSE eigensolvers (Data on GPU)
*/
magma_int_t
magma_slobpcg(          magma_s_sparse_matrix A,
                        magma_s_solver_par *solver_par );




/*/////////////////////////////////////////////////////////////////////////////
    -- MAGMA_SPARSE preconditioners (Data on GPU)
*/

magma_int_t
magma_sjacobisetup(     magma_s_sparse_matrix A, 
                        magma_s_vector b, 
                        magma_s_sparse_matrix *M, 
                        magma_s_vector *c );
magma_int_t
magma_sjacobisetup_matrix(  magma_s_sparse_matrix A, 
                            magma_s_sparse_matrix *M, 
                            magma_s_vector *d );
magma_int_t
magma_sjacobisetup_vector(  magma_s_vector b,  
                            magma_s_vector d, 
                            magma_s_vector *c );

magma_int_t
magma_sjacobiiter(      magma_s_sparse_matrix M, 
                        magma_s_vector c, 
                        magma_s_vector *x,  
                        magma_s_solver_par *solver_par );

magma_int_t
magma_sjacobiiter_precond(      
                        magma_s_sparse_matrix M, 
                        magma_s_vector *x,  
                        magma_s_solver_par *solver_par, 
                        magma_s_preconditioner *precond );

magma_int_t
magma_spastixsetup(     magma_s_sparse_matrix A, magma_s_vector b,
                        magma_s_preconditioner *precond );


magma_int_t
magma_sapplypastix(     magma_s_vector b, magma_s_vector *x, 
                        magma_s_preconditioner *precond );


// CUSPARSE preconditioner

magma_int_t
magma_scuilusetup( magma_s_sparse_matrix A, magma_s_preconditioner *precond );

magma_int_t
magma_sapplycuilu_l( magma_s_vector b, magma_s_vector *x, 
                    magma_s_preconditioner *precond );
magma_int_t
magma_sapplycuilu_r( magma_s_vector b, magma_s_vector *x, 
                    magma_s_preconditioner *precond );


magma_int_t
magma_scuiccsetup( magma_s_sparse_matrix A, magma_s_preconditioner *precond );

magma_int_t
magma_sapplycuicc_l( magma_s_vector b, magma_s_vector *x, 
                    magma_s_preconditioner *precond );
magma_int_t
magma_sapplycuicc_r( magma_s_vector b, magma_s_vector *x, 
                    magma_s_preconditioner *precond );


// block-asynchronous iteration

magma_int_t
magma_sbajac_csr(   magma_int_t localiters,
                    magma_s_sparse_matrix D,
                    magma_s_sparse_matrix R,
                    magma_s_vector b,
                    magma_s_vector *x );

/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_SPARSE utility function definitions
*/

magma_int_t
magma_s_spmv(           float alpha, 
                        magma_s_sparse_matrix A, 
                        magma_s_vector x, 
                        float beta, 
                        magma_s_vector y );

magma_int_t
magma_s_spmv_shift(     float alpha, 
                        magma_s_sparse_matrix A, 
                        float lambda,
                        magma_s_vector x, 
                        float beta, 
                        magma_int_t offset,     
                        magma_int_t blocksize,
                        magma_index_t *add_vecs, 
                        magma_s_vector y );

magma_int_t
magma_scuspmm(          magma_s_sparse_matrix A, 
                        magma_s_sparse_matrix B, 
                        magma_s_sparse_matrix *AB );

magma_int_t
magma_scuspaxpy(        float *alpha, magma_s_sparse_matrix A, 
                        float *beta, magma_s_sparse_matrix B, 
                        magma_s_sparse_matrix *AB );

magma_int_t
magma_s_precond(        magma_s_sparse_matrix A, 
                        magma_s_vector b, magma_s_vector *x,
                        magma_s_preconditioner *precond );

magma_int_t
magma_s_solver(         magma_s_sparse_matrix A, magma_s_vector b, 
                        magma_s_vector *x, magma_sopts *zopts );

magma_int_t
magma_s_precondsetup( magma_s_sparse_matrix A, magma_s_vector b, 
                      magma_s_preconditioner *precond );

magma_int_t
magma_s_applyprecond( magma_s_sparse_matrix A, magma_s_vector b, 
                      magma_s_vector *x, magma_s_preconditioner *precond );


magma_int_t
magma_s_applyprecond_left( magma_s_sparse_matrix A, magma_s_vector b, 
                      magma_s_vector *x, magma_s_preconditioner *precond );


magma_int_t
magma_s_applyprecond_right( magma_s_sparse_matrix A, magma_s_vector b, 
                      magma_s_vector *x, magma_s_preconditioner *precond );

magma_int_t
magma_s_initP2P(        magma_int_t *bandwidth_benchmark,
                        magma_int_t *num_gpus );




/* ////////////////////////////////////////////////////////////////////////////
   -- MAGMA_SPARSE BLAS function definitions
*/
magma_int_t 
magma_sgecsrmv(        magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       float alpha,
                       float *d_val,
                       magma_index_t *d_rowptr,
                       magma_index_t *d_colind,
                       float *d_x,
                       float beta,
                       float *d_y );

magma_int_t 
magma_sgecsrmv_shift(  magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       float alpha,
                       float lambda,
                       float *d_val,
                       magma_index_t *d_rowptr,
                       magma_index_t *d_colind,
                       float *d_x,
                       float beta,
                       int offset,
                       int blocksize,
                       magma_index_t *add_rows,
                       float *d_y );

magma_int_t 
magma_smgecsrmv(       magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       magma_int_t num_vecs,
                       float alpha,
                       float *d_val,
                       magma_index_t *d_rowptr,
                       magma_index_t *d_colind,
                       float *d_x,
                       float beta,
                       float *d_y );

magma_int_t 
magma_sgeellmv(        magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       magma_int_t nnz_per_row,
                       float alpha,
                       float *d_val,
                       magma_index_t *d_colind,
                       float *d_x,
                       float beta,
                       float *d_y );

magma_int_t 
magma_sgeellmv_shift(  magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       magma_int_t nnz_per_row,
                       float alpha,
                       float lambda,
                       float *d_val,
                       magma_index_t *d_colind,
                       float *d_x,
                       float beta,
                       int offset,
                       int blocksize,
                       magma_index_t *add_rows,
                       float *d_y );


magma_int_t 
magma_smgeellmv(       magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       magma_int_t num_vecs,
                       magma_int_t nnz_per_row,
                       float alpha,
                       float *d_val,
                       magma_index_t *d_colind,
                       float *d_x,
                       float beta,
                       float *d_y );


magma_int_t 
magma_sgeelltmv(       magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       magma_int_t nnz_per_row,
                       float alpha,
                       float *d_val,
                       magma_index_t *d_colind,
                       float *d_x,
                       float beta,
                       float *d_y );

magma_int_t 
magma_sgeelltmv_shift( magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       magma_int_t nnz_per_row,
                       float alpha,
                       float lambda,
                       float *d_val,
                       magma_index_t *d_colind,
                       float *d_x,
                       float beta,
                       int offset,
                       int blocksize,
                       magma_index_t *add_rows,
                       float *d_y );


magma_int_t 
magma_smgeelltmv(      magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       magma_int_t num_vecs,
                       magma_int_t nnz_per_row,
                       float alpha,
                       float *d_val,
                       magma_index_t *d_colind,
                       float *d_x,
                       float beta,
                       float *d_y );

magma_int_t 
magma_sgeellrtmv(      magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       magma_int_t nnz_per_row,
                       float alpha,
                       float *d_val,
                       magma_index_t *d_colind,
                       magma_index_t *d_rowlength,
                       float *d_x,
                       float beta,
                       float *d_y,
                       magma_int_t num_threads,
                       magma_int_t threads_per_row );

magma_int_t 
magma_sgesellcmv(      magma_trans_t transA,
                       magma_int_t m, magma_int_t n,
                       magma_int_t blocksize,
                       magma_int_t slices,
                       magma_int_t alignment,
                       float alpha,
                       float *d_val,
                       magma_index_t *d_colind,
                       magma_index_t *d_rowptr,
                       float *d_x,
                       float beta,
                       float *d_y );

magma_int_t
magma_sgesellpmv(  magma_trans_t transA,
                    magma_int_t m, magma_int_t n,
                    magma_int_t blocksize,
                    magma_int_t slices,
                    magma_int_t alignment,
                    float alpha,
                    float *d_val,
                    magma_index_t *d_colind,
                    magma_index_t *d_rowptr,
                    float *d_x,
                    float beta,
                    float *d_y );

magma_int_t
magma_smgesellpmv( magma_trans_t transA,
                    magma_int_t m, magma_int_t n,
                    magma_int_t num_vecs,
                    magma_int_t blocksize,
                    magma_int_t slices,
                    magma_int_t alignment,
                    float alpha,
                    float *d_val,
                    magma_index_t *d_colind,
                    magma_index_t *d_rowptr,
                    float *d_x,
                    float beta,
                    float *d_y );

magma_int_t
magma_smgesellpmv_blocked( magma_trans_t transA,
                    magma_int_t m, magma_int_t n,
                    magma_int_t num_vecs,
                    magma_int_t blocksize,
                    magma_int_t slices,
                    magma_int_t alignment,
                    float alpha,
                    float *d_val,
                    magma_index_t *d_colind,
                    magma_index_t *d_rowptr,
                    float *d_x,
                    float beta,
                    float *d_y );


magma_int_t
magma_smergedgs(        magma_int_t n, 
                        magma_int_t ldh,
                        magma_int_t k, 
                        float *v, 
                        float *r,
                        float *skp );

magma_int_t
magma_scopyscale(       int n, 
                        int k,
                        float *r, 
                        float *v,
                        float *skp );

magma_int_t
magma_snrm2scale(      int m, 
                        float *r, int lddr, 
                        float *drnorm);


magma_int_t
magma_sjacobisetup_vector_gpu( int num_rows, 
                               float *b, 
                               float *d, 
                               float *c,
                               float *x );

magma_int_t
magma_sjacobi_diagscal(         int num_rows, 
                                float *b, 
                                float *d, 
                                float *c);

magma_int_t
magma_sjacobisetup_diagscal( magma_s_sparse_matrix A, magma_s_vector *d );


magma_int_t
magma_sbicgmerge1(  int n, 
                    float *skp,
                    float *v, 
                    float *r, 
                    float *p );


magma_int_t
magma_sbicgmerge2(  int n, 
                    float *skp, 
                    float *r,
                    float *v, 
                    float *s );

magma_int_t
magma_sbicgmerge3(  int n, 
                    float *skp, 
                    float *p,
                    float *s,
                    float *t,
                    float *x, 
                    float *r );
magma_int_t
magma_sbicgmerge4(  int type, 
                    float *skp );

magma_int_t
magma_scgmerge_spmv1(  
                 magma_s_sparse_matrix A,
                 float *d1,
                 float *d2,
                 float *d_d,
                 float *d_z,
                 float *skp );

magma_int_t
magma_scgmerge_xrbeta(  
                 int n,
                 float *d1,
                 float *d2,
                 float *x,
                 float *r,
                 float *d,
                 float *z, 
                 float *skp );

magma_int_t
magma_smdotc(       magma_int_t n, 
                    magma_int_t k, 
                    float *v, 
                    float *r,
                    float *d1,
                    float *d2,
                    float *skp );

magma_int_t
magma_sgemvmdot(    int n, 
                    int k, 
                    float *v, 
                    float *r,
                    float *d1,
                    float *d2,
                    float *skp );

magma_int_t
magma_sbicgmerge_spmv1(  
                 magma_s_sparse_matrix A,
                 float *d1,
                 float *d2,
                 float *d_p,
                 float *d_r,
                 float *d_v,
                 float *skp );

magma_int_t
magma_sbicgmerge_spmv2(  
                 magma_s_sparse_matrix A,
                 float *d1,
                 float *d2,
                 float *d_s,
                 float *d_t,
                 float *skp );

magma_int_t
magma_sbicgmerge_xrbeta(  
                 int n,
                 float *d1,
                 float *d2,
                 float *rr,
                 float *r,
                 float *p,
                 float *s,
                 float *t,
                 float *x, 
                 float *skp );

magma_int_t
magma_sbcsrswp(  magma_int_t n,
                 magma_int_t size_b, 
                 magma_int_t *ipiv,
                 float *x );

magma_int_t
magma_sbcsrtrsv( magma_uplo_t uplo,
                 magma_int_t r_blocks,
                 magma_int_t c_blocks,
                 magma_int_t size_b, 
                 float *A,
                 magma_index_t *blockinfo,   
                 float *x );

magma_int_t
magma_sbcsrvalcpy(  magma_int_t size_b, 
                    magma_int_t num_blocks, 
                    magma_int_t num_zero_blocks, 
                    float **Aval, 
                    float **Bval,
                    float **Bval2 );

magma_int_t
magma_sbcsrluegemm( magma_int_t size_b, 
                    magma_int_t num_block_rows,
                    magma_int_t kblocks,
                    float **dA,  
                    float **dB,  
                    float **dC );

magma_int_t
magma_sbcsrlupivloc( magma_int_t size_b, 
                    magma_int_t kblocks,
                    float **dA,  
                    magma_int_t *ipiv );

magma_int_t
magma_sbcsrblockinfo5(  magma_int_t lustep,
                        magma_int_t num_blocks, 
                        magma_int_t c_blocks, 
                        magma_int_t size_b,
                        magma_index_t *blockinfo,
                        float *val,
                        float **AII );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif /* MAGMASPARSE_S_H */
