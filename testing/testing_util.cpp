/*
 *  -- MAGMA (version 1.3.0) --
 *     Univ. of Tennessee, Knoxville
 *     Univ. of California, Berkeley
 *     Univ. of Colorado, Denver
 *     November 2012
 *
 * @precisions normal z -> c s d
 *
 * Utilities for testing.
 * @author Mark Gates
 **/

#include "testings.h"

// --------------------
// If condition is false, print error message and exit.
// Error message is formatted using printf, using any additional arguments.
extern "C"
void magma_assert( bool condition, const char* msg, ... )
{
    if ( ! condition ) {
        va_list va;
        va_start( va, msg );
        vprintf( msg, va );
        exit(1);
    }
}
