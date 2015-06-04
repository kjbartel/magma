#ifndef AUXILIARY_H
#define AUXILIARY_H

#ifndef max
#define max(a,b) ( ( (a) > (b) ) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ( ( (a) < (b) ) ? (a) : (b))
#endif

#include "zauxiliary.h"
#include "cauxiliary.h"
#include "dauxiliary.h"
#include "sauxiliary.h"

extern int IONE;
extern int ISEED[4];

extern int format[6];
extern int trans[3];
extern int uplo[2];
extern int side[2];
extern int diag[2];
extern char *formatstr[6];
extern char *transstr[3];
extern char *uplostr[2];
extern char *sidestr[2];
extern char *diagstr[2];

extern void *formatmap[6];

int map_CM  (int m, int n, int mb, int nb, int i, int j);
int map_CCRB(int m, int n, int mb, int nb, int i, int j);
int map_CRRB(int m, int n, int mb, int nb, int i, int j);
int map_RCRB(int m, int n, int mb, int nb, int i, int j);
int map_RRRB(int m, int n, int mb, int nb, int i, int j);
int map_RM  (int m, int n, int mb, int nb, int i, int j);

#endif /* AUXILIARY_H */
