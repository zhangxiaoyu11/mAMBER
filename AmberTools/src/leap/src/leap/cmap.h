typedef struct CMNT_t {
    char *record;
    struct CMNT_t *next;
} CMNT;

typedef char WRD[8];

typedef struct {
    char title[256];
    WRD *reslist;
    int nres;
    CMNT *cmnt;
    int resolution;
    double *map;
    int active;
    int mapidx;
    int residx[5];
    WRD atmname[5];
} CMAP;

typedef struct {
    WRD res;
    int atoms[5];
    int mapid;
} PHIPSI;

typedef struct CMAPLST_t {
    CMAP *cmap;
    struct CMAPLST_t *next;
} CMAPLST;

extern CMAP *cmap;
extern CMAPLST *cmaplst;
extern int mapnum;
