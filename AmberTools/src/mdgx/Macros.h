#ifndef MacroDefinitions
#define MacroDefinitions

/***=======================================================================***/
/*** Useful functions                                                      ***/
/***=======================================================================***/
#define SIGN(x)            ((x >= 0.0) ? 1.0 : -1.0)

#define MAX(a, b)          ((a) > (b) ? (a) : (b))

#define MIN(a, b)          ((a) > (b) ? (b) : (a))

#define SWAP(x, y, tmp)    (tmp) = (x); (x) = (y); (y) = (tmp)

/***=======================================================================***/
/*** Powers of 5, for integer encoding of critical values                  ***/
/***=======================================================================***/
#define POWFIVE1           5
#define POWFIVE2           25
#define POWFIVE3           125
#define POWFIVE4           625
#define POWFIVE5           3125
#define POWFIVE6           15625
#define POWFIVE7           78125
#define POWFIVE8           390625
#define POWFIVE9           1953125
#define POWFIVE10          9765625
#define POWFIVE11          48828125
#define POWFIVE12          244140625

#endif
