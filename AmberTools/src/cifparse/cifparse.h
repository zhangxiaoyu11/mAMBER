/*
 *   MxNameLen applies to block, category, item and keyword names...
 *   
 */
#define MxNameLen    		 80 
/*
 *   This is the size of the buffer length allocated to individual
 *   values.
 */

#define CifDefaultSpace           2

#undef   YYLMAX
#define  YYLMAX       1024
#if 0
#undef   YYMAXDEPTH  
#define  YYMAXDEPTH  20000
#undef   YYINITDEPTH 
#define  YYINITDEPTH  1000
#define  YYPRINT(file, type, value) cifpprint(file, type, value)
#endif

#define MAXVALUELENGTH   8182

#define TRUE                   1
#define FALSE                  0

typedef struct NdbCifRowFormat {
  char **columns;
} NdbCifRowFormat;

typedef struct NdbCifCategoryFormat {
  int numCol;
  int allCol;
  int curCol;
  int allRow;
  int numRow;
  int curRow;
  char categoryName[MxNameLen];
  char **colNames;
  NdbCifRowFormat *rows;
} NdbCifCategoryFormat;

typedef struct NdbCifDatablockFormat {
  int numCategory; /* Number of categories in this datablock */
  int allCategory; /* Allocated space */
  int curCategory; /* index of the current category */
  char datablockName[MxNameLen]; 
  NdbCifCategoryFormat *categories;
} NdbCifDatablockFormat;  

typedef struct NdbCifDatablocksFormat {
  int numDatablock; /* Number of datablocks in this structure */
  int curDatablock; /* index of the current datablock */
  int allDatablock; /* Allocated space */
  NdbCifDatablockFormat *datablocks;
} NdbCifDatablocksFormat;


char *ndb_cif_copy_item_value();
char ndb_cif_set_null_char();

#ifdef CIF_GLOBAL
	FILE *cifpin;
	char TempKeyword[MxNameLen+1], TempValue[MAXVALUELENGTH+1];
	NdbCifDatablocksFormat cifFiles;
	int  lineNo;      
#else
	extern char TempKeyword[MxNameLen+1], TempValue[MAXVALUELENGTH+1];
	extern FILE *cifpin;
	extern int  lineNo;
	extern NdbCifDatablocksFormat cifFiles;
#endif
