
!  --- following should not need to be modified unless you are adding
!      more variables to the "locmem" style of memory management

!       BC_MEMORY is the size of the MEMORY common block:
! MFC changing indices into IX and X to more rational names
!       I12 = Iibh
!       I14 = Ijbh
!       I16 = Iicbh
!       I18 = Iiba
!       I20 = Ijba
!       I22 = Iicba
#define BC_MEMORY 181

integer       natom,nres,nbonh,nbona,ntheth,ntheta,nphih, &
      nphia,nnb,ntypes,nconp,maxmem,nwdvar,maxnb,nparm, &
      natc,ibelly,natbel,ishake,nmxrs, &
      mxsub,natyp,npdec,i02,i04,i06,i08,i10, &
      iibh,ijbh,iicbh,iiba,ijba,iicba, &
      i24,i26,i28,i30,i32,i34,i36,i38,i40, &
      i42,i44,i46,i48,i50,i52,i54,i56,i58,ibellygp, &
      icnstrgp,i64,i65,i66,i68, &
      i70,i72,i74,i76,i78,i79,i80,i82,i84,i86, &
      icpstinf,icpresst,icptrsct, icpptcnt, &
      l05,l10,l15,lwinv,lpol,lcrd,lforce,l36,lvel,lvel2,l45,l50, &
      lcrdr,l60,l65,lmass,l75,l80,l85,l90,l95,l96,l97,l98,lfrctmp, &
      l105,l110,l115,l120,l125,l130,l135,l140,l145,l150, &
      l165,l170,l175,l180,l185,l186,l187,l188,l189,l190, &
      lcpcrg,lcpene, &
      m02,m04,m06,m08,m10,m12,m14,m16,m18,i01, &
      iifstwt,iifstwr,nrealb,nintb,nholb,npairb,lastr,lasti,lasth, &
      lastpr,lastrst,lastist,nbper,ngper,ndper,ifpert,lpolp, ncopy, &
      imask1,imask2,numadjst,mxadjmsk

! 1     2      3      4      5      6      7      8      9      10
common/memory/ &
 natom ,nres  ,nbonh ,nbona ,ntheth,ntheta,nphih ,                       & ! 7
 nphia ,nnb   ,ntypes,nconp ,maxmem,nwdvar,maxnb ,nparm ,                & !15
 natc  ,ibelly,natbel,ishake,nmxrs ,                                     & !20
 mxsub ,natyp ,npdec ,i02   ,i04   ,i06   ,i08   ,i10,                   & !28
 iibh  ,ijbh  ,iicbh ,iiba  ,ijba  ,iicba ,                              & !34
 i24   ,i26   ,i28   ,i30   ,i32   ,i34   ,i36   ,i38   ,i40   ,         & !43
 i42   ,i44   ,i46   ,i48   ,i50   ,i52   ,i54   ,i56   ,i58   ,ibellygp,& !53
 icnstrgp,i64 ,i65   ,i66   ,i68   ,                                     & !58
 i70   ,i72   ,i74   ,i76   ,i78   ,i79   ,i80   ,i82   ,                & !66
 i84   ,i86   ,                                                          & !68
 icpstinf,icpresst,icptrsct, icpptcnt,                                   & !72
 l05   ,l10   ,l15   ,lwinv ,lpol  ,lcrd  ,lforce,l36   ,lvel  ,lvel2 ,  & !82
 l45   ,l50   ,                                                          & !84
 lcrdr ,l60   ,l65   ,lmass ,l75   ,l80   ,l85   ,l90   ,l95   ,l96   ,  & !94
 l97   ,l98   ,lfrctmp,                                                  & !97
 l105  ,l110  ,l115  ,l120  ,l125  ,l130  ,l135  ,l140  ,l145  ,l150  ,  & !107
 l165  ,l170  ,l175  ,l180  ,l185  ,l186  ,l187  ,l188  ,l189  ,l190  ,  & !117
 lcpcrg,lcpene,                                                          & !119
 m02   ,m04   ,m06   ,m08   ,m10   ,m12   ,m14   ,m16   ,m18   ,i01   ,  & !129
 iifstwt,iifstwr,nrealb,nintb,nholb,npairb,lastr ,lasti ,lasth ,         & !138
 lastpr ,lastrst,lastist,nbper,ngper,ndper,ifpert,lpolp ,ncopy,          & !147
 imask1 ,imask2 ,numadjst,mxadjmsk                                         !151
