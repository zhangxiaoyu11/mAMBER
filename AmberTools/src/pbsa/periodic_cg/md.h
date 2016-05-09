#include "dprec.h"
!-------------BEGIN    md.h  ------------------------------------------------

integer BC_MDI  ! size in integers of common block mdi
integer BC_MDR  ! size in _REAL_s  of common block mdr

! ... integer variables:

integer nrp,nspm,ig,ntx,ntcx,           &!5
      ntxo,ntt,ntp,ntr,init,             &!10
      ntcm,nscm,isolvp,nsolut,klambda,   &!15
      ntc,ntcc,ntf,ntid,ntn,             &!20
      ntnb,nsnb,ndfmin,nstlim,nrc,       &!25
      ntrx,npscal,imin,maxcyc,ncyc,      &!30
      ntmin,irest,jfastw,                &!33
      ibgwat,ienwat,iorwat,              &!36
      iwatpr,nsolw,igb,iyammp,           &!40
      gbsa,vrand,iwrap,nrespa,irespa,nrespai,icfe,  &!47
      rbornstat,ivcap,iconstreff,        &!50
      idecomp,icnstph,ntcnstph,maxdup     !54
parameter (BC_MDI=54)

common/mdi/nrp,nspm,ig, &
      ntx,ntcx,ntxo,ntt,ntp,ntr,init,ntcm,nscm, &
      isolvp,nsolut,ntc,ntcc,ntf,ntid,ntn,ntnb,nsnb,ndfmin, &
      nstlim,nrc,ntrx,npscal,imin,maxcyc,ncyc,ntmin, &
      irest,jfastw,ibgwat,ienwat,iorwat, &
      iwatpr,nsolw,igb,iyammp,gbsa,vrand, &
      iwrap,nrespa,irespa,nrespai,icfe,rbornstat, &
      ivcap,iconstreff,idecomp,klambda,icnstph,ntcnstph,maxdup

! ... floats:

_REAL_ t,dt,temp0,tautp,pres0,comp,taup,temp,tempi,        & !9
      tol,taur,dx0,drms,vlimit,rbtarg(8),tmass,tmassinv,   & !24
      kappa,offset,surften,gamma_ln,extdiel,intdiel,rdt,   & !31
      gbalpha,gbbeta,gbgamma,cut_inner,clambda,saltcon,    & !37
      solvph,rgbmax,fsmax,restraint_wt                       !41 
parameter (BC_MDR=41)

common/mdr/t,dt,temp0,tautp,pres0,comp,taup,temp,tempi, &
      tol,taur,dx0,drms,vlimit,rbtarg,tmass,tmassinv, &
      kappa,offset,surften,gamma_ln,extdiel,intdiel,rdt, &
      gbalpha,gbbeta,gbgamma,cut_inner,clambda,saltcon, &
      solvph,rgbmax,fsmax,restraint_wt 


! ... strings:

character(len=4) iwtnm,iowtnm,ihwtnm
character(len=256) restraintmask,bellymask
common/mds/ restraintmask,bellymask,iwtnm,iowtnm,ihwtnm(2)

!-------------END    md.h  ------------------------------------------------

