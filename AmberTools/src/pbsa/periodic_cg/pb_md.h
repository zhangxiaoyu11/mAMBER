!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ pb_md.h PB variables shared with calling routines

   logical pbverbose
   logical pbprint
   logical pbgrid
   logical pbinit
   common /pb_mdl/ pbverbose
   common /pb_mdl/ pbprint
   common /pb_mdl/ pbgrid
   common /pb_mdl/ pbinit

   integer npbstep
   integer nsaslag
   integer npbgrid
   integer nsnbr
   integer nsnba
   integer ntnbr
   integer ntnba
   integer dofd
   integer ndofd
   integer dosas
   integer ndosas
   integer dosar
   integer ndosar
   common /pb_mdi/ npbstep
   common /pb_mdi/ nsaslag
   common /pb_mdi/ npbgrid
   common /pb_mdi/ nsnbr
   common /pb_mdi/ nsnba
   common /pb_mdi/ ntnbr
   common /pb_mdi/ ntnba
   common /pb_mdi/ dofd
   common /pb_mdi/ ndofd
   common /pb_mdi/ dosas
   common /pb_mdi/ ndosas
   common /pb_mdi/ dosar
   common /pb_mdi/ ndosar

!+ End of pb_md.h
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
