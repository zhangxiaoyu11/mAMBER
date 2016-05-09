/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /home/case/cvsroot/amber11/AmberTools/src/netcdf/src/nctest/add.h,v 9.1 2007/11/15 21:44:49 jmongan Exp $
 *********************************************************************/

#undef PROTO
#ifndef NO_HAVE_PROTOTYPES 
#   define	PROTO(x)	x
#else
#   define	PROTO(x)	()
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* add a dimension to an in-memory netcdf structure */
extern void	add_dim		PROTO((
				       struct netcdf *,
				       struct cdfdim *
				       ));

/* add a variable var to an in-memory netcdf structure */
extern void	add_var		PROTO((
				       struct netcdf *,
				       struct cdfvar *
				       ));

/* add an attribute att to an in-memory netcdf structure */
extern void	add_att		PROTO((
				       struct netcdf *,
				       int,
				       struct cdfatt *
				       ));

/* reset in-memory netcdf structure to empty */
extern void	add_reset	PROTO((
				       struct netcdf *
				       ));

/* delete an attribute att from an in-memory netcdf structure */
extern void	del_att		PROTO((
				       struct netcdf *,
				       int,
				       struct cdfatt *
				       ));

/* keep max record written updated */
extern void	add_data	PROTO((
				       struct netcdf *,
				       int,
				       long [],
				       long []
				       ));

/* display info about variable, for error messages */
extern void	errvar		PROTO((
				       struct netcdf *,
				       struct cdfvar *
				       ));

#ifdef __cplusplus
}
#endif
