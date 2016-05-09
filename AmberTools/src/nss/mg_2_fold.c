#include <stdio.h>
#include <avs/avs.h>
#include <avs/port.h>
#include <avs/field.h>

#include "nabcode.h"
 
#define	MAT_ALLOC(n)	( ( MATRIX_T * )malloc( ( n )*sizeof( MATRIX_T ) ) )
#define	BPN	14

char	*MAT_getsyminfo();

static	int	pid = -1;
static	char	mname[ 256 ] = "";

/* *****************************************/
/*  Module Description                     */
/* *****************************************/
int mg_2_fold_desc()
{

	int in_port, out_port, param, iresult;
	extern int mg_2_fold_compute();

	AVSset_module_name( "mg_2_fold", MODULE_DATA );

	/* Input Port Specifications               */
	in_port = AVScreate_input_port( "InMatrices", 
		"field 1D 1-space 1-vector uniform byte", OPTIONAL );

	/* Output Port Specifications              */
	out_port = AVScreate_output_port( "OutMatrices", 
		"field 1D 1-space 1-vector uniform byte" );

	out_port = AVScreate_output_port( "o_noID", "integer" );
	AVSset_output_flags( out_port, INVISIBLE );

	out_port = AVScreate_output_port( "o_Center", "string" );
	AVSset_output_flags( out_port, INVISIBLE );

	out_port = AVScreate_output_port( "o_Axis", "string" );
	AVSset_output_flags( out_port, INVISIBLE );

	param = AVSadd_parameter( "noID", "boolean", 0, 0, 1 );

	param = AVSadd_parameter( "Axestype", "choice", "relative axes",
		"relative axes:absolute axes", ":" );
	AVSconnect_widget( param, "radio_buttons" );

	param = AVSadd_parameter( "Center", "string", "0 0 0", "", ":" );
	AVSadd_parameter_prop( param, "width", "integer", 4 );
	AVSconnect_widget( param, "typein" );

	param = AVSadd_parameter( "Axis", "string", "1 0 0", "", ":" );
	AVSadd_parameter_prop( param, "width", "integer", 4 );
	AVSconnect_widget( param, "typein" );

	param = AVSadd_parameter( "Name", "string", "", "", "" );
	AVSadd_parameter_prop( param, "width", "integer", 4 );
	AVSconnect_widget( param, "typein" );

	AVSset_compute_proc( mg_2_fold_compute );
	return( 1 );
}
 
/* *****************************************/
/* Module Compute Routine                  */
/* *****************************************/
int mg_2_fold_compute( InMatrices, OutMatrices,
	o_noID, o_Center, o_Axis,
	noID, Axestype, Center, Axis, Name
)
AVSfield_char	*InMatrices;
AVSfield_char	**OutMatrices;
int		*o_noID;
char		**o_Center;
char		**o_Axis;
int		noID;
char		*Axestype;
char		*Center;
char		*Axis;
char		*Name;
{
	int	noid, create;
	int	i, j, k;
	POINT_T	pts[ 4 ];
	REAL_T	ang;
	int	cnt;
	int	m, ms, mi, mo;
	int	nmats, simats, nimats, nomats;
	MATRIX_T	*mats, *imats, *omats;
	int dims0[ 1 ];
	char	*sip;
	char	m_si[ 1000 ];
	int	s_sip, s_o_sip;
	char	*dp;

	if( pid == -1 )
		pid = getpid();

	if( !*Name ){
		if( !*mname )
			sprintf( mname, "m%d", pid );
	}else
		strcpy( mname, Name );
	AVSmodify_parameter( "Name", AVS_VALUE, mname, "", "" );

	if( *OutMatrices )
		AVSfield_free(*OutMatrices);

	if( *o_Center )
		free( *o_Center );
	if( *o_Axis )
		free( *o_Axis );

	noid = noID;
	create = 1;
	ang = 180.0;
	cnt = 2;

	nmats = nimats = nomats = 0;
	mats = imats = omats = NULL;

	for( i = 0; i < 4; i++ )
		for( j = 0; j < 3; j++ )
			pts[ i ][ j ] = 0.0; 

	if( *Center == '\0' )
		return( 0 );
	if(sscanf(Center, "%lf %lf %lf",&pts[0][0],&pts[0][1],&pts[0][2]) != 3){
		AVSerror( "Center must be specified." );
		return( 0 );
	}

	if( *Axis == '\0' )
		return( 0 );
	if( sscanf(Axis, "%lf %lf %lf",&pts[1][0],&pts[1][1],&pts[1][2]) != 3 ){
		AVSerror( "Axis must be specified." );
		return( 0 );
	}

	nmats = cnt;
	mats = MAT_ALLOC( nmats );
	if( mats == NULL ){
		AVSerror( "Allocation of mats failed." );
		return( 0 );
	}
	if( !strcmp( Axestype, "relative axes" ) )
		NAB_ptadd( pts[1], pts[0], pts[1] );
	MAT_cyclic( pts, &ang, &cnt, mats );

	if( !InMatrices ){
		nimats = 0;
		sip = NULL;
		s_sip = 0;
	}else{
		simats = MAT_count( InMatrices->data );
		imats = MAT_ALLOC( simats );
		if( imats == NULL ){
			AVSerror( "Allocation of mats failed." );
			free( mats );
			return( 0 );
		}
		if( ( nimats = MAT_sscan(InMatrices->data,simats,imats))==0 ){
			AVSerror( "InMatrices is not a valid matrix field." );
			free( imats );
			free( mats );
			return( 0 );
		}
		sip = MAT_getsyminfo( InMatrices->data );
		s_sip = strlen( sip );
		create = 0;
	}

	dp = m_si;
	sprintf( dp, "#S{ cyclic_2 %d %s\n", pid, mname );
	dp += strlen( dp );
	sprintf( dp, "#S+   noid     %s\n", noid ? "true" : "false" );
	dp += strlen( dp );
	sprintf( dp, "#S+   axestype %s\n",
		!strcmp(Axestype, "relative axes") ? "relative" : "absolute" );
	dp += strlen( dp );
	sprintf( dp, "#S+   center   %e %e %e\n",
		pts[0][0], pts[0][1], pts[0][2] );
	dp += strlen( dp );
	sprintf( dp, "#S+   axis     %e %e %e\n",
		pts[1][0], pts[1][1], pts[1][2] );
	s_o_sip = strlen( m_si ) + strlen( "#S}\n" ) + s_sip;

	if( create || nimats == 0 ){
		if( noid ){
			ms = 1;
			nomats = nmats - 1;
		}else{
			ms = 0;
			nomats = nmats;
		}
		omats = ( MATRIX_T * )mats[ ms ];
	}else{
		if( noid ){
			ms = 1;
			nomats = nimats * ( nmats - 1 );
		}else{
			ms = 0;
			nomats = nimats * nmats;
		}
		omats = MAT_ALLOC( nomats );
		if( omats == NULL ){
			AVSerror( "Allocation of omats failed." );
			free( imats );
			free( mats );
			return( 0 );
		}
		mo = 0;
		for( m = ms; m < nmats; m++ ){
			for( mi = 0; mi < nimats; mi++, mo++ ){
				NAB_matcpy( omats[ mo ],
					MAT_concat( imats[ mi ], mats[ m ] ) );
			}
		}
	}

	*o_noID = noID;
	*o_Center = ( char * )malloc( strlen( Center ) + 1 );
	strcpy( *o_Center, Center );
	*o_Axis = ( char * )malloc( strlen( Axis ) + 1 );
	strcpy( *o_Axis, Axis );

	dims0[0] = 16 * BPN * nomats + s_o_sip; 
	*OutMatrices = (AVSfield_char *) AVSdata_alloc(
		"field 1D 1-space 1-vector uniform byte", 
		dims0);
	if (*OutMatrices == NULL) {
		AVSerror( "Allocation of output field failed." );
		if( !create && nimats != 0 )
			free( omats );
		if( imats )
			free( imats );
		free( mats );
		return( 0 );
	}

	dp = ( *OutMatrices )->data;
	strcpy( dp, m_si );
	dp += strlen( m_si );
	if( sip != NULL ){
		strcpy( dp, sip );
		dp += s_sip;
	}
	strcpy( dp, "#S}\n" );
	dp += strlen( "#S}\n" );

	MAT_sprint( dp, nomats, omats );

	if( !create && nimats != 0 )
		free( omats );
	if( imats )
		free( imats );
	free( mats );
 
	return(1);
}
 
/* ***********************************************************************/
/* Initialization for modules contained in this file.                    */
/* ***********************************************************************/
static int ((*mod_list[])()) = {
	mg_2_fold_desc
};
#define NMODS (sizeof(mod_list) / sizeof(char *))

AVSinit_modules()
{
	AVSinit_from_module_list(mod_list, NMODS);
}
