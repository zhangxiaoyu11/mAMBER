#! /bin/bash																												

########################################################################################
#
# Calculate average electron density map and structure factors from an MD trajectory.
# Returns:
#			md_avg.map - average ASU electron density map in ccp4 map format
#			md_avg.mtz - structure factors of md_avg.mtz
#
# Author v1:James Holton
#        v2: Pawel Janowski
#
########################################################################################


# SET VARIABLES
# all pdb files (including experimental) need to have CRYST record  
sim_pdb_dir=/net/casegroup2/u2/pjanowsk/Case/4lzt/RunSi/average_density/PDBData/  #directory with simulation frames (pdb format, must have element column, )
exp_pdb=solvpep1.pdb								  #experimental structure (pdb format, cryst, remove EPW, etc)
cell=( 27.240 31.870 34.230 88.52 108.53 111.89 )  #experimental unit cell
MD_mult=( 3 2 2 )								  #supercell multiplicity (x y z)
SG="P1"											  #ccp4 space group symbol
GRID="GRID 96 108 120"							  #grid spacing for SFALL map 
B=15											  #flat B-factor added to frames to avoid singularity in SFALL map calculation
align=0											  #set to 1 if request correlation alignment
#~ sed -i 's/CRSYT/CRYST/g' ${sim_pdb_dir}/*  #Bug in ptraj prints CRSYT instead of CRYST. Comment out to save time if not necessary.




#======================================================================#
#======================================================================#

function MakeMap {
	####################################################################
	#                                                                  #
	# CALCULATE ED MAP                                                 #
	#                                                                  #
	####################################################################
	
	#INPUT VARIABLES
	local pdb2map=$1 				#input pdb file name
	local mapname=$2				#output map name

	# MAKE PARAMETER FILE
	cat << EOF > params.txt
MDMULT ${MD_mult[*]}
SHIFT ${shift}
SCALE 1 1 1
CELL ${cell[*]}
BFAC $B
EOF
	# REFORMAT PDB FILE (unit cell cryst and space group, applied shift and scale to 
	#          coordinates, set b-factore, set occupancy=1, reformat atom
	#          names so there is no chance SFALL could confuse hydrogen with mercury)
	cat params.txt $pdb2map |\
	awk 'BEGIN{sx=sy=sz=1;B=20} \
		   /^MDMULT/{na=$2;nb=$3;nc=$4;next}\
		   /^SHIFT/{dx=$2;dy=$3;dz=$4;next}\
		   /^SCALE/ && $2*$3*$4>0{sx=$2;sy=$3;sz=$4;next}\
		   /^BFAC/{B=$2;next}\
		   /^CELL/{a=$2;b=$3;c=$4;al=$5;be=$6;ga=$7;next}\
		   /^CRYST/{a=$2/na*sx;b=$3/nb*sy;c=$4/nc*sz;al=$5;be=$6;ga=$7;sg=(substr($0,56,15));\
			 printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %15s\n",a,b,c,al,be,ga,sg}\
		  /^ATOM/{\
			RES= substr($0, 18, 9);\
			X = (substr($0, 31, 8)+dx)*sx;\
			Y = (substr($0, 39, 8)+dy)*sy;\
			Z = (substr($0, 47, 8)+dz)*sz;\
			Ee = $NF;\
		printf("ATOM %6d %2s   %9s    %8.3f%8.3f%8.3f  1.00%6.2f%12s\n",++n,Ee,RES,X,Y,Z,B,Ee);}\
		END{print "END"}' |\
		cat > sfallme.pdb

	# GET UNIT CELL PARAMETERS
	pcell=( `head -1 sfallme.pdb | awk '{print $2,$3,$4,$5,$6,$7}'` )
	echo cell is ${pcell[*]}
	
	# ADD PDB SCALE MATRIX INFORMATION
	pdbset xyzin sfallme.pdb xyzout new.pdb << EOF > /dev/null
SPACE $SG
CELL ${pcell[*]}
EOF
	mv new.pdb sfallme.pdb

#~ ###########
	#~ ####################################
	#~ tmp=${pdb2map##$sim_pdb_dir}
	#~ tmp2=${tmp%.pdb}_SC.pdb
	#~ cp sfallme.pdb SC_shift/$tmp2
	#~ #####################################
#~ ############


	# COLLAPSE ALL UNIT CELLS INTO ONE
	# Note: this is not necessary as sfall will do it automatically if
	#       given the unit cell and space group
	# convert to fractionals
	#~ coordconv xyzin sfallme.pdb xyzout new.xyz << EOF > /dev/null
#~ INPUT PDB
#~ OUTPUT FRAC
#~ END
#~ EOF
	#~ # Change fractionals to be in range [0,1]
	#~ cat new.xyz |\
	  #~ awk '{fx=$2;while(fx<0)++fx;while(fx>1)--fx;\
			#~ fy=$3;while(fy<0)++fy;while(fy>1)--fy;\
			#~ fz=$4;while(fz<0)++fz;while(fz>1)--fz;\
			#~ printf("%5d%10.5f%10.5f%10.5f%s\n",\
			#~ $1,fx,fy,fz,substr($0,36))}' |\
	  #~ cat > sfallme.xyz
	#~ # Convert back to cartesian
	#~ coordconv xyzin sfallme.xyz xyzout sfallme.pdb << EOF > /dev/null
#~ CELL ${pcell[*]}
#~ INPUT FRAC
#~ OUTPUT PDB
#~ END
#~ EOF
	
	
#~ ###########
	#~ ####################################
	#~ tmp=${pdb2map##$sim_pdb_dir}
	#~ tmp2=${tmp%.pdb}_ASU.pdb
	#~ cp sfallme.pdb ASU_shift/$tmp2
	#~ #####################################
#~ ############


	# CALCULATE MAP
	sfall xyzin sfallme.pdb mapout ${mapname} << EOF > /dev/null
mode atmmap
CELL ${pcell[*]}
SYMM $SG
FORMFAC NGAUSS 5
$GRID
EOF
	#CLEAN
	rm -f params.txt 
	rm -f new.xyz 
	rm -f sfallme.xyz 
	#~ rm -f sfallme.pdb 
}


function AlignMaps {
	####################################################################
	#																   #
	# ALIGN TWO ED MAPS BY CONVOLUTION IN RECIPROCAL SPACE			   #
	#																   #
	####################################################################

	# The optimal x y z shift will be stored in $shift environmental variable
	
	# input variables
	right_map=$1 						# reference map (usually experimental)
	wrong_map=$2						# map to be shifted
	reso=1								# resolution limit for SF calculation
	tempfile="$(mktemp -u tempXXXX_)"   # tempfile prefix
	
	# calculate SF from experimental map
	sfall mapin $right_map hklout ${tempfile}right.mtz << EOF > /dev/null
mode sfcalc mapin
resolution $reso
EOF

	# calculate SF from simulation map
	sfall mapin $wrong_map hklout ${tempfile}wrong.mtz << EOF > /dev/null
mode sfcalc mapin
resolution $reso
EOF

	# use deconvolution to find optimal shift
	# calculate Fexp/Fsim which are the SF's of the correlation function
	rm -f ${tempfile}del.mtz
	sftools << EOF > /dev/null
read ${tempfile}right.mtz
read ${tempfile}wrong.mtz
set labels
Fright
PHIright
Fwrong
PHIwrong
calc ( COL Fq PHIdel ) = ( COL Fright PHIright ) ( COL Fwrong PHIwrong ) /
calc COL W = COL Fq
select COL Fq > 1
calc COL W = 1 COL Fq /
select all
calc F COL Fdel = COL W 0.5 **
write ${tempfile}del.mtz col Fdel PHIdel
y
stop
EOF

	#calculate correlation map from correlation SF's
	fft hklin ${tempfile}del.mtz mapout ${tempfile}del.map << EOF > ${tempfile}.log
labin F1=Fdel PHI=PHIdel
reso $reso
EOF

	# make sure that we define "sigma" for the unmasked map
	echo "scale sigma 1 0" |\
	mapmask mapin ${tempfile}del.map mapout ${tempfile}zscore.map  > /dev/null

	#find highest peak on correlation map
	peakmax mapin ${tempfile}zscore.map xyzout ${tempfile}peak.pdb << EOF > ${tempfile}.log
numpeaks 10
EOF

	#print results
	cat ${tempfile}peak.pdb | awk '/^ATOM/{\
	 print substr($0,31,8),substr($0,39,8),substr($0,47,8),"   ",substr($0,61)+0}' |\
	 awk 'NR==1{max=$4} $4>max/3' | cat > ${tempfile}best_shift.txt
	zscore=`awk '{print $4;exit}' ${tempfile}best_shift.txt`
	echo "z-score: $zscore"
	shift=`awk '{print $1,$2,$3;exit}' ${tempfile}best_shift.txt`

	#clean
	rm -f ${tempfile}best_shift.txt
	rm -f ${tempfile}zscore.map
	rm -f ${tempfile}peak.pdb
	rm -f ${tempfile}del.map
	rm -f ${tempfile}.log
	rm -f ${tempfile}right.mtz
	rm -f ${tempfile}wrong.mtz
	rm -f ${tempfile}del.mtz
}


function AverageDensity {
	####################################################################
	#																   #
	# CALCULATE AVERAGE DENSITY MAP AND ITS STRUCTURE FACTORS   	   #
	#																   #
	####################################################################
		
	# SUM ALL MAPS
	echo "Summing maps"
	rm -f sum.map
	for ((frame=1; frame<=$frames; frame++)); do
		frame_map=`ls -1 maps | sed $frame'q;d'`
		if [ ! -e sum.map ]; then
			cp maps/${frame_map} sum.map

		else
			mapmask mapin1 sum.map mapin2 maps/${frame_map} mapout new.map <<EOF | grep -e "Mean density" |head -1
maps add
EOF
			mv new.map sum.map
		fi
	done
	
	#NORMALIZE SUM MAP
	scale=`echo $frames ${MD_mult[*]} $symops | awk '{print $1*$2*$3*$4*$5}'`	
	mapmask mapin sum.map mapout md_avg.map <<EOF | grep -e "Mean density"
scale factor `echo "1/$scale" | bc -l` 0
EOF
	echo -e "=======================\n"
	echo "Scaling map by 1/$scale"
	echo | mapdump mapin md_avg.map | egrep dens
	
	# CALCULATE SF FROM AVERAGE DENISTY MAP
	sfall mapin md_avg.map hklout md_avg_1.mtz << EOF > /dev/null
MODE SFCALC MAPIN
EOF
	# EDIT MTZ FILE TO INCLUDE SIGFP AND FOM COLUMNS
	sftools << EOF > /dev/null
read md_avg_1.mtz
set labels
FP
PHI
calc Q COL SIGFP = 0.1
calc W COL FOM = 0.9
write md_avg_2.mtz col FP SIGFP PHI FOM
y
stop
EOF
	# REMOVE THE BFACTOR THAT WAS ADDED FOR SF CALCULATION TO AVOID SINGULARITY
	cad hklin1 md_avg_2.mtz hklout md_avg.mtz << EOF >/dev/null
scale file 1 1 -$B
labin file 1 all
EOF
	# CLEAN
	rm -f md_avg_1.mtz 
	rm -f md_avg_2.mtz
	rm -f sfall.mtz
	rm -f sum.map
}


########################################################################
#																	   #
# 	MAIN															   #
#																	   #
########################################################################

# SET-UP
frames=`ls -1 $sim_pdb_dir | wc -l`								  #no. of frames
symops=`awk -v SG=$SG '$4==SG{print $2;exit}' ${CLIBD}/symop.lib` #no. of sym operations
shift=" 0 0 0 "															
rm -rf shifts.txt maps
mkdir -p maps


# CALCULATE EXPERIMENTAL MAP
if [ $align == 1 ]; then
	echo "Calculating experimental map."
	MakeMap $exp_pdb right.map
	echo -e "=======================\n"
fi

# FOR EACH TRAJECTORY FRAME
start_time=`echo "puts [clock clicks -milliseconds]" | tclsh`
for ((frame=1; frame<=$frames; frame++)); do
    shift=""  #set shift to 0 0 0
	sim_pdb=`ls -1 $sim_pdb_dir | sed $frame'q;d'`
	
	# CALCULATE OPTIMAL SHIFT BY FOURIER CONVOLUTION
	if [ $align == 1 ]; then
		echo "Aligning ${sim_pdb_dir}${sim_pdb}"
		MakeMap ${sim_pdb_dir}${sim_pdb} wrong.map
		AlignMaps right.map wrong.map
		echo "$sim_pdb $shift" | awk '{printf("%-15s  %10s  %10s  %10s\n",$1,$2,$3,$4) }' | tee -a shifts.txt
	fi
	
	# CALCULATE MAP WITH SHIFT ($shift is now set to the optimal value)
	echo "Calculating map: ${sim_pdb_dir}${sim_pdb}"
	MakeMap ${sim_pdb_dir}${sim_pdb} maps/${sim_pdb%.*}.map
	
	# REPORT RUNNING TIME
	now_time=`echo "puts [clock clicks -milliseconds]" | tclsh`
	sofar=`echo $now_time $start_time | awk '{print ($1-$2)/1000}'`
	togo=`echo $sofar $frame $frames | awk '{print int($1/$2*$3-$1)}'`
	finish=`echo "puts [clock format [expr [clock seconds] + $togo]]" | tclsh`
	echo "Expect to finish at: $finish"
	echo -e "=======================\n"
done 

# CALCULATE SIMULATION AVERAGE ELECTRON DENSITY
AverageDensity





