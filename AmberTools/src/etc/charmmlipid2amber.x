#!/bin/csh 

###############################################################
#                                                             #
#                    charmmlipid2amber.x                      #
#      Charmm Lipid Builder to AMBER Lipid12  Convertor       #
#                         v1.3                                #
#                          by                                 #
#                    Ross Walker (SDSC)                       #
#                    Age Skjevik (UiB)                        #
#                    Ben Madej (UCSD, SDSC)                   #
#                                                             #
###############################################################

#This script is designed to convert pdb output from the Charmm
#Lipid builder to a pdb that is readable by AMBER's Leap using
#the Lipid11 force field.
#
#Charmm Lipid Builder: http://www.charmm-gui.org/?doc=input/membrane
#
#Revision History
#----------------
#Jun 17 2010 - Version 1.0
#Mar 16 2012 - Version 1.1 - Added all compatible lipids CHARMM-GUI -> LIPID11
#May 02 2012 - Version 1.1.1 - Fixed bug in ion conversion
#Feb 01 2013 - Version 1.2 - Added DLPC and DMPC. Added DMPG. Other LA and MY phospholipids to follow in next version
#Feb 28 2013 - Version 1.3 - Added all LA and MY phospholipids.


#Usage:
#./charmmlipid2amber.x foo_input.pdb foo_output.pdb
#

if ($#argv != 2) then
   echo "Incorrect number of arguments."
   echo "usage:  `basename $0` foo_input.pdb foo_output.pdb"
   exit(1)
endif

set input_filename = $1
set output_filename = $2

if (! -e $input_filename) then
  echo "Input file: $input_filename does not exist."
  exit(1)
endif

echo "               charmmlipid2amber.x              "
echo "Charmm Lipid Builder to AMBER Lipid 12 Convertor"
echo "------------------------------------------------"
echo "                                                "
echo "                    v1.3                        "
echo "                     by                         "
echo "               Ross Walker (SDSC)               "
echo "               Age Skjevik (UiB)                "
echo "               Ben Madej (UCSD, SDSC)           "
echo "                                                "
echo "------------------------------------------------"
echo " "
echo "Input PDB file: $input_filename"
echo "Output PDB file: $output_filename"
echo " "

echo " *** STAGE 1 : STRIP REMARKS AND END STATEMENTS *** "
echo " "
grep -v "REMARK" $input_filename >lipid_tmp
grep -v "END" lipid_tmp >lipid_tmp2
mv lipid_tmp2 lipid_tmp


echo " *** STAGE 2 : STRIP AND SAVE WATER *** "
echo " "
echo " Currently supported water: TIP3"
echo " "

#We need to strip off the water and save it in AMBER format. We will
#re-append it later to the output file.

rm -f water_tmp
grep " TIP3" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found TIP3 Water"
  grep " TIP3" lipid_tmp >water_tmp
  grep -v " TIP3"  lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
  set water_found = 1
else
  echo " Water was not found."
  set water_found = 0
endif
echo " "

if ( $water_found == 1 ) then
  sed 's/ TIP3/ WAT /' water_tmp > water_tmp2
  sed 's/OH2/O  /' water_tmp2 >water_tmp
  rm -f water_tmp2
endif

echo " *** STAGE 3 : STRIP AND SAVE IONS *** "
echo " "
echo " Currently supported ions: CLA, POT, SOD"
echo " "

#Process Ions
#First check for ions in the file.

rm -f ions_tmp
touch ions_tmp

grep "CLA CLA" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found Chlorine Ions"
  grep "CLA CLA" lipid_tmp >>ions_tmp
  grep -v " CLA" lipid_tmp > lipid_tmp2
  mv lipid_tmp2 lipid_tmp
endif

grep "SOD SOD" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found Sodium Ions"
  grep "SOD SOD" lipid_tmp >>ions_tmp
  grep -v " SOD" lipid_tmp > lipid_tmp2
  mv lipid_tmp2 lipid_tmp
endif

grep "POT POT" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found Potassium Ions"
  grep "POT POT" lipid_tmp >>ions_tmp
  grep -v " POT" lipid_tmp > lipid_tmp2
  mv lipid_tmp2 lipid_tmp
endif


if (-z ions_tmp) then
  echo " No ions found."
  set ions_found = 0
  rm -f ions_tmp
else 
  echo " Ions were found."
  set ions_found = 1

  sed 's/CLA CLA/Cl- Cl-/' ions_tmp >ions_tmp2
  sed 's/SOD SOD/Na+ Na+/' ions_tmp2 >ions_tmp
  sed 's/POT POT/K+  K+ /' ions_tmp >ions_tmp2
  mv ions_tmp2 ions_tmp
endif
echo " "

echo " *** STAGE 4 : STRIP INDIVIDUAL LIPIDS TO SEPARATE FILES *** "
echo " "
echo " Currently supported lipids: DPPC DPPE DPPS DPPG DPPA"
echo "                             DOPC DOPE DOPS DOPG DOPA"
echo "                             POPC POPE POPS POPG POPA"
echo "                             DLPC DLPE DLPS DLPG DLPA"
echo "                             DMPC DMPE DMPS DMPG DMPA"
echo "                             CHL1"
echo " "

set lipids_found = 0

#DPPC
grep " DPPC" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found DPPC"
  set dppc_found = 1
  set lipids_found = 1
  grep " DPPC" lipid_tmp >DPPC_tmp
  grep -v " DPPC" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dppc_found = 0
endif

#DPPE
grep " DPPE" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set dppe_found = 1
  set lipids_found = 1
  grep " DPPE" lipid_tmp >DPPE_tmp
  grep -v " DPPE" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dppe_found = 0
endif

#DPPS
grep " DPPS" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set dpps_found = 1
  set lipids_found = 1
  grep " DPPS" lipid_tmp >DPPS_tmp
  grep -v " DPPS" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dpps_found = 0
endif

#DPPG
grep " DPPG" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set dppg_found = 1
  set lipids_found = 1
  grep " DPPG" lipid_tmp >DPPG_tmp
  grep -v " DPPG" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dppg_found = 0
endif

#DPPA
grep " DPPA" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found DPPA"
  set dppa_found = 1
  set lipids_found = 1
  grep " DPPA" lipid_tmp >DPPA_tmp
  grep -v " DPPA" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dppa_found = 0
endif

#DOPC
grep " DOPC" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set dopc_found = 1
  set lipids_found = 1
  grep " DOPC" lipid_tmp >DOPC_tmp
  grep -v " DOPC" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dopc_found = 0
endif

#DOPE
grep " DOPE" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set dope_found = 1
  set lipids_found = 1
  grep " DOPE" lipid_tmp >DOPE_tmp
  grep -v " DOPE" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dope_found = 0
endif

#DOPS
grep " DOPS" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set dops_found = 1
  set lipids_found = 1
  grep " DOPS" lipid_tmp >DOPS_tmp
  grep -v " DOPS" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dops_found = 0
endif

#DOPG
grep " DOPG" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set dopg_found = 1
  set lipids_found = 1
  grep " DOPG" lipid_tmp >DOPG_tmp
  grep -v " DOPG" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dopg_found = 0
endif

#DOPA
grep " DOPA" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set dopa_found = 1
  set lipids_found = 1
  grep " DOPA" lipid_tmp >DOPA_tmp
  grep -v " DOPA" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dopa_found = 0
endif

#POPC
grep " POPC" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set popc_found = 1
  set lipids_found = 1
  grep " POPC" lipid_tmp >POPC_tmp
  grep -v " POPC" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set popc_found = 0
endif

#POPE
grep " POPE" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set pope_found = 1
  set lipids_found = 1
  grep " POPE" lipid_tmp >POPE_tmp
  grep -v " POPE" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set pope_found = 0
endif

#POPS
grep " POPS" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set pops_found = 1
  set lipids_found = 1
  echo " Found POPS"
  grep " POPS" lipid_tmp >POPS_tmp
  grep -v " POPS" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set pops_found = 0
endif

#POPG
grep " POPG" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set popg_found = 1
  set lipids_found = 1
  grep " POPG" lipid_tmp >POPG_tmp
  grep -v " POPG" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set popg_found = 0
endif

#POPA
grep " POPA" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set popa_found = 1
  set lipids_found = 1
  grep " POPA" lipid_tmp >POPA_tmp
  grep -v " POPA" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set popa_found = 0
endif

#DLPC
grep " DLPC" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found DLPC"
  set dlpc_found = 1
  set lipids_found = 1
  grep " DLPC" lipid_tmp >DLPC_tmp
  grep -v " DLPC" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dlpc_found = 0
endif

#DLPE
grep " DLPE" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found DLPE"
  set dlpe_found = 1
  set lipids_found = 1
  grep " DLPE" lipid_tmp >DLPE_tmp
  grep -v " DLPE" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dlpe_found = 0
endif

#DLPS
grep " DLPS" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set dlps_found = 1
  set lipids_found = 1
  grep " DLPS" lipid_tmp >DLPS_tmp
  grep -v " DLPS" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dlps_found = 0
endif

#DLPG
grep " DLPG" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set dlpg_found = 1
  set lipids_found = 1
  grep " DLPG" lipid_tmp >DLPG_tmp
  grep -v " DLPG" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dlpg_found = 0
endif

#DLPA
grep " DLPA" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found DLPA"
  set dlpa_found = 1
  set lipids_found = 1
  grep " DLPA" lipid_tmp >DLPA_tmp
  grep -v " DLPA" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dlpa_found = 0
endif

#DMPC
grep " DMPC" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found DMPC"
  set dmpc_found = 1
  set lipids_found = 1
  grep " DMPC" lipid_tmp >DMPC_tmp
  grep -v " DMPC" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dmpc_found = 0
endif

#DMPE
grep " DMPE" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found DMPE"
  set dmpe_found = 1
  set lipids_found = 1
  grep " DMPE" lipid_tmp >DMPE_tmp
  grep -v " DMPE" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dmpe_found = 0
endif

#DMPS
grep " DMPS" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found DMPS"
  set dmps_found = 1
  set lipids_found = 1
  grep " DMPS" lipid_tmp >DMPS_tmp
  grep -v " DMPS" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dmps_found = 0
endif

#DMPG
grep " DMPG" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found DMPG"
  set dmpg_found = 1
  set lipids_found = 1
  grep " DMPG" lipid_tmp >DMPG_tmp
  grep -v " DMPG" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dmpg_found = 0
endif

#DMPA
grep " DMPA" lipid_tmp >/dev/null
if ( $status == 0 ) then
  echo " Found DMPA"
  set dmpa_found = 1
  set lipids_found = 1
  grep " DMPA" lipid_tmp >DMPA_tmp
  grep -v " DMPA" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set dmpa_found = 0
endif

#CHL
grep " CHL1" lipid_tmp >/dev/null
if ( $status == 0 ) then
  set chl_found = 1
  set lipids_found = 1
  grep " CHL1" lipid_tmp >CHL_tmp
  grep -v " CHL1" lipid_tmp >lipid_tmp2
  mv lipid_tmp2 lipid_tmp
else
  set chl_found = 0
endif

echo " "

#Lipid tmp file should now just contain protein.
if ( -z lipid_tmp ) then
  echo " No protein atoms detected in input pdb file."
  rm -f lipid_tmp
  set proteins_found = 0
else
  echo " Protein atoms were detected in the input pdb file."
  echo -n " Line count for protein section of pdb = "
  wc -l lipid_tmp | awk '{print $1}'
  mv lipid_tmp protein_tmp
  set proteins_found = 1
endif
echo " "

if ( $lipids_found == 0 ) then
  echo " ERROR - No supported lipids found in pdb file."
  exit (1)
endif

echo " *** STAGE 5 : CONVERT LIPID FILES TO AMBER LIPID 11 FORMAT *** "
echo " "
echo " Currently supported lipids: DPPC DPPE DPPS DPPG DPPA"
echo "                             DOPC DOPE DOPS DOPG DOPA"
echo "                             POPC POPE POPS POPG POPA"
echo "                             DLPC DLPE DLPS DLPG DLPA"
echo "                             DMPC DMPE DMPS DMPG DMPA"
echo "                             CHL1"
echo " "

#The next stage is to take the individual residules making up the various lipid files and rearrange
#things / modify the lipid atom name to give what AMBER lipid11 expects.

touch ./lipid_tmp_final

#DPPC
if ( $dppc_found == 1 ) then
  #DPPC is 130 atoms per lipid. 
  set atoms_per_lipid = 130

  echo " Processing: DPPC"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DPPC_tmp | awk '{print$1}'`
  echo "  DPPC Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DPPC Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DPPC Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DPPC_tmp
  rm -f ../DPPC_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DPPC Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DPPC" $i > tail1
    grep -A42 "C33 DPPC" $i >> tail1

    grep -A2 "C22 DPPC" $i > tail2
    grep -A42 "C23 DPPC" $i >> tail2

    head -n32 $i > head
    grep -A5 "C3  DPPC" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DPPC/ C12 PA  /' tail1 >tail1_tmp
    sed 's/ C33 DPPC/ C13 PA  /' tail1_tmp >tail1
    sed 's/ C34 DPPC/ C14 PA  /' tail1 >tail1_tmp
    sed 's/ C35 DPPC/ C15 PA  /' tail1_tmp >tail1
    sed 's/ C36 DPPC/ C16 PA  /' tail1 >tail1_tmp
    sed 's/ C37 DPPC/ C17 PA  /' tail1_tmp >tail1
    sed 's/ C38 DPPC/ C18 PA  /' tail1 >tail1_tmp
    sed 's/ C39 DPPC/ C19 PA  /' tail1_tmp >tail1
    sed 's/C310 DPPC/C110 PA  /' tail1 >tail1_tmp
    sed 's/C311 DPPC/C111 PA  /' tail1_tmp >tail1
    sed 's/C312 DPPC/C112 PA  /' tail1 >tail1_tmp
    sed 's/C313 DPPC/C113 PA  /' tail1_tmp >tail1
    sed 's/C314 DPPC/C114 PA  /' tail1 >tail1_tmp
    sed 's/C315 DPPC/C115 PA  /' tail1_tmp >tail1
    sed 's/C316 DPPC/C116 PA  /' tail1 >tail1_tmp
    sed 's/ H2X DPPC/ H2R PA  /' tail1_tmp >tail1
    sed 's/ H2Y DPPC/ H2S PA  /' tail1 >tail1_tmp
    sed 's/ H3X DPPC/ H3R PA  /' tail1_tmp >tail1
    sed 's/ H3Y DPPC/ H3S PA  /' tail1 >tail1_tmp
    sed 's/ H4X DPPC/ H4R PA  /' tail1_tmp >tail1
    sed 's/ H4Y DPPC/ H4S PA  /' tail1 >tail1_tmp
    sed 's/ H5X DPPC/ H5R PA  /' tail1_tmp >tail1
    sed 's/ H5Y DPPC/ H5S PA  /' tail1 >tail1_tmp
    sed 's/ H6X DPPC/ H6R PA  /' tail1_tmp >tail1
    sed 's/ H6Y DPPC/ H6S PA  /' tail1 >tail1_tmp
    sed 's/ H7X DPPC/ H7R PA  /' tail1_tmp >tail1
    sed 's/ H7Y DPPC/ H7S PA  /' tail1 >tail1_tmp
    sed 's/ H8X DPPC/ H8R PA  /' tail1_tmp >tail1
    sed 's/ H8Y DPPC/ H8S PA  /' tail1 >tail1_tmp
    sed 's/ H9X DPPC/ H9R PA  /' tail1_tmp >tail1
    sed 's/ H9Y DPPC/ H9S PA  /' tail1 >tail1_tmp
    sed 's/H10X DPPC/H10R PA  /' tail1_tmp >tail1
    sed 's/H10Y DPPC/H10S PA  /' tail1 >tail1_tmp
    sed 's/H11X DPPC/H11R PA  /' tail1_tmp >tail1
    sed 's/H11Y DPPC/H11S PA  /' tail1 >tail1_tmp
    sed 's/H12X DPPC/H12R PA  /' tail1_tmp >tail1
    sed 's/H12Y DPPC/H12S PA  /' tail1 >tail1_tmp
    sed 's/H13X DPPC/H13R PA  /' tail1_tmp >tail1
    sed 's/H13Y DPPC/H13S PA  /' tail1 >tail1_tmp
    sed 's/H14X DPPC/H14R PA  /' tail1_tmp >tail1
    sed 's/H14Y DPPC/H14S PA  /' tail1 >tail1_tmp
    sed 's/H15X DPPC/H15R PA  /' tail1_tmp >tail1
    sed 's/H15Y DPPC/H15S PA  /' tail1 >tail1_tmp
    sed 's/H16X DPPC/H16R PA  /' tail1_tmp >tail1
    sed 's/H16Y DPPC/H16S PA  /' tail1 >tail1_tmp
    sed 's/H16Z DPPC/H16T PA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DPPC/ C12 PA  /' tail2 >tail2_tmp
    sed 's/ C23 DPPC/ C13 PA  /' tail2_tmp >tail2
    sed 's/ C24 DPPC/ C14 PA  /' tail2 >tail2_tmp
    sed 's/ C25 DPPC/ C15 PA  /' tail2_tmp >tail2
    sed 's/ C26 DPPC/ C16 PA  /' tail2 >tail2_tmp
    sed 's/ C27 DPPC/ C17 PA  /' tail2_tmp >tail2
    sed 's/ C28 DPPC/ C18 PA  /' tail2 >tail2_tmp
    sed 's/ C29 DPPC/ C19 PA  /' tail2_tmp >tail2
    sed 's/C210 DPPC/C110 PA  /' tail2 >tail2_tmp
    sed 's/C211 DPPC/C111 PA  /' tail2_tmp >tail2
    sed 's/C212 DPPC/C112 PA  /' tail2 >tail2_tmp
    sed 's/C213 DPPC/C113 PA  /' tail2_tmp >tail2
    sed 's/C214 DPPC/C114 PA  /' tail2 >tail2_tmp
    sed 's/C215 DPPC/C115 PA  /' tail2_tmp >tail2
    sed 's/C216 DPPC/C116 PA  /' tail2 >tail2_tmp
    sed 's/ H2R DPPC/ H2R PA  /' tail2_tmp >tail2
    sed 's/ H2S DPPC/ H2S PA  /' tail2 >tail2_tmp
    sed 's/ H3R DPPC/ H3R PA  /' tail2_tmp >tail2
    sed 's/ H3S DPPC/ H3S PA  /' tail2 >tail2_tmp
    sed 's/ H4R DPPC/ H4R PA  /' tail2_tmp >tail2
    sed 's/ H4S DPPC/ H4S PA  /' tail2 >tail2_tmp
    sed 's/ H5R DPPC/ H5R PA  /' tail2_tmp >tail2
    sed 's/ H5S DPPC/ H5S PA  /' tail2 >tail2_tmp
    sed 's/ H6R DPPC/ H6R PA  /' tail2_tmp >tail2
    sed 's/ H6S DPPC/ H6S PA  /' tail2 >tail2_tmp
    sed 's/ H7R DPPC/ H7R PA  /' tail2_tmp >tail2
    sed 's/ H7S DPPC/ H7S PA  /' tail2 >tail2_tmp
    sed 's/ H8R DPPC/ H8R PA  /' tail2_tmp >tail2
    sed 's/ H8S DPPC/ H8S PA  /' tail2 >tail2_tmp
    sed 's/ H9R DPPC/ H9R PA  /' tail2_tmp >tail2
    sed 's/ H9S DPPC/ H9S PA  /' tail2 >tail2_tmp
    sed 's/H10R DPPC/H10R PA  /' tail2_tmp >tail2
    sed 's/H10S DPPC/H10S PA  /' tail2 >tail2_tmp
    sed 's/H11R DPPC/H11R PA  /' tail2_tmp >tail2
    sed 's/H11S DPPC/H11S PA  /' tail2 >tail2_tmp
    sed 's/H12R DPPC/H12R PA  /' tail2_tmp >tail2
    sed 's/H12S DPPC/H12S PA  /' tail2 >tail2_tmp
    sed 's/H13R DPPC/H13R PA  /' tail2_tmp >tail2
    sed 's/H13S DPPC/H13S PA  /' tail2 >tail2_tmp
    sed 's/H14R DPPC/H14R PA  /' tail2_tmp >tail2
    sed 's/H14S DPPC/H14S PA  /' tail2 >tail2_tmp
    sed 's/H15R DPPC/H15R PA  /' tail2_tmp >tail2
    sed 's/H15S DPPC/H15S PA  /' tail2 >tail2_tmp
    sed 's/H16R DPPC/H16R PA  /' tail2_tmp >tail2
    sed 's/H16S DPPC/H16S PA  /' tail2 >tail2_tmp
    sed 's/H16T DPPC/H16T PA  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DPPC/ N31 PC  /' head > head_tmp
    sed 's/ C13 DPPC/ C33 PC  /' head_tmp > head
    sed 's/ C14 DPPC/ C34 PC  /' head > head_tmp
    sed 's/ C15 DPPC/ C35 PC  /' head_tmp > head
    sed 's/ C12 DPPC/ C32 PC  /' head > head_tmp
    sed 's/ C11 DPPC/ C31 PC  /' head_tmp > head
    sed 's/ O12 DPPC/ O32 PC  /' head > head_tmp
    sed 's/ P   DPPC/ P31 PC  /' head_tmp > head
    sed 's/ O11 DPPC/ O31 PC  /' head > head_tmp
    sed 's/ O13 DPPC/ O33 PC  /' head_tmp > head
    sed 's/ O14 DPPC/ O34 PC  /' head > head_tmp
    sed 's/ C1  DPPC/ C3  PC  /' head_tmp > head
    sed 's/ C2  DPPC/ C2  PC  /' head > head_tmp
    sed 's/ C3  DPPC/ C1  PC  /' head_tmp > head
    sed 's/ O31 DPPC/ O11 PC  /' head > head_tmp
    sed 's/ C31 DPPC/ C11 PC  /' head_tmp > head
    sed 's/ O32 DPPC/ O12 PC  /' head > head_tmp
    sed 's/ O21 DPPC/ O21 PC  /' head_tmp > head
    sed 's/ C21 DPPC/ C21 PC  /' head > head_tmp
    sed 's/ O22 DPPC/ O22 PC  /' head_tmp > head
    sed 's/ HX  DPPC/ HR  PC  /' head > head_tmp
    sed 's/ HY  DPPC/ HS  PC  /' head_tmp > head
    sed 's/ HS  DPPC/ HX  PC  /' head > head_tmp
    sed 's/ HA  DPPC/ HA  PC  /' head_tmp > head
    sed 's/ HB  DPPC/ HB  PC  /' head > head_tmp
    sed 's/H11A DPPC/ H1A PC  /' head_tmp > head
    sed 's/H11B DPPC/ H1B PC  /' head > head_tmp
    sed 's/H12A DPPC/ H2A PC  /' head_tmp > head
    sed 's/H12B DPPC/ H2B PC  /' head > head_tmp
    sed 's/H13A DPPC/ H3A PC  /' head_tmp > head
    sed 's/H13B DPPC/ H3B PC  /' head > head_tmp
    sed 's/H13C DPPC/ H3C PC  /' head_tmp > head
    sed 's/H14A DPPC/ H4A PC  /' head > head_tmp
    sed 's/H14B DPPC/ H4B PC  /' head_tmp > head
    sed 's/H14C DPPC/ H4C PC  /' head > head_tmp
    sed 's/H15A DPPC/ H5A PC  /' head_tmp > head
    sed 's/H15B DPPC/ H5B PC  /' head > head_tmp
    sed 's/H15C DPPC/ H5C PC  /' head_tmp > head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DPPE
if ( $dppe_found == 1 ) then
  #DPPE is 121 atoms per lipid. 
  set atoms_per_lipid = 121

  echo " Processing: DPPE"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DPPE_tmp | awk '{print$1}'`
  echo "  DPPE Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DPPE Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DPPE Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DPPE_tmp
  rm -f ../DPPE_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DPPE Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DPPE" $i > tail1
    grep -A42 "C33 DPPE" $i >> tail1

    grep -A2 "C22 DPPE" $i > tail2
    grep -A42 "C23 DPPE" $i >> tail2

    head -n23 $i > head
    grep -A5 "C3  DPPE" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DPPE/ C12 PA  /' tail1 >tail1_tmp
    sed 's/ C33 DPPE/ C13 PA  /' tail1_tmp >tail1
    sed 's/ C34 DPPE/ C14 PA  /' tail1 >tail1_tmp
    sed 's/ C35 DPPE/ C15 PA  /' tail1_tmp >tail1
    sed 's/ C36 DPPE/ C16 PA  /' tail1 >tail1_tmp
    sed 's/ C37 DPPE/ C17 PA  /' tail1_tmp >tail1
    sed 's/ C38 DPPE/ C18 PA  /' tail1 >tail1_tmp
    sed 's/ C39 DPPE/ C19 PA  /' tail1_tmp >tail1
    sed 's/C310 DPPE/C110 PA  /' tail1 >tail1_tmp
    sed 's/C311 DPPE/C111 PA  /' tail1_tmp >tail1
    sed 's/C312 DPPE/C112 PA  /' tail1 >tail1_tmp
    sed 's/C313 DPPE/C113 PA  /' tail1_tmp >tail1
    sed 's/C314 DPPE/C114 PA  /' tail1 >tail1_tmp
    sed 's/C315 DPPE/C115 PA  /' tail1_tmp >tail1
    sed 's/C316 DPPE/C116 PA  /' tail1 >tail1_tmp
    sed 's/ H2X DPPE/ H2R PA  /' tail1_tmp >tail1
    sed 's/ H2Y DPPE/ H2S PA  /' tail1 >tail1_tmp
    sed 's/ H3X DPPE/ H3R PA  /' tail1_tmp >tail1
    sed 's/ H3Y DPPE/ H3S PA  /' tail1 >tail1_tmp
    sed 's/ H4X DPPE/ H4R PA  /' tail1_tmp >tail1
    sed 's/ H4Y DPPE/ H4S PA  /' tail1 >tail1_tmp
    sed 's/ H5X DPPE/ H5R PA  /' tail1_tmp >tail1
    sed 's/ H5Y DPPE/ H5S PA  /' tail1 >tail1_tmp
    sed 's/ H6X DPPE/ H6R PA  /' tail1_tmp >tail1
    sed 's/ H6Y DPPE/ H6S PA  /' tail1 >tail1_tmp
    sed 's/ H7X DPPE/ H7R PA  /' tail1_tmp >tail1
    sed 's/ H7Y DPPE/ H7S PA  /' tail1 >tail1_tmp
    sed 's/ H8X DPPE/ H8R PA  /' tail1_tmp >tail1
    sed 's/ H8Y DPPE/ H8S PA  /' tail1 >tail1_tmp
    sed 's/ H9X DPPE/ H9R PA  /' tail1_tmp >tail1
    sed 's/ H9Y DPPE/ H9S PA  /' tail1 >tail1_tmp
    sed 's/H10X DPPE/H10R PA  /' tail1_tmp >tail1
    sed 's/H10Y DPPE/H10S PA  /' tail1 >tail1_tmp
    sed 's/H11X DPPE/H11R PA  /' tail1_tmp >tail1
    sed 's/H11Y DPPE/H11S PA  /' tail1 >tail1_tmp
    sed 's/H12X DPPE/H12R PA  /' tail1_tmp >tail1
    sed 's/H12Y DPPE/H12S PA  /' tail1 >tail1_tmp
    sed 's/H13X DPPE/H13R PA  /' tail1_tmp >tail1
    sed 's/H13Y DPPE/H13S PA  /' tail1 >tail1_tmp
    sed 's/H14X DPPE/H14R PA  /' tail1_tmp >tail1
    sed 's/H14Y DPPE/H14S PA  /' tail1 >tail1_tmp
    sed 's/H15X DPPE/H15R PA  /' tail1_tmp >tail1
    sed 's/H15Y DPPE/H15S PA  /' tail1 >tail1_tmp
    sed 's/H16X DPPE/H16R PA  /' tail1_tmp >tail1
    sed 's/H16Y DPPE/H16S PA  /' tail1 >tail1_tmp
    sed 's/H16Z DPPE/H16T PA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DPPE/ C12 PA  /' tail2 >tail2_tmp
    sed 's/ C23 DPPE/ C13 PA  /' tail2_tmp >tail2
    sed 's/ C24 DPPE/ C14 PA  /' tail2 >tail2_tmp
    sed 's/ C25 DPPE/ C15 PA  /' tail2_tmp >tail2
    sed 's/ C26 DPPE/ C16 PA  /' tail2 >tail2_tmp
    sed 's/ C27 DPPE/ C17 PA  /' tail2_tmp >tail2
    sed 's/ C28 DPPE/ C18 PA  /' tail2 >tail2_tmp
    sed 's/ C29 DPPE/ C19 PA  /' tail2_tmp >tail2
    sed 's/C210 DPPE/C110 PA  /' tail2 >tail2_tmp
    sed 's/C211 DPPE/C111 PA  /' tail2_tmp >tail2
    sed 's/C212 DPPE/C112 PA  /' tail2 >tail2_tmp
    sed 's/C213 DPPE/C113 PA  /' tail2_tmp >tail2
    sed 's/C214 DPPE/C114 PA  /' tail2 >tail2_tmp
    sed 's/C215 DPPE/C115 PA  /' tail2_tmp >tail2
    sed 's/C216 DPPE/C116 PA  /' tail2 >tail2_tmp
    sed 's/ H2R DPPE/ H2R PA  /' tail2_tmp >tail2
    sed 's/ H2S DPPE/ H2S PA  /' tail2 >tail2_tmp
    sed 's/ H3R DPPE/ H3R PA  /' tail2_tmp >tail2
    sed 's/ H3S DPPE/ H3S PA  /' tail2 >tail2_tmp
    sed 's/ H4R DPPE/ H4R PA  /' tail2_tmp >tail2
    sed 's/ H4S DPPE/ H4S PA  /' tail2 >tail2_tmp
    sed 's/ H5R DPPE/ H5R PA  /' tail2_tmp >tail2
    sed 's/ H5S DPPE/ H5S PA  /' tail2 >tail2_tmp
    sed 's/ H6R DPPE/ H6R PA  /' tail2_tmp >tail2
    sed 's/ H6S DPPE/ H6S PA  /' tail2 >tail2_tmp
    sed 's/ H7R DPPE/ H7R PA  /' tail2_tmp >tail2
    sed 's/ H7S DPPE/ H7S PA  /' tail2 >tail2_tmp
    sed 's/ H8R DPPE/ H8R PA  /' tail2_tmp >tail2
    sed 's/ H8S DPPE/ H8S PA  /' tail2 >tail2_tmp
    sed 's/ H9R DPPE/ H9R PA  /' tail2_tmp >tail2
    sed 's/ H9S DPPE/ H9S PA  /' tail2 >tail2_tmp
    sed 's/H10R DPPE/H10R PA  /' tail2_tmp >tail2
    sed 's/H10S DPPE/H10S PA  /' tail2 >tail2_tmp
    sed 's/H11R DPPE/H11R PA  /' tail2_tmp >tail2
    sed 's/H11S DPPE/H11S PA  /' tail2 >tail2_tmp
    sed 's/H12R DPPE/H12R PA  /' tail2_tmp >tail2
    sed 's/H12S DPPE/H12S PA  /' tail2 >tail2_tmp
    sed 's/H13R DPPE/H13R PA  /' tail2_tmp >tail2
    sed 's/H13S DPPE/H13S PA  /' tail2 >tail2_tmp
    sed 's/H14R DPPE/H14R PA  /' tail2_tmp >tail2
    sed 's/H14S DPPE/H14S PA  /' tail2 >tail2_tmp
    sed 's/H15R DPPE/H15R PA  /' tail2_tmp >tail2
    sed 's/H15S DPPE/H15S PA  /' tail2 >tail2_tmp
    sed 's/H16R DPPE/H16R PA  /' tail2_tmp >tail2
    sed 's/H16S DPPE/H16S PA  /' tail2 >tail2_tmp
    sed 's/H16T DPPE/H16T PA  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DPPE/ N31 PE  /' head > head_tmp
    sed 's/ C12 DPPE/ C32 PE  /' head_tmp > head
    sed 's/ C11 DPPE/ C31 PE  /' head > head_tmp
    sed 's/ O12 DPPE/ O32 PE  /' head_tmp > head
    sed 's/ P   DPPE/ P31 PE  /' head > head_tmp
    sed 's/ O11 DPPE/ O31 PE  /' head_tmp > head
    sed 's/ O13 DPPE/ O33 PE  /' head > head_tmp
    sed 's/ O14 DPPE/ O34 PE  /' head_tmp > head
    sed 's/ C1  DPPE/ C3  PE  /' head > head_tmp
    sed 's/ C2  DPPE/ C2  PE  /' head_tmp > head
    sed 's/ C3  DPPE/ C1  PE  /' head > head_tmp
    sed 's/ O31 DPPE/ O11 PE  /' head_tmp > head
    sed 's/ C31 DPPE/ C11 PE  /' head > head_tmp
    sed 's/ O32 DPPE/ O12 PE  /' head_tmp > head
    sed 's/ O21 DPPE/ O21 PE  /' head > head_tmp
    sed 's/ C21 DPPE/ C21 PE  /' head_tmp > head
    sed 's/ O22 DPPE/ O22 PE  /' head > head_tmp
    sed 's/ HX  DPPE/ HR  PE  /' head_tmp > head
    sed 's/ HY  DPPE/ HS  PE  /' head > head_tmp
    sed 's/ HS  DPPE/ HX  PE  /' head_tmp > head
    sed 's/ HA  DPPE/ HA  PE  /' head > head_tmp
    sed 's/ HB  DPPE/ HB  PE  /' head_tmp > head
    sed 's/H11A DPPE/ H1A PE  /' head > head_tmp
    sed 's/H11B DPPE/ H1B PE  /' head_tmp > head
    sed 's/H12A DPPE/ H2A PE  /' head > head_tmp
    sed 's/H12B DPPE/ H2B PE  /' head_tmp > head
    sed 's/ HN1 DPPE/HN1A PE  /' head > head_tmp
    sed 's/ HN2 DPPE/HN1B PE  /' head_tmp > head
    sed 's/ HN3 DPPE/HN1C PE  /' head > head_tmp
    mv head_tmp head  
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DPPS
if ( $dpps_found == 1 ) then
  #DPPS is 123 atoms per lipid. 
  set atoms_per_lipid = 123

  echo " Processing: DPPS"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DPPS_tmp | awk '{print$1}'`
  echo "  DPPS Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DPPS Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DPPS Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DPPS_tmp
  rm -f ../DPPS_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DPPS Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DPPS" $i > tail1
    grep -A42 "C33 DPPS" $i >> tail1

    grep -A2 "C22 DPPS" $i > tail2
    grep -A42 "C23 DPPS" $i >> tail2

    head -n25 $i > head
    grep -A5 "C3  DPPS" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DPPS/ C12 PA  /' tail1 >tail1_tmp
    sed 's/ C33 DPPS/ C13 PA  /' tail1_tmp >tail1
    sed 's/ C34 DPPS/ C14 PA  /' tail1 >tail1_tmp
    sed 's/ C35 DPPS/ C15 PA  /' tail1_tmp >tail1
    sed 's/ C36 DPPS/ C16 PA  /' tail1 >tail1_tmp
    sed 's/ C37 DPPS/ C17 PA  /' tail1_tmp >tail1
    sed 's/ C38 DPPS/ C18 PA  /' tail1 >tail1_tmp
    sed 's/ C39 DPPS/ C19 PA  /' tail1_tmp >tail1
    sed 's/C310 DPPS/C110 PA  /' tail1 >tail1_tmp
    sed 's/C311 DPPS/C111 PA  /' tail1_tmp >tail1
    sed 's/C312 DPPS/C112 PA  /' tail1 >tail1_tmp
    sed 's/C313 DPPS/C113 PA  /' tail1_tmp >tail1
    sed 's/C314 DPPS/C114 PA  /' tail1 >tail1_tmp
    sed 's/C315 DPPS/C115 PA  /' tail1_tmp >tail1
    sed 's/C316 DPPS/C116 PA  /' tail1 >tail1_tmp
    sed 's/ H2X DPPS/ H2R PA  /' tail1_tmp >tail1
    sed 's/ H2Y DPPS/ H2S PA  /' tail1 >tail1_tmp
    sed 's/ H3X DPPS/ H3R PA  /' tail1_tmp >tail1
    sed 's/ H3Y DPPS/ H3S PA  /' tail1 >tail1_tmp
    sed 's/ H4X DPPS/ H4R PA  /' tail1_tmp >tail1
    sed 's/ H4Y DPPS/ H4S PA  /' tail1 >tail1_tmp
    sed 's/ H5X DPPS/ H5R PA  /' tail1_tmp >tail1
    sed 's/ H5Y DPPS/ H5S PA  /' tail1 >tail1_tmp
    sed 's/ H6X DPPS/ H6R PA  /' tail1_tmp >tail1
    sed 's/ H6Y DPPS/ H6S PA  /' tail1 >tail1_tmp
    sed 's/ H7X DPPS/ H7R PA  /' tail1_tmp >tail1
    sed 's/ H7Y DPPS/ H7S PA  /' tail1 >tail1_tmp
    sed 's/ H8X DPPS/ H8R PA  /' tail1_tmp >tail1
    sed 's/ H8Y DPPS/ H8S PA  /' tail1 >tail1_tmp
    sed 's/ H9X DPPS/ H9R PA  /' tail1_tmp >tail1
    sed 's/ H9Y DPPS/ H9S PA  /' tail1 >tail1_tmp
    sed 's/H10X DPPS/H10R PA  /' tail1_tmp >tail1
    sed 's/H10Y DPPS/H10S PA  /' tail1 >tail1_tmp
    sed 's/H11X DPPS/H11R PA  /' tail1_tmp >tail1
    sed 's/H11Y DPPS/H11S PA  /' tail1 >tail1_tmp
    sed 's/H12X DPPS/H12R PA  /' tail1_tmp >tail1
    sed 's/H12Y DPPS/H12S PA  /' tail1 >tail1_tmp
    sed 's/H13X DPPS/H13R PA  /' tail1_tmp >tail1
    sed 's/H13Y DPPS/H13S PA  /' tail1 >tail1_tmp
    sed 's/H14X DPPS/H14R PA  /' tail1_tmp >tail1
    sed 's/H14Y DPPS/H14S PA  /' tail1 >tail1_tmp
    sed 's/H15X DPPS/H15R PA  /' tail1_tmp >tail1
    sed 's/H15Y DPPS/H15S PA  /' tail1 >tail1_tmp
    sed 's/H16X DPPS/H16R PA  /' tail1_tmp >tail1
    sed 's/H16Y DPPS/H16S PA  /' tail1 >tail1_tmp
    sed 's/H16Z DPPS/H16T PA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DPPS/ C12 PA  /' tail2 >tail2_tmp
    sed 's/ C23 DPPS/ C13 PA  /' tail2_tmp >tail2
    sed 's/ C24 DPPS/ C14 PA  /' tail2 >tail2_tmp
    sed 's/ C25 DPPS/ C15 PA  /' tail2_tmp >tail2
    sed 's/ C26 DPPS/ C16 PA  /' tail2 >tail2_tmp
    sed 's/ C27 DPPS/ C17 PA  /' tail2_tmp >tail2
    sed 's/ C28 DPPS/ C18 PA  /' tail2 >tail2_tmp
    sed 's/ C29 DPPS/ C19 PA  /' tail2_tmp >tail2
    sed 's/C210 DPPS/C110 PA  /' tail2 >tail2_tmp
    sed 's/C211 DPPS/C111 PA  /' tail2_tmp >tail2
    sed 's/C212 DPPS/C112 PA  /' tail2 >tail2_tmp
    sed 's/C213 DPPS/C113 PA  /' tail2_tmp >tail2
    sed 's/C214 DPPS/C114 PA  /' tail2 >tail2_tmp
    sed 's/C215 DPPS/C115 PA  /' tail2_tmp >tail2
    sed 's/C216 DPPS/C116 PA  /' tail2 >tail2_tmp
    sed 's/ H2R DPPS/ H2R PA  /' tail2_tmp >tail2
    sed 's/ H2S DPPS/ H2S PA  /' tail2 >tail2_tmp
    sed 's/ H3R DPPS/ H3R PA  /' tail2_tmp >tail2
    sed 's/ H3S DPPS/ H3S PA  /' tail2 >tail2_tmp
    sed 's/ H4R DPPS/ H4R PA  /' tail2_tmp >tail2
    sed 's/ H4S DPPS/ H4S PA  /' tail2 >tail2_tmp
    sed 's/ H5R DPPS/ H5R PA  /' tail2_tmp >tail2
    sed 's/ H5S DPPS/ H5S PA  /' tail2 >tail2_tmp
    sed 's/ H6R DPPS/ H6R PA  /' tail2_tmp >tail2
    sed 's/ H6S DPPS/ H6S PA  /' tail2 >tail2_tmp
    sed 's/ H7R DPPS/ H7R PA  /' tail2_tmp >tail2
    sed 's/ H7S DPPS/ H7S PA  /' tail2 >tail2_tmp
    sed 's/ H8R DPPS/ H8R PA  /' tail2_tmp >tail2
    sed 's/ H8S DPPS/ H8S PA  /' tail2 >tail2_tmp
    sed 's/ H9R DPPS/ H9R PA  /' tail2_tmp >tail2
    sed 's/ H9S DPPS/ H9S PA  /' tail2 >tail2_tmp
    sed 's/H10R DPPS/H10R PA  /' tail2_tmp >tail2
    sed 's/H10S DPPS/H10S PA  /' tail2 >tail2_tmp
    sed 's/H11R DPPS/H11R PA  /' tail2_tmp >tail2
    sed 's/H11S DPPS/H11S PA  /' tail2 >tail2_tmp
    sed 's/H12R DPPS/H12R PA  /' tail2_tmp >tail2
    sed 's/H12S DPPS/H12S PA  /' tail2 >tail2_tmp
    sed 's/H13R DPPS/H13R PA  /' tail2_tmp >tail2
    sed 's/H13S DPPS/H13S PA  /' tail2 >tail2_tmp
    sed 's/H14R DPPS/H14R PA  /' tail2_tmp >tail2
    sed 's/H14S DPPS/H14S PA  /' tail2 >tail2_tmp
    sed 's/H15R DPPS/H15R PA  /' tail2_tmp >tail2
    sed 's/H15S DPPS/H15S PA  /' tail2 >tail2_tmp
    sed 's/H16R DPPS/H16R PA  /' tail2_tmp >tail2
    sed 's/H16S DPPS/H16S PA  /' tail2 >tail2_tmp
    sed 's/H16T DPPS/H16T PA  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DPPS/ N31 PS  /' head > head_tmp
    sed 's/ C12 DPPS/ C32 PS  /' head_tmp > head
    sed 's/ C13 DPPS/ C33 PS  /' head > head_tmp
    sed 's/O13A DPPS/ O35 PS  /' head_tmp > head
    sed 's/O13B DPPS/ O36 PS  /' head > head_tmp
    sed 's/ O12 DPPS/ O32 PS  /' head_tmp > head
    sed 's/ P   DPPS/ P31 PS  /' head > head_tmp
    sed 's/ O11 DPPS/ O31 PS  /' head_tmp > head
    sed 's/ O13 DPPS/ O33 PS  /' head > head_tmp
    sed 's/ C11 DPPS/ C31 PS  /' head_tmp > head
    sed 's/ O14 DPPS/ O34 PS  /' head > head_tmp
    sed 's/ C1  DPPS/ C3  PS  /' head_tmp > head
    sed 's/ C2  DPPS/ C2  PS  /' head > head_tmp
    sed 's/ C3  DPPS/ C1  PS  /' head_tmp > head
    sed 's/ O31 DPPS/ O11 PS  /' head > head_tmp
    sed 's/ C31 DPPS/ C11 PS  /' head_tmp > head
    sed 's/ O32 DPPS/ O12 PS  /' head > head_tmp
    sed 's/ O21 DPPS/ O21 PS  /' head_tmp > head
    sed 's/ C21 DPPS/ C21 PS  /' head > head_tmp
    sed 's/ O22 DPPS/ O22 PS  /' head_tmp > head
    sed 's/ HX  DPPS/ HR  PS  /' head > head_tmp
    sed 's/ HY  DPPS/ HS  PS  /' head_tmp > head
    sed 's/ HS  DPPS/ HX  PS  /' head > head_tmp
    sed 's/ HA  DPPS/ HA  PS  /' head_tmp > head
    sed 's/ HB  DPPS/ HB  PS  /' head > head_tmp
    sed 's/H11A DPPS/ H1A PS  /' head_tmp > head
    sed 's/H11B DPPS/ H1B PS  /' head > head_tmp
    sed 's/H12A DPPS/ H2A PS  /' head_tmp > head
    sed 's/ HN1 DPPS/HN1A PS  /' head > head_tmp
    sed 's/ HN2 DPPS/HN1B PS  /' head_tmp > head
    sed 's/ HN3 DPPS/HN1C PS  /' head > head_tmp
    mv head_tmp head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DPPG
if ( $dppg_found == 1 ) then
  #DPPG is 123 atoms per lipid. 
  set atoms_per_lipid = 123

  echo " Processing: DPPG"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DPPG_tmp | awk '{print$1}'`
  echo "  DPPG Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DPPG Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DPPG Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DPPG_tmp
  rm -f ../DPPG_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DPPG Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DPPG" $i > tail1
    grep -A42 "C33 DPPG" $i >> tail1

    grep -A2 "C22 DPPG" $i > tail2
    grep -A42 "C23 DPPG" $i >> tail2

    head -n25 $i > head
    grep -A5 "C3  DPPG" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DPPG/ C12 PA  /' tail1 >tail1_tmp
    sed 's/ C33 DPPG/ C13 PA  /' tail1_tmp >tail1
    sed 's/ C34 DPPG/ C14 PA  /' tail1 >tail1_tmp
    sed 's/ C35 DPPG/ C15 PA  /' tail1_tmp >tail1
    sed 's/ C36 DPPG/ C16 PA  /' tail1 >tail1_tmp
    sed 's/ C37 DPPG/ C17 PA  /' tail1_tmp >tail1
    sed 's/ C38 DPPG/ C18 PA  /' tail1 >tail1_tmp
    sed 's/ C39 DPPG/ C19 PA  /' tail1_tmp >tail1
    sed 's/C310 DPPG/C110 PA  /' tail1 >tail1_tmp
    sed 's/C311 DPPG/C111 PA  /' tail1_tmp >tail1
    sed 's/C312 DPPG/C112 PA  /' tail1 >tail1_tmp
    sed 's/C313 DPPG/C113 PA  /' tail1_tmp >tail1
    sed 's/C314 DPPG/C114 PA  /' tail1 >tail1_tmp
    sed 's/C315 DPPG/C115 PA  /' tail1_tmp >tail1
    sed 's/C316 DPPG/C116 PA  /' tail1 >tail1_tmp
    sed 's/ H2X DPPG/ H2R PA  /' tail1_tmp >tail1
    sed 's/ H2Y DPPG/ H2S PA  /' tail1 >tail1_tmp
    sed 's/ H3X DPPG/ H3R PA  /' tail1_tmp >tail1
    sed 's/ H3Y DPPG/ H3S PA  /' tail1 >tail1_tmp
    sed 's/ H4X DPPG/ H4R PA  /' tail1_tmp >tail1
    sed 's/ H4Y DPPG/ H4S PA  /' tail1 >tail1_tmp
    sed 's/ H5X DPPG/ H5R PA  /' tail1_tmp >tail1
    sed 's/ H5Y DPPG/ H5S PA  /' tail1 >tail1_tmp
    sed 's/ H6X DPPG/ H6R PA  /' tail1_tmp >tail1
    sed 's/ H6Y DPPG/ H6S PA  /' tail1 >tail1_tmp
    sed 's/ H7X DPPG/ H7R PA  /' tail1_tmp >tail1
    sed 's/ H7Y DPPG/ H7S PA  /' tail1 >tail1_tmp
    sed 's/ H8X DPPG/ H8R PA  /' tail1_tmp >tail1
    sed 's/ H8Y DPPG/ H8S PA  /' tail1 >tail1_tmp
    sed 's/ H9X DPPG/ H9R PA  /' tail1_tmp >tail1
    sed 's/ H9Y DPPG/ H9S PA  /' tail1 >tail1_tmp
    sed 's/H10X DPPG/H10R PA  /' tail1_tmp >tail1
    sed 's/H10Y DPPG/H10S PA  /' tail1 >tail1_tmp
    sed 's/H11X DPPG/H11R PA  /' tail1_tmp >tail1
    sed 's/H11Y DPPG/H11S PA  /' tail1 >tail1_tmp
    sed 's/H12X DPPG/H12R PA  /' tail1_tmp >tail1
    sed 's/H12Y DPPG/H12S PA  /' tail1 >tail1_tmp
    sed 's/H13X DPPG/H13R PA  /' tail1_tmp >tail1
    sed 's/H13Y DPPG/H13S PA  /' tail1 >tail1_tmp
    sed 's/H14X DPPG/H14R PA  /' tail1_tmp >tail1
    sed 's/H14Y DPPG/H14S PA  /' tail1 >tail1_tmp
    sed 's/H15X DPPG/H15R PA  /' tail1_tmp >tail1
    sed 's/H15Y DPPG/H15S PA  /' tail1 >tail1_tmp
    sed 's/H16X DPPG/H16R PA  /' tail1_tmp >tail1
    sed 's/H16Y DPPG/H16S PA  /' tail1 >tail1_tmp
    sed 's/H16Z DPPG/H16T PA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DPPG/ C12 PA  /' tail2 >tail2_tmp
    sed 's/ C23 DPPG/ C13 PA  /' tail2_tmp >tail2
    sed 's/ C24 DPPG/ C14 PA  /' tail2 >tail2_tmp
    sed 's/ C25 DPPG/ C15 PA  /' tail2_tmp >tail2
    sed 's/ C26 DPPG/ C16 PA  /' tail2 >tail2_tmp
    sed 's/ C27 DPPG/ C17 PA  /' tail2_tmp >tail2
    sed 's/ C28 DPPG/ C18 PA  /' tail2 >tail2_tmp
    sed 's/ C29 DPPG/ C19 PA  /' tail2_tmp >tail2
    sed 's/C210 DPPG/C110 PA  /' tail2 >tail2_tmp
    sed 's/C211 DPPG/C111 PA  /' tail2_tmp >tail2
    sed 's/C212 DPPG/C112 PA  /' tail2 >tail2_tmp
    sed 's/C213 DPPG/C113 PA  /' tail2_tmp >tail2
    sed 's/C214 DPPG/C114 PA  /' tail2 >tail2_tmp
    sed 's/C215 DPPG/C115 PA  /' tail2_tmp >tail2
    sed 's/C216 DPPG/C116 PA  /' tail2 >tail2_tmp
    sed 's/ H2R DPPG/ H2R PA  /' tail2_tmp >tail2
    sed 's/ H2S DPPG/ H2S PA  /' tail2 >tail2_tmp
    sed 's/ H3R DPPG/ H3R PA  /' tail2_tmp >tail2
    sed 's/ H3S DPPG/ H3S PA  /' tail2 >tail2_tmp
    sed 's/ H4R DPPG/ H4R PA  /' tail2_tmp >tail2
    sed 's/ H4S DPPG/ H4S PA  /' tail2 >tail2_tmp
    sed 's/ H5R DPPG/ H5R PA  /' tail2_tmp >tail2
    sed 's/ H5S DPPG/ H5S PA  /' tail2 >tail2_tmp
    sed 's/ H6R DPPG/ H6R PA  /' tail2_tmp >tail2
    sed 's/ H6S DPPG/ H6S PA  /' tail2 >tail2_tmp
    sed 's/ H7R DPPG/ H7R PA  /' tail2_tmp >tail2
    sed 's/ H7S DPPG/ H7S PA  /' tail2 >tail2_tmp
    sed 's/ H8R DPPG/ H8R PA  /' tail2_tmp >tail2
    sed 's/ H8S DPPG/ H8S PA  /' tail2 >tail2_tmp
    sed 's/ H9R DPPG/ H9R PA  /' tail2_tmp >tail2
    sed 's/ H9S DPPG/ H9S PA  /' tail2 >tail2_tmp
    sed 's/H10R DPPG/H10R PA  /' tail2_tmp >tail2
    sed 's/H10S DPPG/H10S PA  /' tail2 >tail2_tmp
    sed 's/H11R DPPG/H11R PA  /' tail2_tmp >tail2
    sed 's/H11S DPPG/H11S PA  /' tail2 >tail2_tmp
    sed 's/H12R DPPG/H12R PA  /' tail2_tmp >tail2
    sed 's/H12S DPPG/H12S PA  /' tail2 >tail2_tmp
    sed 's/H13R DPPG/H13R PA  /' tail2_tmp >tail2
    sed 's/H13S DPPG/H13S PA  /' tail2 >tail2_tmp
    sed 's/H14R DPPG/H14R PA  /' tail2_tmp >tail2
    sed 's/H14S DPPG/H14S PA  /' tail2 >tail2_tmp
    sed 's/H15R DPPG/H15R PA  /' tail2_tmp >tail2
    sed 's/H15S DPPG/H15S PA  /' tail2 >tail2_tmp
    sed 's/H16R DPPG/H16R PA  /' tail2_tmp >tail2
    sed 's/H16S DPPG/H16S PA  /' tail2 >tail2_tmp
    sed 's/H16T DPPG/H16T PA  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ C13 DPPG/ C33 PGR /' head > head_tmp
    sed 's/H13A DPPG/ H3A PGR /' head_tmp > head
    sed 's/H13B DPPG/ H3B PGR /' head > head_tmp
    sed 's/ OC3 DPPG/ O36 PGR /' head_tmp > head
    sed 's/ HO3 DPPG/HO6A PGR /' head > head_tmp
    sed 's/ C12 DPPG/ C32 PGR /' head_tmp > head
    sed 's/H12A DPPG/ H2A PGR /' head > head_tmp
    sed 's/ OC2 DPPG/ O35 PGR /' head_tmp > head
    sed 's/ HO2 DPPG/HO5A PGR /' head > head_tmp
    sed 's/ C11 DPPG/ C31 PGR /' head_tmp > head
    sed 's/H11A DPPG/ H1A PGR /' head > head_tmp
    sed 's/H11B DPPG/ H1B PGR /' head_tmp > head
    sed 's/ P   DPPG/ P31 PGR /' head > head_tmp
    sed 's/ O13 DPPG/ O33 PGR /' head_tmp > head
    sed 's/ O14 DPPG/ O34 PGR /' head > head_tmp
    sed 's/ O12 DPPG/ O32 PGR /' head_tmp > head
    sed 's/ O11 DPPG/ O31 PGR /' head > head_tmp
    sed 's/ C1  DPPG/ C3  PGR /' head_tmp > head
    sed 's/ HA  DPPG/ HA  PGR /' head > head_tmp
    sed 's/ HB  DPPG/ HB  PGR /' head_tmp > head
    sed 's/ C2  DPPG/ C2  PGR /' head > head_tmp
    sed 's/ HS  DPPG/ HX  PGR /' head_tmp > head
    sed 's/ O21 DPPG/ O21 PGR /' head > head_tmp
    sed 's/ C21 DPPG/ C21 PGR /' head_tmp > head
    sed 's/ O22 DPPG/ O22 PGR /' head > head_tmp
    sed 's/ C3  DPPG/ C1  PGR /' head_tmp > head
    sed 's/ HX  DPPG/ HR  PGR /' head > head_tmp
    sed 's/ HY  DPPG/ HS  PGR /' head_tmp > head
    sed 's/ O31 DPPG/ O11 PGR /' head > head_tmp
    sed 's/ C31 DPPG/ C11 PGR /' head_tmp > head
    sed 's/ O32 DPPG/ O12 PGR /' head > head_tmp
    mv head_tmp head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DPPA
if ( $dppa_found == 1 ) then
  #DPPA is 112 atoms per lipid. 
  set atoms_per_lipid = 112

  echo " Processing: DPPA"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DPPA_tmp | awk '{print$1}'`
  echo "  DPPA Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DPPA Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DPPA Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DPPA_tmp
  rm -f ../DPPA_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DPPA Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DPPA" $i > tail1
    grep -A42 "C33 DPPA" $i >> tail1

    grep -A2 "C22 DPPA" $i > tail2
    grep -A42 "C23 DPPA" $i >> tail2

    head -n14 $i > head
    grep -A5 "C3  DPPA" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DPPA/ C12 PA  /' tail1 >tail1_tmp
    sed 's/ C33 DPPA/ C13 PA  /' tail1_tmp >tail1
    sed 's/ C34 DPPA/ C14 PA  /' tail1 >tail1_tmp
    sed 's/ C35 DPPA/ C15 PA  /' tail1_tmp >tail1
    sed 's/ C36 DPPA/ C16 PA  /' tail1 >tail1_tmp
    sed 's/ C37 DPPA/ C17 PA  /' tail1_tmp >tail1
    sed 's/ C38 DPPA/ C18 PA  /' tail1 >tail1_tmp
    sed 's/ C39 DPPA/ C19 PA  /' tail1_tmp >tail1
    sed 's/C310 DPPA/C110 PA  /' tail1 >tail1_tmp
    sed 's/C311 DPPA/C111 PA  /' tail1_tmp >tail1
    sed 's/C312 DPPA/C112 PA  /' tail1 >tail1_tmp
    sed 's/C313 DPPA/C113 PA  /' tail1_tmp >tail1
    sed 's/C314 DPPA/C114 PA  /' tail1 >tail1_tmp
    sed 's/C315 DPPA/C115 PA  /' tail1_tmp >tail1
    sed 's/C316 DPPA/C116 PA  /' tail1 >tail1_tmp
    sed 's/ H2X DPPA/ H2R PA  /' tail1_tmp >tail1
    sed 's/ H2Y DPPA/ H2S PA  /' tail1 >tail1_tmp
    sed 's/ H3X DPPA/ H3R PA  /' tail1_tmp >tail1
    sed 's/ H3Y DPPA/ H3S PA  /' tail1 >tail1_tmp
    sed 's/ H4X DPPA/ H4R PA  /' tail1_tmp >tail1
    sed 's/ H4Y DPPA/ H4S PA  /' tail1 >tail1_tmp
    sed 's/ H5X DPPA/ H5R PA  /' tail1_tmp >tail1
    sed 's/ H5Y DPPA/ H5S PA  /' tail1 >tail1_tmp
    sed 's/ H6X DPPA/ H6R PA  /' tail1_tmp >tail1
    sed 's/ H6Y DPPA/ H6S PA  /' tail1 >tail1_tmp
    sed 's/ H7X DPPA/ H7R PA  /' tail1_tmp >tail1
    sed 's/ H7Y DPPA/ H7S PA  /' tail1 >tail1_tmp
    sed 's/ H8X DPPA/ H8R PA  /' tail1_tmp >tail1
    sed 's/ H8Y DPPA/ H8S PA  /' tail1 >tail1_tmp
    sed 's/ H9X DPPA/ H9R PA  /' tail1_tmp >tail1
    sed 's/ H9Y DPPA/ H9S PA  /' tail1 >tail1_tmp
    sed 's/H10X DPPA/H10R PA  /' tail1_tmp >tail1
    sed 's/H10Y DPPA/H10S PA  /' tail1 >tail1_tmp
    sed 's/H11X DPPA/H11R PA  /' tail1_tmp >tail1
    sed 's/H11Y DPPA/H11S PA  /' tail1 >tail1_tmp
    sed 's/H12X DPPA/H12R PA  /' tail1_tmp >tail1
    sed 's/H12Y DPPA/H12S PA  /' tail1 >tail1_tmp
    sed 's/H13X DPPA/H13R PA  /' tail1_tmp >tail1
    sed 's/H13Y DPPA/H13S PA  /' tail1 >tail1_tmp
    sed 's/H14X DPPA/H14R PA  /' tail1_tmp >tail1
    sed 's/H14Y DPPA/H14S PA  /' tail1 >tail1_tmp
    sed 's/H15X DPPA/H15R PA  /' tail1_tmp >tail1
    sed 's/H15Y DPPA/H15S PA  /' tail1 >tail1_tmp
    sed 's/H16X DPPA/H16R PA  /' tail1_tmp >tail1
    sed 's/H16Y DPPA/H16S PA  /' tail1 >tail1_tmp
    sed 's/H16Z DPPA/H16T PA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DPPA/ C12 PA  /' tail2 >tail2_tmp
    sed 's/ C23 DPPA/ C13 PA  /' tail2_tmp >tail2
    sed 's/ C24 DPPA/ C14 PA  /' tail2 >tail2_tmp
    sed 's/ C25 DPPA/ C15 PA  /' tail2_tmp >tail2
    sed 's/ C26 DPPA/ C16 PA  /' tail2 >tail2_tmp
    sed 's/ C27 DPPA/ C17 PA  /' tail2_tmp >tail2
    sed 's/ C28 DPPA/ C18 PA  /' tail2 >tail2_tmp
    sed 's/ C29 DPPA/ C19 PA  /' tail2_tmp >tail2
    sed 's/C210 DPPA/C110 PA  /' tail2 >tail2_tmp
    sed 's/C211 DPPA/C111 PA  /' tail2_tmp >tail2
    sed 's/C212 DPPA/C112 PA  /' tail2 >tail2_tmp
    sed 's/C213 DPPA/C113 PA  /' tail2_tmp >tail2
    sed 's/C214 DPPA/C114 PA  /' tail2 >tail2_tmp
    sed 's/C215 DPPA/C115 PA  /' tail2_tmp >tail2
    sed 's/C216 DPPA/C116 PA  /' tail2 >tail2_tmp
    sed 's/ H2R DPPA/ H2R PA  /' tail2_tmp >tail2
    sed 's/ H2S DPPA/ H2S PA  /' tail2 >tail2_tmp
    sed 's/ H3R DPPA/ H3R PA  /' tail2_tmp >tail2
    sed 's/ H3S DPPA/ H3S PA  /' tail2 >tail2_tmp
    sed 's/ H4R DPPA/ H4R PA  /' tail2_tmp >tail2
    sed 's/ H4S DPPA/ H4S PA  /' tail2 >tail2_tmp
    sed 's/ H5R DPPA/ H5R PA  /' tail2_tmp >tail2
    sed 's/ H5S DPPA/ H5S PA  /' tail2 >tail2_tmp
    sed 's/ H6R DPPA/ H6R PA  /' tail2_tmp >tail2
    sed 's/ H6S DPPA/ H6S PA  /' tail2 >tail2_tmp
    sed 's/ H7R DPPA/ H7R PA  /' tail2_tmp >tail2
    sed 's/ H7S DPPA/ H7S PA  /' tail2 >tail2_tmp
    sed 's/ H8R DPPA/ H8R PA  /' tail2_tmp >tail2
    sed 's/ H8S DPPA/ H8S PA  /' tail2 >tail2_tmp
    sed 's/ H9R DPPA/ H9R PA  /' tail2_tmp >tail2
    sed 's/ H9S DPPA/ H9S PA  /' tail2 >tail2_tmp
    sed 's/H10R DPPA/H10R PA  /' tail2_tmp >tail2
    sed 's/H10S DPPA/H10S PA  /' tail2 >tail2_tmp
    sed 's/H11R DPPA/H11R PA  /' tail2_tmp >tail2
    sed 's/H11S DPPA/H11S PA  /' tail2 >tail2_tmp
    sed 's/H12R DPPA/H12R PA  /' tail2_tmp >tail2
    sed 's/H12S DPPA/H12S PA  /' tail2 >tail2_tmp
    sed 's/H13R DPPA/H13R PA  /' tail2_tmp >tail2
    sed 's/H13S DPPA/H13S PA  /' tail2 >tail2_tmp
    sed 's/H14R DPPA/H14R PA  /' tail2_tmp >tail2
    sed 's/H14S DPPA/H14S PA  /' tail2 >tail2_tmp
    sed 's/H15R DPPA/H15R PA  /' tail2_tmp >tail2
    sed 's/H15S DPPA/H15S PA  /' tail2 >tail2_tmp
    sed 's/H16R DPPA/H16R PA  /' tail2_tmp >tail2
    sed 's/H16S DPPA/H16S PA  /' tail2 >tail2_tmp
    sed 's/H16T DPPA/H16T PA  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ P   DPPA/ P31 PH- /' head > head_tmp
    sed 's/ O13 DPPA/ O33 PH- /' head_tmp > head
    sed 's/ O14 DPPA/ O34 PH- /' head > head_tmp
    sed 's/ O12 DPPA/ O32 PH- /' head_tmp > head
    sed 's/ H12 DPPA/HO2A PH- /' head > head_tmp
    sed 's/ O11 DPPA/ O31 PH- /' head_tmp > head
    sed 's/ C1  DPPA/ C3  PH- /' head > head_tmp
    sed 's/ HA  DPPA/ HA  PH- /' head_tmp > head
    sed 's/ HB  DPPA/ HB  PH- /' head > head_tmp
    sed 's/ C2  DPPA/ C2  PH- /' head_tmp > head
    sed 's/ HS  DPPA/ HX  PH- /' head > head_tmp
    sed 's/ O21 DPPA/ O21 PH- /' head_tmp > head
    sed 's/ C21 DPPA/ C21 PH- /' head > head_tmp
    sed 's/ O22 DPPA/ O22 PH- /' head_tmp > head
    sed 's/ C3  DPPA/ C1  PH- /' head > head_tmp
    sed 's/ HX  DPPA/ HR  PH- /' head_tmp > head
    sed 's/ HY  DPPA/ HS  PH- /' head > head_tmp
    sed 's/ O31 DPPA/ O11 PH- /' head_tmp > head
    sed 's/ C31 DPPA/ C11 PH- /' head > head_tmp
    sed 's/ O32 DPPA/ O12 PH- /' head_tmp > head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DOPC
if ( $dopc_found == 1 ) then
  #DOPC is 138 atoms per lipid. 
  set atoms_per_lipid = 138

  echo " Processing: DOPC"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DOPC_tmp | awk '{print$1}'`
  echo "  DOPC Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DOPC Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DOPC Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DOPC_tmp
  rm -f ../DOPC_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DOPC Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DOPC" $i > tail1
    grep -A46 "C33 DOPC" $i >> tail1

    grep -A2 "C22 DOPC" $i > tail2
    grep -A46 "C23 DOPC" $i >> tail2

    head -n32 $i > head
    grep -A5 "C3  DOPC" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DOPC/ C12 OL  /' tail1 >tail1_tmp
    sed 's/ C33 DOPC/ C13 OL  /' tail1_tmp >tail1
    sed 's/ C34 DOPC/ C14 OL  /' tail1 >tail1_tmp
    sed 's/ C35 DOPC/ C15 OL  /' tail1_tmp >tail1
    sed 's/ C36 DOPC/ C16 OL  /' tail1 >tail1_tmp
    sed 's/ C37 DOPC/ C17 OL  /' tail1_tmp >tail1
    sed 's/ C38 DOPC/ C18 OL  /' tail1 >tail1_tmp
    sed 's/ C39 DOPC/ C19 OL  /' tail1_tmp >tail1
    sed 's/C310 DOPC/C110 OL  /' tail1 >tail1_tmp
    sed 's/C311 DOPC/C111 OL  /' tail1_tmp >tail1
    sed 's/C312 DOPC/C112 OL  /' tail1 >tail1_tmp
    sed 's/C313 DOPC/C113 OL  /' tail1_tmp >tail1
    sed 's/C314 DOPC/C114 OL  /' tail1 >tail1_tmp
    sed 's/C315 DOPC/C115 OL  /' tail1_tmp >tail1
    sed 's/C316 DOPC/C116 OL  /' tail1 >tail1_tmp
    sed 's/C317 DOPC/C117 OL  /' tail1_tmp >tail1
    sed 's/C318 DOPC/C118 OL  /' tail1 >tail1_tmp
    sed 's/ H2X DOPC/ H2R OL  /' tail1_tmp >tail1
    sed 's/ H2Y DOPC/ H2S OL  /' tail1 >tail1_tmp
    sed 's/ H3X DOPC/ H3R OL  /' tail1_tmp >tail1
    sed 's/ H3Y DOPC/ H3S OL  /' tail1 >tail1_tmp
    sed 's/ H4X DOPC/ H4R OL  /' tail1_tmp >tail1
    sed 's/ H4Y DOPC/ H4S OL  /' tail1 >tail1_tmp
    sed 's/ H5X DOPC/ H5R OL  /' tail1_tmp >tail1
    sed 's/ H5Y DOPC/ H5S OL  /' tail1 >tail1_tmp
    sed 's/ H6X DOPC/ H6R OL  /' tail1_tmp >tail1
    sed 's/ H6Y DOPC/ H6S OL  /' tail1 >tail1_tmp
    sed 's/ H7X DOPC/ H7R OL  /' tail1_tmp >tail1
    sed 's/ H7Y DOPC/ H7S OL  /' tail1 >tail1_tmp
    sed 's/ H8X DOPC/ H8R OL  /' tail1_tmp >tail1
    sed 's/ H8Y DOPC/ H8S OL  /' tail1 >tail1_tmp
    sed 's/ H9X DOPC/ H9R OL  /' tail1_tmp >tail1
    sed 's/H10X DOPC/H10R OL  /' tail1 >tail1_tmp
    sed 's/H11X DOPC/H11R OL  /' tail1_tmp >tail1
    sed 's/H11Y DOPC/H11S OL  /' tail1 >tail1_tmp
    sed 's/H12X DOPC/H12R OL  /' tail1_tmp >tail1
    sed 's/H12Y DOPC/H12S OL  /' tail1 >tail1_tmp
    sed 's/H13X DOPC/H13R OL  /' tail1_tmp >tail1
    sed 's/H13Y DOPC/H13S OL  /' tail1 >tail1_tmp
    sed 's/H14X DOPC/H14R OL  /' tail1_tmp >tail1
    sed 's/H14Y DOPC/H14S OL  /' tail1 >tail1_tmp
    sed 's/H15X DOPC/H15R OL  /' tail1_tmp >tail1
    sed 's/H15Y DOPC/H15S OL  /' tail1 >tail1_tmp
    sed 's/H16X DOPC/H16R OL  /' tail1_tmp >tail1
    sed 's/H16Y DOPC/H16S OL  /' tail1 >tail1_tmp
    sed 's/H17X DOPC/H17R OL  /' tail1_tmp >tail1
    sed 's/H17Y DOPC/H17S OL  /' tail1 >tail1_tmp
    sed 's/H18X DOPC/H18R OL  /' tail1_tmp >tail1
    sed 's/H18Y DOPC/H18S OL  /' tail1 >tail1_tmp
    sed 's/H18Z DOPC/H18T OL  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DOPC/ C12 OL  /' tail2 >tail2_tmp
    sed 's/ C23 DOPC/ C13 OL  /' tail2_tmp >tail2
    sed 's/ C24 DOPC/ C14 OL  /' tail2 >tail2_tmp
    sed 's/ C25 DOPC/ C15 OL  /' tail2_tmp >tail2
    sed 's/ C26 DOPC/ C16 OL  /' tail2 >tail2_tmp
    sed 's/ C27 DOPC/ C17 OL  /' tail2_tmp >tail2
    sed 's/ C28 DOPC/ C18 OL  /' tail2 >tail2_tmp
    sed 's/ C29 DOPC/ C19 OL  /' tail2_tmp >tail2
    sed 's/C210 DOPC/C110 OL  /' tail2 >tail2_tmp
    sed 's/C211 DOPC/C111 OL  /' tail2_tmp >tail2
    sed 's/C212 DOPC/C112 OL  /' tail2 >tail2_tmp
    sed 's/C213 DOPC/C113 OL  /' tail2_tmp >tail2
    sed 's/C214 DOPC/C114 OL  /' tail2 >tail2_tmp
    sed 's/C215 DOPC/C115 OL  /' tail2_tmp >tail2
    sed 's/C216 DOPC/C116 OL  /' tail2 >tail2_tmp
    sed 's/C217 DOPC/C117 OL  /' tail2_tmp >tail2
    sed 's/C218 DOPC/C118 OL  /' tail2 >tail2_tmp
    sed 's/ H2R DOPC/ H2R OL  /' tail2_tmp >tail2
    sed 's/ H2S DOPC/ H2S OL  /' tail2 >tail2_tmp
    sed 's/ H3R DOPC/ H3R OL  /' tail2_tmp >tail2
    sed 's/ H3S DOPC/ H3S OL  /' tail2 >tail2_tmp
    sed 's/ H4R DOPC/ H4R OL  /' tail2_tmp >tail2
    sed 's/ H4S DOPC/ H4S OL  /' tail2 >tail2_tmp
    sed 's/ H5R DOPC/ H5R OL  /' tail2_tmp >tail2
    sed 's/ H5S DOPC/ H5S OL  /' tail2 >tail2_tmp
    sed 's/ H6R DOPC/ H6R OL  /' tail2_tmp >tail2
    sed 's/ H6S DOPC/ H6S OL  /' tail2 >tail2_tmp
    sed 's/ H7R DOPC/ H7R OL  /' tail2_tmp >tail2
    sed 's/ H7S DOPC/ H7S OL  /' tail2 >tail2_tmp
    sed 's/ H8R DOPC/ H8R OL  /' tail2_tmp >tail2
    sed 's/ H8S DOPC/ H8S OL  /' tail2 >tail2_tmp
    sed 's/ H9R DOPC/ H9R OL  /' tail2_tmp >tail2
    sed 's/H10R DOPC/H10R OL  /' tail2 >tail2_tmp
    sed 's/H11R DOPC/H11R OL  /' tail2_tmp >tail2
    sed 's/H11S DOPC/H11S OL  /' tail2 >tail2_tmp
    sed 's/H12R DOPC/H12R OL  /' tail2_tmp >tail2
    sed 's/H12S DOPC/H12S OL  /' tail2 >tail2_tmp
    sed 's/H13R DOPC/H13R OL  /' tail2_tmp >tail2
    sed 's/H13S DOPC/H13S OL  /' tail2 >tail2_tmp
    sed 's/H14R DOPC/H14R OL  /' tail2_tmp >tail2
    sed 's/H14S DOPC/H14S OL  /' tail2 >tail2_tmp
    sed 's/H15R DOPC/H15R OL  /' tail2_tmp >tail2
    sed 's/H15S DOPC/H15S OL  /' tail2 >tail2_tmp
    sed 's/H16R DOPC/H16R OL  /' tail2_tmp >tail2
    sed 's/H16S DOPC/H16S OL  /' tail2 >tail2_tmp
    sed 's/H17R DOPC/H17R OL  /' tail2_tmp >tail2
    sed 's/H17S DOPC/H17S OL  /' tail2 >tail2_tmp
    sed 's/H18R DOPC/H18R OL  /' tail2_tmp >tail2
    sed 's/H18S DOPC/H18S OL  /' tail2 >tail2_tmp
    sed 's/H18T DOPC/H18T OL  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DOPC/ N31 PC  /' head > head_tmp
    sed 's/ C13 DOPC/ C33 PC  /' head_tmp > head
    sed 's/ C14 DOPC/ C34 PC  /' head > head_tmp
    sed 's/ C15 DOPC/ C35 PC  /' head_tmp > head
    sed 's/ C12 DOPC/ C32 PC  /' head > head_tmp
    sed 's/ C11 DOPC/ C31 PC  /' head_tmp > head
    sed 's/ O12 DOPC/ O32 PC  /' head > head_tmp
    sed 's/ P   DOPC/ P31 PC  /' head_tmp > head
    sed 's/ O11 DOPC/ O31 PC  /' head > head_tmp
    sed 's/ O13 DOPC/ O33 PC  /' head_tmp > head
    sed 's/ O14 DOPC/ O34 PC  /' head > head_tmp
    sed 's/ C1  DOPC/ C3  PC  /' head_tmp > head
    sed 's/ C2  DOPC/ C2  PC  /' head > head_tmp
    sed 's/ C3  DOPC/ C1  PC  /' head_tmp > head
    sed 's/ O31 DOPC/ O11 PC  /' head > head_tmp
    sed 's/ C31 DOPC/ C11 PC  /' head_tmp > head
    sed 's/ O32 DOPC/ O12 PC  /' head > head_tmp
    sed 's/ O21 DOPC/ O21 PC  /' head_tmp > head
    sed 's/ C21 DOPC/ C21 PC  /' head > head_tmp
    sed 's/ O22 DOPC/ O22 PC  /' head_tmp > head
    sed 's/ HX  DOPC/ HR  PC  /' head > head_tmp
    sed 's/ HY  DOPC/ HS  PC  /' head_tmp > head
    sed 's/ HS  DOPC/ HX  PC  /' head > head_tmp
    sed 's/ HA  DOPC/ HA  PC  /' head_tmp > head
    sed 's/ HB  DOPC/ HB  PC  /' head > head_tmp
    sed 's/H11A DOPC/ H1A PC  /' head_tmp > head
    sed 's/H11B DOPC/ H1B PC  /' head > head_tmp
    sed 's/H12A DOPC/ H2A PC  /' head_tmp > head
    sed 's/H12B DOPC/ H2B PC  /' head > head_tmp
    sed 's/H13A DOPC/ H3A PC  /' head_tmp > head
    sed 's/H13B DOPC/ H3B PC  /' head > head_tmp
    sed 's/H13C DOPC/ H3C PC  /' head_tmp > head
    sed 's/H14A DOPC/ H4A PC  /' head > head_tmp
    sed 's/H14B DOPC/ H4B PC  /' head_tmp > head
    sed 's/H14C DOPC/ H4C PC  /' head > head_tmp
    sed 's/H15A DOPC/ H5A PC  /' head_tmp > head
    sed 's/H15B DOPC/ H5B PC  /' head > head_tmp
    sed 's/H15C DOPC/ H5C PC  /' head_tmp > head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DOPE
if ( $dope_found == 1 ) then
  #DOPE is 129 atoms per lipid. 
  set atoms_per_lipid = 129

  echo " Processing: DOPE"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DOPE_tmp | awk '{print$1}'`
  echo "  DOPE Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DOPE Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DOPE Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DOPE_tmp
  rm -f ../DOPE_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DOPE Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DOPE" $i > tail1
    grep -A46 "C33 DOPE" $i >> tail1

    grep -A2 "C22 DOPE" $i > tail2
    grep -A46 "C23 DOPE" $i >> tail2

    head -n23 $i > head
    grep -A5 "C3  DOPE" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DOPE/ C12 OL  /' tail1 >tail1_tmp
    sed 's/ C33 DOPE/ C13 OL  /' tail1_tmp >tail1
    sed 's/ C34 DOPE/ C14 OL  /' tail1 >tail1_tmp
    sed 's/ C35 DOPE/ C15 OL  /' tail1_tmp >tail1
    sed 's/ C36 DOPE/ C16 OL  /' tail1 >tail1_tmp
    sed 's/ C37 DOPE/ C17 OL  /' tail1_tmp >tail1
    sed 's/ C38 DOPE/ C18 OL  /' tail1 >tail1_tmp
    sed 's/ C39 DOPE/ C19 OL  /' tail1_tmp >tail1
    sed 's/C310 DOPE/C110 OL  /' tail1 >tail1_tmp
    sed 's/C311 DOPE/C111 OL  /' tail1_tmp >tail1
    sed 's/C312 DOPE/C112 OL  /' tail1 >tail1_tmp
    sed 's/C313 DOPE/C113 OL  /' tail1_tmp >tail1
    sed 's/C314 DOPE/C114 OL  /' tail1 >tail1_tmp
    sed 's/C315 DOPE/C115 OL  /' tail1_tmp >tail1
    sed 's/C316 DOPE/C116 OL  /' tail1 >tail1_tmp
    sed 's/C317 DOPE/C117 OL  /' tail1_tmp >tail1
    sed 's/C318 DOPE/C118 OL  /' tail1 >tail1_tmp
    sed 's/ H2X DOPE/ H2R OL  /' tail1_tmp >tail1
    sed 's/ H2Y DOPE/ H2S OL  /' tail1 >tail1_tmp
    sed 's/ H3X DOPE/ H3R OL  /' tail1_tmp >tail1
    sed 's/ H3Y DOPE/ H3S OL  /' tail1 >tail1_tmp
    sed 's/ H4X DOPE/ H4R OL  /' tail1_tmp >tail1
    sed 's/ H4Y DOPE/ H4S OL  /' tail1 >tail1_tmp
    sed 's/ H5X DOPE/ H5R OL  /' tail1_tmp >tail1
    sed 's/ H5Y DOPE/ H5S OL  /' tail1 >tail1_tmp
    sed 's/ H6X DOPE/ H6R OL  /' tail1_tmp >tail1
    sed 's/ H6Y DOPE/ H6S OL  /' tail1 >tail1_tmp
    sed 's/ H7X DOPE/ H7R OL  /' tail1_tmp >tail1
    sed 's/ H7Y DOPE/ H7S OL  /' tail1 >tail1_tmp
    sed 's/ H8X DOPE/ H8R OL  /' tail1_tmp >tail1
    sed 's/ H8Y DOPE/ H8S OL  /' tail1 >tail1_tmp
    sed 's/ H9X DOPE/ H9R OL  /' tail1_tmp >tail1
    sed 's/H10X DOPE/H10R OL  /' tail1 >tail1_tmp
    sed 's/H11X DOPE/H11R OL  /' tail1_tmp >tail1
    sed 's/H11Y DOPE/H11S OL  /' tail1 >tail1_tmp
    sed 's/H12X DOPE/H12R OL  /' tail1_tmp >tail1
    sed 's/H12Y DOPE/H12S OL  /' tail1 >tail1_tmp
    sed 's/H13X DOPE/H13R OL  /' tail1_tmp >tail1
    sed 's/H13Y DOPE/H13S OL  /' tail1 >tail1_tmp
    sed 's/H14X DOPE/H14R OL  /' tail1_tmp >tail1
    sed 's/H14Y DOPE/H14S OL  /' tail1 >tail1_tmp
    sed 's/H15X DOPE/H15R OL  /' tail1_tmp >tail1
    sed 's/H15Y DOPE/H15S OL  /' tail1 >tail1_tmp
    sed 's/H16X DOPE/H16R OL  /' tail1_tmp >tail1
    sed 's/H16Y DOPE/H16S OL  /' tail1 >tail1_tmp
    sed 's/H17X DOPE/H17R OL  /' tail1_tmp >tail1
    sed 's/H17Y DOPE/H17S OL  /' tail1 >tail1_tmp
    sed 's/H18X DOPE/H18R OL  /' tail1_tmp >tail1
    sed 's/H18Y DOPE/H18S OL  /' tail1 >tail1_tmp
    sed 's/H18Z DOPE/H18T OL  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DOPE/ C12 OL  /' tail2 >tail2_tmp
    sed 's/ C23 DOPE/ C13 OL  /' tail2_tmp >tail2
    sed 's/ C24 DOPE/ C14 OL  /' tail2 >tail2_tmp
    sed 's/ C25 DOPE/ C15 OL  /' tail2_tmp >tail2
    sed 's/ C26 DOPE/ C16 OL  /' tail2 >tail2_tmp
    sed 's/ C27 DOPE/ C17 OL  /' tail2_tmp >tail2
    sed 's/ C28 DOPE/ C18 OL  /' tail2 >tail2_tmp
    sed 's/ C29 DOPE/ C19 OL  /' tail2_tmp >tail2
    sed 's/C210 DOPE/C110 OL  /' tail2 >tail2_tmp
    sed 's/C211 DOPE/C111 OL  /' tail2_tmp >tail2
    sed 's/C212 DOPE/C112 OL  /' tail2 >tail2_tmp
    sed 's/C213 DOPE/C113 OL  /' tail2_tmp >tail2
    sed 's/C214 DOPE/C114 OL  /' tail2 >tail2_tmp
    sed 's/C215 DOPE/C115 OL  /' tail2_tmp >tail2
    sed 's/C216 DOPE/C116 OL  /' tail2 >tail2_tmp
    sed 's/C217 DOPE/C117 OL  /' tail2_tmp >tail2
    sed 's/C218 DOPE/C118 OL  /' tail2 >tail2_tmp
    sed 's/ H2R DOPE/ H2R OL  /' tail2_tmp >tail2
    sed 's/ H2S DOPE/ H2S OL  /' tail2 >tail2_tmp
    sed 's/ H3R DOPE/ H3R OL  /' tail2_tmp >tail2
    sed 's/ H3S DOPE/ H3S OL  /' tail2 >tail2_tmp
    sed 's/ H4R DOPE/ H4R OL  /' tail2_tmp >tail2
    sed 's/ H4S DOPE/ H4S OL  /' tail2 >tail2_tmp
    sed 's/ H5R DOPE/ H5R OL  /' tail2_tmp >tail2
    sed 's/ H5S DOPE/ H5S OL  /' tail2 >tail2_tmp
    sed 's/ H6R DOPE/ H6R OL  /' tail2_tmp >tail2
    sed 's/ H6S DOPE/ H6S OL  /' tail2 >tail2_tmp
    sed 's/ H7R DOPE/ H7R OL  /' tail2_tmp >tail2
    sed 's/ H7S DOPE/ H7S OL  /' tail2 >tail2_tmp
    sed 's/ H8R DOPE/ H8R OL  /' tail2_tmp >tail2
    sed 's/ H8S DOPE/ H8S OL  /' tail2 >tail2_tmp
    sed 's/ H9R DOPE/ H9R OL  /' tail2_tmp >tail2
    sed 's/H10R DOPE/H10R OL  /' tail2 >tail2_tmp
    sed 's/H11R DOPE/H11R OL  /' tail2_tmp >tail2
    sed 's/H11S DOPE/H11S OL  /' tail2 >tail2_tmp
    sed 's/H12R DOPE/H12R OL  /' tail2_tmp >tail2
    sed 's/H12S DOPE/H12S OL  /' tail2 >tail2_tmp
    sed 's/H13R DOPE/H13R OL  /' tail2_tmp >tail2
    sed 's/H13S DOPE/H13S OL  /' tail2 >tail2_tmp
    sed 's/H14R DOPE/H14R OL  /' tail2_tmp >tail2
    sed 's/H14S DOPE/H14S OL  /' tail2 >tail2_tmp
    sed 's/H15R DOPE/H15R OL  /' tail2_tmp >tail2
    sed 's/H15S DOPE/H15S OL  /' tail2 >tail2_tmp
    sed 's/H16R DOPE/H16R OL  /' tail2_tmp >tail2
    sed 's/H16S DOPE/H16S OL  /' tail2 >tail2_tmp
    sed 's/H17R DOPE/H17R OL  /' tail2_tmp >tail2
    sed 's/H17S DOPE/H17S OL  /' tail2 >tail2_tmp
    sed 's/H18R DOPE/H18R OL  /' tail2_tmp >tail2
    sed 's/H18S DOPE/H18S OL  /' tail2 >tail2_tmp
    sed 's/H18T DOPE/H18T OL  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DOPE/ N31 PE  /' head > head_tmp
    sed 's/ C12 DOPE/ C32 PE  /' head_tmp > head
    sed 's/ C11 DOPE/ C31 PE  /' head > head_tmp
    sed 's/ O12 DOPE/ O32 PE  /' head_tmp > head
    sed 's/ P   DOPE/ P31 PE  /' head > head_tmp
    sed 's/ O11 DOPE/ O31 PE  /' head_tmp > head
    sed 's/ O13 DOPE/ O33 PE  /' head > head_tmp
    sed 's/ O14 DOPE/ O34 PE  /' head_tmp > head
    sed 's/ C1  DOPE/ C3  PE  /' head > head_tmp
    sed 's/ C2  DOPE/ C2  PE  /' head_tmp > head
    sed 's/ C3  DOPE/ C1  PE  /' head > head_tmp
    sed 's/ O31 DOPE/ O11 PE  /' head_tmp > head
    sed 's/ C31 DOPE/ C11 PE  /' head > head_tmp
    sed 's/ O32 DOPE/ O12 PE  /' head_tmp > head
    sed 's/ O21 DOPE/ O21 PE  /' head > head_tmp
    sed 's/ C21 DOPE/ C21 PE  /' head_tmp > head
    sed 's/ O22 DOPE/ O22 PE  /' head > head_tmp
    sed 's/ HX  DOPE/ HR  PE  /' head_tmp > head
    sed 's/ HY  DOPE/ HS  PE  /' head > head_tmp
    sed 's/ HS  DOPE/ HX  PE  /' head_tmp > head
    sed 's/ HA  DOPE/ HA  PE  /' head > head_tmp
    sed 's/ HB  DOPE/ HB  PE  /' head_tmp > head
    sed 's/H11A DOPE/ H1A PE  /' head > head_tmp
    sed 's/H11B DOPE/ H1B PE  /' head_tmp > head
    sed 's/H12A DOPE/ H2A PE  /' head > head_tmp
    sed 's/H12B DOPE/ H2B PE  /' head_tmp > head
    sed 's/ HN1 DOPE/HN1A PE  /' head > head_tmp
    sed 's/ HN2 DOPE/HN1B PE  /' head_tmp > head
    sed 's/ HN3 DOPE/HN1C PE  /' head > head_tmp
    mv head_tmp head  
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DOPS
if ( $dops_found == 1 ) then
  #DOPS is 131 atoms per lipid. 
  set atoms_per_lipid = 131

  echo " Processing: DOPS"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DOPS_tmp | awk '{print$1}'`
  echo "  DOPS Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DOPS Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DOPS Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DOPS_tmp
  rm -f ../DOPS_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DOPS Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DOPS" $i > tail1
    grep -A46 "C33 DOPS" $i >> tail1

    grep -A2 "C22 DOPS" $i > tail2
    grep -A46 "C23 DOPS" $i >> tail2

    head -n25 $i > head
    grep -A5 "C3  DOPS" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DOPS/ C12 OL  /' tail1 >tail1_tmp
    sed 's/ C33 DOPS/ C13 OL  /' tail1_tmp >tail1
    sed 's/ C34 DOPS/ C14 OL  /' tail1 >tail1_tmp
    sed 's/ C35 DOPS/ C15 OL  /' tail1_tmp >tail1
    sed 's/ C36 DOPS/ C16 OL  /' tail1 >tail1_tmp
    sed 's/ C37 DOPS/ C17 OL  /' tail1_tmp >tail1
    sed 's/ C38 DOPS/ C18 OL  /' tail1 >tail1_tmp
    sed 's/ C39 DOPS/ C19 OL  /' tail1_tmp >tail1
    sed 's/C310 DOPS/C110 OL  /' tail1 >tail1_tmp
    sed 's/C311 DOPS/C111 OL  /' tail1_tmp >tail1
    sed 's/C312 DOPS/C112 OL  /' tail1 >tail1_tmp
    sed 's/C313 DOPS/C113 OL  /' tail1_tmp >tail1
    sed 's/C314 DOPS/C114 OL  /' tail1 >tail1_tmp
    sed 's/C315 DOPS/C115 OL  /' tail1_tmp >tail1
    sed 's/C316 DOPS/C116 OL  /' tail1 >tail1_tmp
    sed 's/C317 DOPS/C117 OL  /' tail1_tmp >tail1
    sed 's/C318 DOPS/C118 OL  /' tail1 >tail1_tmp
    sed 's/ H2X DOPS/ H2R OL  /' tail1_tmp >tail1
    sed 's/ H2Y DOPS/ H2S OL  /' tail1 >tail1_tmp
    sed 's/ H3X DOPS/ H3R OL  /' tail1_tmp >tail1
    sed 's/ H3Y DOPS/ H3S OL  /' tail1 >tail1_tmp
    sed 's/ H4X DOPS/ H4R OL  /' tail1_tmp >tail1
    sed 's/ H4Y DOPS/ H4S OL  /' tail1 >tail1_tmp
    sed 's/ H5X DOPS/ H5R OL  /' tail1_tmp >tail1
    sed 's/ H5Y DOPS/ H5S OL  /' tail1 >tail1_tmp
    sed 's/ H6X DOPS/ H6R OL  /' tail1_tmp >tail1
    sed 's/ H6Y DOPS/ H6S OL  /' tail1 >tail1_tmp
    sed 's/ H7X DOPS/ H7R OL  /' tail1_tmp >tail1
    sed 's/ H7Y DOPS/ H7S OL  /' tail1 >tail1_tmp
    sed 's/ H8X DOPS/ H8R OL  /' tail1_tmp >tail1
    sed 's/ H8Y DOPS/ H8S OL  /' tail1 >tail1_tmp
    sed 's/ H9X DOPS/ H9R OL  /' tail1_tmp >tail1
    sed 's/H10X DOPS/H10R OL  /' tail1 >tail1_tmp
    sed 's/H11X DOPS/H11R OL  /' tail1_tmp >tail1
    sed 's/H11Y DOPS/H11S OL  /' tail1 >tail1_tmp
    sed 's/H12X DOPS/H12R OL  /' tail1_tmp >tail1
    sed 's/H12Y DOPS/H12S OL  /' tail1 >tail1_tmp
    sed 's/H13X DOPS/H13R OL  /' tail1_tmp >tail1
    sed 's/H13Y DOPS/H13S OL  /' tail1 >tail1_tmp
    sed 's/H14X DOPS/H14R OL  /' tail1_tmp >tail1
    sed 's/H14Y DOPS/H14S OL  /' tail1 >tail1_tmp
    sed 's/H15X DOPS/H15R OL  /' tail1_tmp >tail1
    sed 's/H15Y DOPS/H15S OL  /' tail1 >tail1_tmp
    sed 's/H16X DOPS/H16R OL  /' tail1_tmp >tail1
    sed 's/H16Y DOPS/H16S OL  /' tail1 >tail1_tmp
    sed 's/H17X DOPS/H17R OL  /' tail1_tmp >tail1
    sed 's/H17Y DOPS/H17S OL  /' tail1 >tail1_tmp
    sed 's/H18X DOPS/H18R OL  /' tail1_tmp >tail1
    sed 's/H18Y DOPS/H18S OL  /' tail1 >tail1_tmp
    sed 's/H18Z DOPS/H18T OL  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DOPS/ C12 OL  /' tail2 >tail2_tmp
    sed 's/ C23 DOPS/ C13 OL  /' tail2_tmp >tail2
    sed 's/ C24 DOPS/ C14 OL  /' tail2 >tail2_tmp
    sed 's/ C25 DOPS/ C15 OL  /' tail2_tmp >tail2
    sed 's/ C26 DOPS/ C16 OL  /' tail2 >tail2_tmp
    sed 's/ C27 DOPS/ C17 OL  /' tail2_tmp >tail2
    sed 's/ C28 DOPS/ C18 OL  /' tail2 >tail2_tmp
    sed 's/ C29 DOPS/ C19 OL  /' tail2_tmp >tail2
    sed 's/C210 DOPS/C110 OL  /' tail2 >tail2_tmp
    sed 's/C211 DOPS/C111 OL  /' tail2_tmp >tail2
    sed 's/C212 DOPS/C112 OL  /' tail2 >tail2_tmp
    sed 's/C213 DOPS/C113 OL  /' tail2_tmp >tail2
    sed 's/C214 DOPS/C114 OL  /' tail2 >tail2_tmp
    sed 's/C215 DOPS/C115 OL  /' tail2_tmp >tail2
    sed 's/C216 DOPS/C116 OL  /' tail2 >tail2_tmp
    sed 's/C217 DOPS/C117 OL  /' tail2_tmp >tail2
    sed 's/C218 DOPS/C118 OL  /' tail2 >tail2_tmp
    sed 's/ H2R DOPS/ H2R OL  /' tail2_tmp >tail2
    sed 's/ H2S DOPS/ H2S OL  /' tail2 >tail2_tmp
    sed 's/ H3R DOPS/ H3R OL  /' tail2_tmp >tail2
    sed 's/ H3S DOPS/ H3S OL  /' tail2 >tail2_tmp
    sed 's/ H4R DOPS/ H4R OL  /' tail2_tmp >tail2
    sed 's/ H4S DOPS/ H4S OL  /' tail2 >tail2_tmp
    sed 's/ H5R DOPS/ H5R OL  /' tail2_tmp >tail2
    sed 's/ H5S DOPS/ H5S OL  /' tail2 >tail2_tmp
    sed 's/ H6R DOPS/ H6R OL  /' tail2_tmp >tail2
    sed 's/ H6S DOPS/ H6S OL  /' tail2 >tail2_tmp
    sed 's/ H7R DOPS/ H7R OL  /' tail2_tmp >tail2
    sed 's/ H7S DOPS/ H7S OL  /' tail2 >tail2_tmp
    sed 's/ H8R DOPS/ H8R OL  /' tail2_tmp >tail2
    sed 's/ H8S DOPS/ H8S OL  /' tail2 >tail2_tmp
    sed 's/ H9R DOPS/ H9R OL  /' tail2_tmp >tail2
    sed 's/H10R DOPS/H10R OL  /' tail2 >tail2_tmp
    sed 's/H11R DOPS/H11R OL  /' tail2_tmp >tail2
    sed 's/H11S DOPS/H11S OL  /' tail2 >tail2_tmp
    sed 's/H12R DOPS/H12R OL  /' tail2_tmp >tail2
    sed 's/H12S DOPS/H12S OL  /' tail2 >tail2_tmp
    sed 's/H13R DOPS/H13R OL  /' tail2_tmp >tail2
    sed 's/H13S DOPS/H13S OL  /' tail2 >tail2_tmp
    sed 's/H14R DOPS/H14R OL  /' tail2_tmp >tail2
    sed 's/H14S DOPS/H14S OL  /' tail2 >tail2_tmp
    sed 's/H15R DOPS/H15R OL  /' tail2_tmp >tail2
    sed 's/H15S DOPS/H15S OL  /' tail2 >tail2_tmp
    sed 's/H16R DOPS/H16R OL  /' tail2_tmp >tail2
    sed 's/H16S DOPS/H16S OL  /' tail2 >tail2_tmp
    sed 's/H17R DOPS/H17R OL  /' tail2_tmp >tail2
    sed 's/H17S DOPS/H17S OL  /' tail2 >tail2_tmp
    sed 's/H18R DOPS/H18R OL  /' tail2_tmp >tail2
    sed 's/H18S DOPS/H18S OL  /' tail2 >tail2_tmp
    sed 's/H18T DOPS/H18T OL  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DOPS/ N31 PS  /' head > head_tmp
    sed 's/ C12 DOPS/ C32 PS  /' head_tmp > head
    sed 's/ C13 DOPS/ C33 PS  /' head > head_tmp
    sed 's/O13A DOPS/ O35 PS  /' head_tmp > head
    sed 's/O13B DOPS/ O36 PS  /' head > head_tmp
    sed 's/ O12 DOPS/ O32 PS  /' head_tmp > head
    sed 's/ P   DOPS/ P31 PS  /' head > head_tmp
    sed 's/ O11 DOPS/ O31 PS  /' head_tmp > head
    sed 's/ O13 DOPS/ O33 PS  /' head > head_tmp
    sed 's/ C11 DOPS/ C31 PS  /' head_tmp > head
    sed 's/ O14 DOPS/ O34 PS  /' head > head_tmp
    sed 's/ C1  DOPS/ C3  PS  /' head_tmp > head
    sed 's/ C2  DOPS/ C2  PS  /' head > head_tmp
    sed 's/ C3  DOPS/ C1  PS  /' head_tmp > head
    sed 's/ O31 DOPS/ O11 PS  /' head > head_tmp
    sed 's/ C31 DOPS/ C11 PS  /' head_tmp > head
    sed 's/ O32 DOPS/ O12 PS  /' head > head_tmp
    sed 's/ O21 DOPS/ O21 PS  /' head_tmp > head
    sed 's/ C21 DOPS/ C21 PS  /' head > head_tmp
    sed 's/ O22 DOPS/ O22 PS  /' head_tmp > head
    sed 's/ HX  DOPS/ HR  PS  /' head > head_tmp
    sed 's/ HY  DOPS/ HS  PS  /' head_tmp > head
    sed 's/ HS  DOPS/ HX  PS  /' head > head_tmp
    sed 's/ HA  DOPS/ HA  PS  /' head_tmp > head
    sed 's/ HB  DOPS/ HB  PS  /' head > head_tmp
    sed 's/H11A DOPS/ H1A PS  /' head_tmp > head
    sed 's/H11B DOPS/ H1B PS  /' head > head_tmp
    sed 's/H12A DOPS/ H2A PS  /' head_tmp > head
    sed 's/ HN1 DOPS/HN1A PS  /' head > head_tmp
    sed 's/ HN2 DOPS/HN1B PS  /' head_tmp > head
    sed 's/ HN3 DOPS/HN1C PS  /' head > head_tmp
    mv head_tmp head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DOPG
if ( $dopg_found == 1 ) then
  #DOPG is 131 atoms per lipid. 
  set atoms_per_lipid = 131

  echo " Processing: DOPG"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DOPG_tmp | awk '{print$1}'`
  echo "  DOPG Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DOPG Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DOPG Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DOPG_tmp
  rm -f ../DOPG_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DOPG Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DOPG" $i > tail1
    grep -A46 "C33 DOPG" $i >> tail1

    grep -A2 "C22 DOPG" $i > tail2
    grep -A46 "C23 DOPG" $i >> tail2

    head -n25 $i > head
    grep -A5 "C3  DOPG" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DOPG/ C12 OL  /' tail1 >tail1_tmp
    sed 's/ C33 DOPG/ C13 OL  /' tail1_tmp >tail1
    sed 's/ C34 DOPG/ C14 OL  /' tail1 >tail1_tmp
    sed 's/ C35 DOPG/ C15 OL  /' tail1_tmp >tail1
    sed 's/ C36 DOPG/ C16 OL  /' tail1 >tail1_tmp
    sed 's/ C37 DOPG/ C17 OL  /' tail1_tmp >tail1
    sed 's/ C38 DOPG/ C18 OL  /' tail1 >tail1_tmp
    sed 's/ C39 DOPG/ C19 OL  /' tail1_tmp >tail1
    sed 's/C310 DOPG/C110 OL  /' tail1 >tail1_tmp
    sed 's/C311 DOPG/C111 OL  /' tail1_tmp >tail1
    sed 's/C312 DOPG/C112 OL  /' tail1 >tail1_tmp
    sed 's/C313 DOPG/C113 OL  /' tail1_tmp >tail1
    sed 's/C314 DOPG/C114 OL  /' tail1 >tail1_tmp
    sed 's/C315 DOPG/C115 OL  /' tail1_tmp >tail1
    sed 's/C316 DOPG/C116 OL  /' tail1 >tail1_tmp
    sed 's/C317 DOPG/C117 OL  /' tail1_tmp >tail1
    sed 's/C318 DOPG/C118 OL  /' tail1 >tail1_tmp
    sed 's/ H2X DOPG/ H2R OL  /' tail1_tmp >tail1
    sed 's/ H2Y DOPG/ H2S OL  /' tail1 >tail1_tmp
    sed 's/ H3X DOPG/ H3R OL  /' tail1_tmp >tail1
    sed 's/ H3Y DOPG/ H3S OL  /' tail1 >tail1_tmp
    sed 's/ H4X DOPG/ H4R OL  /' tail1_tmp >tail1
    sed 's/ H4Y DOPG/ H4S OL  /' tail1 >tail1_tmp
    sed 's/ H5X DOPG/ H5R OL  /' tail1_tmp >tail1
    sed 's/ H5Y DOPG/ H5S OL  /' tail1 >tail1_tmp
    sed 's/ H6X DOPG/ H6R OL  /' tail1_tmp >tail1
    sed 's/ H6Y DOPG/ H6S OL  /' tail1 >tail1_tmp
    sed 's/ H7X DOPG/ H7R OL  /' tail1_tmp >tail1
    sed 's/ H7Y DOPG/ H7S OL  /' tail1 >tail1_tmp
    sed 's/ H8X DOPG/ H8R OL  /' tail1_tmp >tail1
    sed 's/ H8Y DOPG/ H8S OL  /' tail1 >tail1_tmp
    sed 's/ H9X DOPG/ H9R OL  /' tail1_tmp >tail1
    sed 's/H10X DOPG/H10R OL  /' tail1 >tail1_tmp
    sed 's/H11X DOPG/H11R OL  /' tail1_tmp >tail1
    sed 's/H11Y DOPG/H11S OL  /' tail1 >tail1_tmp
    sed 's/H12X DOPG/H12R OL  /' tail1_tmp >tail1
    sed 's/H12Y DOPG/H12S OL  /' tail1 >tail1_tmp
    sed 's/H13X DOPG/H13R OL  /' tail1_tmp >tail1
    sed 's/H13Y DOPG/H13S OL  /' tail1 >tail1_tmp
    sed 's/H14X DOPG/H14R OL  /' tail1_tmp >tail1
    sed 's/H14Y DOPG/H14S OL  /' tail1 >tail1_tmp
    sed 's/H15X DOPG/H15R OL  /' tail1_tmp >tail1
    sed 's/H15Y DOPG/H15S OL  /' tail1 >tail1_tmp
    sed 's/H16X DOPG/H16R OL  /' tail1_tmp >tail1
    sed 's/H16Y DOPG/H16S OL  /' tail1 >tail1_tmp
    sed 's/H17X DOPG/H17R OL  /' tail1_tmp >tail1
    sed 's/H17Y DOPG/H17S OL  /' tail1 >tail1_tmp
    sed 's/H18X DOPG/H18R OL  /' tail1_tmp >tail1
    sed 's/H18Y DOPG/H18S OL  /' tail1 >tail1_tmp
    sed 's/H18Z DOPG/H18T OL  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DOPG/ C12 OL  /' tail2 >tail2_tmp
    sed 's/ C23 DOPG/ C13 OL  /' tail2_tmp >tail2
    sed 's/ C24 DOPG/ C14 OL  /' tail2 >tail2_tmp
    sed 's/ C25 DOPG/ C15 OL  /' tail2_tmp >tail2
    sed 's/ C26 DOPG/ C16 OL  /' tail2 >tail2_tmp
    sed 's/ C27 DOPG/ C17 OL  /' tail2_tmp >tail2
    sed 's/ C28 DOPG/ C18 OL  /' tail2 >tail2_tmp
    sed 's/ C29 DOPG/ C19 OL  /' tail2_tmp >tail2
    sed 's/C210 DOPG/C110 OL  /' tail2 >tail2_tmp
    sed 's/C211 DOPG/C111 OL  /' tail2_tmp >tail2
    sed 's/C212 DOPG/C112 OL  /' tail2 >tail2_tmp
    sed 's/C213 DOPG/C113 OL  /' tail2_tmp >tail2
    sed 's/C214 DOPG/C114 OL  /' tail2 >tail2_tmp
    sed 's/C215 DOPG/C115 OL  /' tail2_tmp >tail2
    sed 's/C216 DOPG/C116 OL  /' tail2 >tail2_tmp
    sed 's/C217 DOPG/C117 OL  /' tail2_tmp >tail2
    sed 's/C218 DOPG/C118 OL  /' tail2 >tail2_tmp
    sed 's/ H2R DOPG/ H2R OL  /' tail2_tmp >tail2
    sed 's/ H2S DOPG/ H2S OL  /' tail2 >tail2_tmp
    sed 's/ H3R DOPG/ H3R OL  /' tail2_tmp >tail2
    sed 's/ H3S DOPG/ H3S OL  /' tail2 >tail2_tmp
    sed 's/ H4R DOPG/ H4R OL  /' tail2_tmp >tail2
    sed 's/ H4S DOPG/ H4S OL  /' tail2 >tail2_tmp
    sed 's/ H5R DOPG/ H5R OL  /' tail2_tmp >tail2
    sed 's/ H5S DOPG/ H5S OL  /' tail2 >tail2_tmp
    sed 's/ H6R DOPG/ H6R OL  /' tail2_tmp >tail2
    sed 's/ H6S DOPG/ H6S OL  /' tail2 >tail2_tmp
    sed 's/ H7R DOPG/ H7R OL  /' tail2_tmp >tail2
    sed 's/ H7S DOPG/ H7S OL  /' tail2 >tail2_tmp
    sed 's/ H8R DOPG/ H8R OL  /' tail2_tmp >tail2
    sed 's/ H8S DOPG/ H8S OL  /' tail2 >tail2_tmp
    sed 's/ H9R DOPG/ H9R OL  /' tail2_tmp >tail2
    sed 's/H10R DOPG/H10R OL  /' tail2 >tail2_tmp
    sed 's/H11R DOPG/H11R OL  /' tail2_tmp >tail2
    sed 's/H11S DOPG/H11S OL  /' tail2 >tail2_tmp
    sed 's/H12R DOPG/H12R OL  /' tail2_tmp >tail2
    sed 's/H12S DOPG/H12S OL  /' tail2 >tail2_tmp
    sed 's/H13R DOPG/H13R OL  /' tail2_tmp >tail2
    sed 's/H13S DOPG/H13S OL  /' tail2 >tail2_tmp
    sed 's/H14R DOPG/H14R OL  /' tail2_tmp >tail2
    sed 's/H14S DOPG/H14S OL  /' tail2 >tail2_tmp
    sed 's/H15R DOPG/H15R OL  /' tail2_tmp >tail2
    sed 's/H15S DOPG/H15S OL  /' tail2 >tail2_tmp
    sed 's/H16R DOPG/H16R OL  /' tail2_tmp >tail2
    sed 's/H16S DOPG/H16S OL  /' tail2 >tail2_tmp
    sed 's/H17R DOPG/H17R OL  /' tail2_tmp >tail2
    sed 's/H17S DOPG/H17S OL  /' tail2 >tail2_tmp
    sed 's/H18R DOPG/H18R OL  /' tail2_tmp >tail2
    sed 's/H18S DOPG/H18S OL  /' tail2 >tail2_tmp
    sed 's/H18T DOPG/H18T OL  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ C13 DOPG/ C33 PGR /' head > head_tmp
    sed 's/H13A DOPG/ H3A PGR /' head_tmp > head
    sed 's/H13B DOPG/ H3B PGR /' head > head_tmp
    sed 's/ OC3 DOPG/ O36 PGR /' head_tmp > head
    sed 's/ HO3 DOPG/HO6A PGR /' head > head_tmp
    sed 's/ C12 DOPG/ C32 PGR /' head_tmp > head
    sed 's/H12A DOPG/ H2A PGR /' head > head_tmp
    sed 's/ OC2 DOPG/ O35 PGR /' head_tmp > head
    sed 's/ HO2 DOPG/HO5A PGR /' head > head_tmp
    sed 's/ C11 DOPG/ C31 PGR /' head_tmp > head
    sed 's/H11A DOPG/ H1A PGR /' head > head_tmp
    sed 's/H11B DOPG/ H1B PGR /' head_tmp > head
    sed 's/ P   DOPG/ P31 PGR /' head > head_tmp
    sed 's/ O13 DOPG/ O33 PGR /' head_tmp > head
    sed 's/ O14 DOPG/ O34 PGR /' head > head_tmp
    sed 's/ O12 DOPG/ O32 PGR /' head_tmp > head
    sed 's/ O11 DOPG/ O31 PGR /' head > head_tmp
    sed 's/ C1  DOPG/ C3  PGR /' head_tmp > head
    sed 's/ HA  DOPG/ HA  PGR /' head > head_tmp
    sed 's/ HB  DOPG/ HB  PGR /' head_tmp > head
    sed 's/ C2  DOPG/ C2  PGR /' head > head_tmp
    sed 's/ HS  DOPG/ HX  PGR /' head_tmp > head
    sed 's/ O21 DOPG/ O21 PGR /' head > head_tmp
    sed 's/ C21 DOPG/ C21 PGR /' head_tmp > head
    sed 's/ O22 DOPG/ O22 PGR /' head > head_tmp
    sed 's/ C3  DOPG/ C1  PGR /' head_tmp > head
    sed 's/ HX  DOPG/ HR  PGR /' head > head_tmp
    sed 's/ HY  DOPG/ HS  PGR /' head_tmp > head
    sed 's/ O31 DOPG/ O11 PGR /' head > head_tmp
    sed 's/ C31 DOPG/ C11 PGR /' head_tmp > head
    sed 's/ O32 DOPG/ O12 PGR /' head > head_tmp
    mv head_tmp head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DOPA
if ( $dopa_found == 1 ) then
  #DOPA is 120 atoms per lipid. 
  set atoms_per_lipid = 120

  echo " Processing: DOPA"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DOPA_tmp | awk '{print$1}'`
  echo "  DOPA Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DOPA Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DOPA Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DOPA_tmp
  rm -f ../DOPA_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DOPA Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DOPA" $i > tail1
    grep -A46 "C33 DOPA" $i >> tail1

    grep -A2 "C22 DOPA" $i > tail2
    grep -A46 "C23 DOPA" $i >> tail2

    head -n14 $i > head
    grep -A5 "C3  DOPA" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DOPA/ C12 OL  /' tail1 >tail1_tmp
    sed 's/ C33 DOPA/ C13 OL  /' tail1_tmp >tail1
    sed 's/ C34 DOPA/ C14 OL  /' tail1 >tail1_tmp
    sed 's/ C35 DOPA/ C15 OL  /' tail1_tmp >tail1
    sed 's/ C36 DOPA/ C16 OL  /' tail1 >tail1_tmp
    sed 's/ C37 DOPA/ C17 OL  /' tail1_tmp >tail1
    sed 's/ C38 DOPA/ C18 OL  /' tail1 >tail1_tmp
    sed 's/ C39 DOPA/ C19 OL  /' tail1_tmp >tail1
    sed 's/C310 DOPA/C110 OL  /' tail1 >tail1_tmp
    sed 's/C311 DOPA/C111 OL  /' tail1_tmp >tail1
    sed 's/C312 DOPA/C112 OL  /' tail1 >tail1_tmp
    sed 's/C313 DOPA/C113 OL  /' tail1_tmp >tail1
    sed 's/C314 DOPA/C114 OL  /' tail1 >tail1_tmp
    sed 's/C315 DOPA/C115 OL  /' tail1_tmp >tail1
    sed 's/C316 DOPA/C116 OL  /' tail1 >tail1_tmp
    sed 's/C317 DOPA/C117 OL  /' tail1_tmp >tail1
    sed 's/C318 DOPA/C118 OL  /' tail1 >tail1_tmp
    sed 's/ H2X DOPA/ H2R OL  /' tail1_tmp >tail1
    sed 's/ H2Y DOPA/ H2S OL  /' tail1 >tail1_tmp
    sed 's/ H3X DOPA/ H3R OL  /' tail1_tmp >tail1
    sed 's/ H3Y DOPA/ H3S OL  /' tail1 >tail1_tmp
    sed 's/ H4X DOPA/ H4R OL  /' tail1_tmp >tail1
    sed 's/ H4Y DOPA/ H4S OL  /' tail1 >tail1_tmp
    sed 's/ H5X DOPA/ H5R OL  /' tail1_tmp >tail1
    sed 's/ H5Y DOPA/ H5S OL  /' tail1 >tail1_tmp
    sed 's/ H6X DOPA/ H6R OL  /' tail1_tmp >tail1
    sed 's/ H6Y DOPA/ H6S OL  /' tail1 >tail1_tmp
    sed 's/ H7X DOPA/ H7R OL  /' tail1_tmp >tail1
    sed 's/ H7Y DOPA/ H7S OL  /' tail1 >tail1_tmp
    sed 's/ H8X DOPA/ H8R OL  /' tail1_tmp >tail1
    sed 's/ H8Y DOPA/ H8S OL  /' tail1 >tail1_tmp
    sed 's/ H9X DOPA/ H9R OL  /' tail1_tmp >tail1
    sed 's/H10X DOPA/H10R OL  /' tail1 >tail1_tmp
    sed 's/H11X DOPA/H11R OL  /' tail1_tmp >tail1
    sed 's/H11Y DOPA/H11S OL  /' tail1 >tail1_tmp
    sed 's/H12X DOPA/H12R OL  /' tail1_tmp >tail1
    sed 's/H12Y DOPA/H12S OL  /' tail1 >tail1_tmp
    sed 's/H13X DOPA/H13R OL  /' tail1_tmp >tail1
    sed 's/H13Y DOPA/H13S OL  /' tail1 >tail1_tmp
    sed 's/H14X DOPA/H14R OL  /' tail1_tmp >tail1
    sed 's/H14Y DOPA/H14S OL  /' tail1 >tail1_tmp
    sed 's/H15X DOPA/H15R OL  /' tail1_tmp >tail1
    sed 's/H15Y DOPA/H15S OL  /' tail1 >tail1_tmp
    sed 's/H16X DOPA/H16R OL  /' tail1_tmp >tail1
    sed 's/H16Y DOPA/H16S OL  /' tail1 >tail1_tmp
    sed 's/H17X DOPA/H17R OL  /' tail1_tmp >tail1
    sed 's/H17Y DOPA/H17S OL  /' tail1 >tail1_tmp
    sed 's/H18X DOPA/H18R OL  /' tail1_tmp >tail1
    sed 's/H18Y DOPA/H18S OL  /' tail1 >tail1_tmp
    sed 's/H18Z DOPA/H18T OL  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DOPA/ C12 OL  /' tail2 >tail2_tmp
    sed 's/ C23 DOPA/ C13 OL  /' tail2_tmp >tail2
    sed 's/ C24 DOPA/ C14 OL  /' tail2 >tail2_tmp
    sed 's/ C25 DOPA/ C15 OL  /' tail2_tmp >tail2
    sed 's/ C26 DOPA/ C16 OL  /' tail2 >tail2_tmp
    sed 's/ C27 DOPA/ C17 OL  /' tail2_tmp >tail2
    sed 's/ C28 DOPA/ C18 OL  /' tail2 >tail2_tmp
    sed 's/ C29 DOPA/ C19 OL  /' tail2_tmp >tail2
    sed 's/C210 DOPA/C110 OL  /' tail2 >tail2_tmp
    sed 's/C211 DOPA/C111 OL  /' tail2_tmp >tail2
    sed 's/C212 DOPA/C112 OL  /' tail2 >tail2_tmp
    sed 's/C213 DOPA/C113 OL  /' tail2_tmp >tail2
    sed 's/C214 DOPA/C114 OL  /' tail2 >tail2_tmp
    sed 's/C215 DOPA/C115 OL  /' tail2_tmp >tail2
    sed 's/C216 DOPA/C116 OL  /' tail2 >tail2_tmp
    sed 's/C217 DOPA/C117 OL  /' tail2_tmp >tail2
    sed 's/C218 DOPA/C118 OL  /' tail2 >tail2_tmp
    sed 's/ H2R DOPA/ H2R OL  /' tail2_tmp >tail2
    sed 's/ H2S DOPA/ H2S OL  /' tail2 >tail2_tmp
    sed 's/ H3R DOPA/ H3R OL  /' tail2_tmp >tail2
    sed 's/ H3S DOPA/ H3S OL  /' tail2 >tail2_tmp
    sed 's/ H4R DOPA/ H4R OL  /' tail2_tmp >tail2
    sed 's/ H4S DOPA/ H4S OL  /' tail2 >tail2_tmp
    sed 's/ H5R DOPA/ H5R OL  /' tail2_tmp >tail2
    sed 's/ H5S DOPA/ H5S OL  /' tail2 >tail2_tmp
    sed 's/ H6R DOPA/ H6R OL  /' tail2_tmp >tail2
    sed 's/ H6S DOPA/ H6S OL  /' tail2 >tail2_tmp
    sed 's/ H7R DOPA/ H7R OL  /' tail2_tmp >tail2
    sed 's/ H7S DOPA/ H7S OL  /' tail2 >tail2_tmp
    sed 's/ H8R DOPA/ H8R OL  /' tail2_tmp >tail2
    sed 's/ H8S DOPA/ H8S OL  /' tail2 >tail2_tmp
    sed 's/ H9R DOPA/ H9R OL  /' tail2_tmp >tail2
    sed 's/H10R DOPA/H10R OL  /' tail2 >tail2_tmp
    sed 's/H11R DOPA/H11R OL  /' tail2_tmp >tail2
    sed 's/H11S DOPA/H11S OL  /' tail2 >tail2_tmp
    sed 's/H12R DOPA/H12R OL  /' tail2_tmp >tail2
    sed 's/H12S DOPA/H12S OL  /' tail2 >tail2_tmp
    sed 's/H13R DOPA/H13R OL  /' tail2_tmp >tail2
    sed 's/H13S DOPA/H13S OL  /' tail2 >tail2_tmp
    sed 's/H14R DOPA/H14R OL  /' tail2_tmp >tail2
    sed 's/H14S DOPA/H14S OL  /' tail2 >tail2_tmp
    sed 's/H15R DOPA/H15R OL  /' tail2_tmp >tail2
    sed 's/H15S DOPA/H15S OL  /' tail2 >tail2_tmp
    sed 's/H16R DOPA/H16R OL  /' tail2_tmp >tail2
    sed 's/H16S DOPA/H16S OL  /' tail2 >tail2_tmp
    sed 's/H17R DOPA/H17R OL  /' tail2_tmp >tail2
    sed 's/H17S DOPA/H17S OL  /' tail2 >tail2_tmp
    sed 's/H18R DOPA/H18R OL  /' tail2_tmp >tail2
    sed 's/H18S DOPA/H18S OL  /' tail2 >tail2_tmp
    sed 's/H18T DOPA/H18T OL  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ P   DOPA/ P31 PH- /' head > head_tmp
    sed 's/ O13 DOPA/ O33 PH- /' head_tmp > head
    sed 's/ O14 DOPA/ O34 PH- /' head > head_tmp
    sed 's/ O12 DOPA/ O32 PH- /' head_tmp > head
    sed 's/ H12 DOPA/HO2A PH- /' head > head_tmp
    sed 's/ O11 DOPA/ O31 PH- /' head_tmp > head
    sed 's/ C1  DOPA/ C3  PH- /' head > head_tmp
    sed 's/ HA  DOPA/ HA  PH- /' head_tmp > head
    sed 's/ HB  DOPA/ HB  PH- /' head > head_tmp
    sed 's/ C2  DOPA/ C2  PH- /' head_tmp > head
    sed 's/ HS  DOPA/ HX  PH- /' head > head_tmp
    sed 's/ O21 DOPA/ O21 PH- /' head_tmp > head
    sed 's/ C21 DOPA/ C21 PH- /' head > head_tmp
    sed 's/ O22 DOPA/ O22 PH- /' head_tmp > head
    sed 's/ C3  DOPA/ C1  PH- /' head > head_tmp
    sed 's/ HX  DOPA/ HR  PH- /' head_tmp > head
    sed 's/ HY  DOPA/ HS  PH- /' head > head_tmp
    sed 's/ O31 DOPA/ O11 PH- /' head_tmp > head
    sed 's/ C31 DOPA/ C11 PH- /' head > head_tmp
    sed 's/ O32 DOPA/ O12 PH- /' head_tmp > head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#POPC
if ( $popc_found == 1 ) then
  #POPC is 134 atoms per lipid. 
  set atoms_per_lipid = 134

  echo " Processing: POPC"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l POPC_tmp | awk '{print$1}'`
  echo "  POPC Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted POPC Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " POPC Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../POPC_tmp
  rm -f ../POPC_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing POPC Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 POPC" $i > tail1
    grep -A42 "C33 POPC" $i >> tail1

    grep -A2 "C22 POPC" $i > tail2
    grep -A46 "C23 POPC" $i >> tail2

    head -n32 $i > head
    grep -A5 "C3  POPC" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 POPC/ C12 PA  /' tail1 >tail1_tmp
    sed 's/ C33 POPC/ C13 PA  /' tail1_tmp >tail1
    sed 's/ C34 POPC/ C14 PA  /' tail1 >tail1_tmp
    sed 's/ C35 POPC/ C15 PA  /' tail1_tmp >tail1
    sed 's/ C36 POPC/ C16 PA  /' tail1 >tail1_tmp
    sed 's/ C37 POPC/ C17 PA  /' tail1_tmp >tail1
    sed 's/ C38 POPC/ C18 PA  /' tail1 >tail1_tmp
    sed 's/ C39 POPC/ C19 PA  /' tail1_tmp >tail1
    sed 's/C310 POPC/C110 PA  /' tail1 >tail1_tmp
    sed 's/C311 POPC/C111 PA  /' tail1_tmp >tail1
    sed 's/C312 POPC/C112 PA  /' tail1 >tail1_tmp
    sed 's/C313 POPC/C113 PA  /' tail1_tmp >tail1
    sed 's/C314 POPC/C114 PA  /' tail1 >tail1_tmp
    sed 's/C315 POPC/C115 PA  /' tail1_tmp >tail1
    sed 's/C316 POPC/C116 PA  /' tail1 >tail1_tmp
    sed 's/C317 POPC/C117 PA  /' tail1_tmp >tail1
    sed 's/C318 POPC/C118 PA  /' tail1 >tail1_tmp
    sed 's/ H2X POPC/ H2R PA  /' tail1_tmp >tail1
    sed 's/ H2Y POPC/ H2S PA  /' tail1 >tail1_tmp
    sed 's/ H3X POPC/ H3R PA  /' tail1_tmp >tail1
    sed 's/ H3Y POPC/ H3S PA  /' tail1 >tail1_tmp
    sed 's/ H4X POPC/ H4R PA  /' tail1_tmp >tail1
    sed 's/ H4Y POPC/ H4S PA  /' tail1 >tail1_tmp
    sed 's/ H5X POPC/ H5R PA  /' tail1_tmp >tail1
    sed 's/ H5Y POPC/ H5S PA  /' tail1 >tail1_tmp
    sed 's/ H6X POPC/ H6R PA  /' tail1_tmp >tail1
    sed 's/ H6Y POPC/ H6S PA  /' tail1 >tail1_tmp
    sed 's/ H7X POPC/ H7R PA  /' tail1_tmp >tail1
    sed 's/ H7Y POPC/ H7S PA  /' tail1 >tail1_tmp
    sed 's/ H8X POPC/ H8R PA  /' tail1_tmp >tail1
    sed 's/ H8Y POPC/ H8S PA  /' tail1 >tail1_tmp
    sed 's/ H9X POPC/ H9R PA  /' tail1_tmp >tail1
    sed 's/ H9Y POPC/ H9S PA  /' tail1 >tail1_tmp
    sed 's/H10X POPC/H10R PA  /' tail1_tmp >tail1
    sed 's/H10Y POPC/H10S PA  /' tail1 >tail1_tmp
    sed 's/H11X POPC/H11R PA  /' tail1_tmp >tail1
    sed 's/H11Y POPC/H11S PA  /' tail1 >tail1_tmp
    sed 's/H12X POPC/H12R PA  /' tail1_tmp >tail1
    sed 's/H12Y POPC/H12S PA  /' tail1 >tail1_tmp
    sed 's/H13X POPC/H13R PA  /' tail1_tmp >tail1
    sed 's/H13Y POPC/H13S PA  /' tail1 >tail1_tmp
    sed 's/H14X POPC/H14R PA  /' tail1_tmp >tail1
    sed 's/H14Y POPC/H14S PA  /' tail1 >tail1_tmp
    sed 's/H15X POPC/H15R PA  /' tail1_tmp >tail1
    sed 's/H15Y POPC/H15S PA  /' tail1 >tail1_tmp
    sed 's/H16X POPC/H16R PA  /' tail1_tmp >tail1
    sed 's/H16Y POPC/H16S PA  /' tail1 >tail1_tmp
    sed 's/H16Z POPC/H16T PA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 POPC/ C12 OL  /' tail2 >tail2_tmp
    sed 's/ C23 POPC/ C13 OL  /' tail2_tmp >tail2
    sed 's/ C24 POPC/ C14 OL  /' tail2 >tail2_tmp
    sed 's/ C25 POPC/ C15 OL  /' tail2_tmp >tail2
    sed 's/ C26 POPC/ C16 OL  /' tail2 >tail2_tmp
    sed 's/ C27 POPC/ C17 OL  /' tail2_tmp >tail2
    sed 's/ C28 POPC/ C18 OL  /' tail2 >tail2_tmp
    sed 's/ C29 POPC/ C19 OL  /' tail2_tmp >tail2
    sed 's/C210 POPC/C110 OL  /' tail2 >tail2_tmp
    sed 's/C211 POPC/C111 OL  /' tail2_tmp >tail2
    sed 's/C212 POPC/C112 OL  /' tail2 >tail2_tmp
    sed 's/C213 POPC/C113 OL  /' tail2_tmp >tail2
    sed 's/C214 POPC/C114 OL  /' tail2 >tail2_tmp
    sed 's/C215 POPC/C115 OL  /' tail2_tmp >tail2
    sed 's/C216 POPC/C116 OL  /' tail2 >tail2_tmp
    sed 's/C217 POPC/C117 OL  /' tail2_tmp >tail2
    sed 's/C218 POPC/C118 OL  /' tail2 >tail2_tmp
    sed 's/ H2R POPC/ H2R OL  /' tail2_tmp >tail2
    sed 's/ H2S POPC/ H2S OL  /' tail2 >tail2_tmp
    sed 's/ H3R POPC/ H3R OL  /' tail2_tmp >tail2
    sed 's/ H3S POPC/ H3S OL  /' tail2 >tail2_tmp
    sed 's/ H4R POPC/ H4R OL  /' tail2_tmp >tail2
    sed 's/ H4S POPC/ H4S OL  /' tail2 >tail2_tmp
    sed 's/ H5R POPC/ H5R OL  /' tail2_tmp >tail2
    sed 's/ H5S POPC/ H5S OL  /' tail2 >tail2_tmp
    sed 's/ H6R POPC/ H6R OL  /' tail2_tmp >tail2
    sed 's/ H6S POPC/ H6S OL  /' tail2 >tail2_tmp
    sed 's/ H7R POPC/ H7R OL  /' tail2_tmp >tail2
    sed 's/ H7S POPC/ H7S OL  /' tail2 >tail2_tmp
    sed 's/ H8R POPC/ H8R OL  /' tail2_tmp >tail2
    sed 's/ H8S POPC/ H8S OL  /' tail2 >tail2_tmp
    sed 's/ H91 POPC/ H9R OL  /' tail2_tmp >tail2
    sed 's/H101 POPC/H10R OL  /' tail2 >tail2_tmp
    sed 's/H11R POPC/H11R OL  /' tail2_tmp >tail2
    sed 's/H11S POPC/H11S OL  /' tail2 >tail2_tmp
    sed 's/H12R POPC/H12R OL  /' tail2_tmp >tail2
    sed 's/H12S POPC/H12S OL  /' tail2 >tail2_tmp
    sed 's/H13R POPC/H13R OL  /' tail2_tmp >tail2
    sed 's/H13S POPC/H13S OL  /' tail2 >tail2_tmp
    sed 's/H14R POPC/H14R OL  /' tail2_tmp >tail2
    sed 's/H14S POPC/H14S OL  /' tail2 >tail2_tmp
    sed 's/H15R POPC/H15R OL  /' tail2_tmp >tail2
    sed 's/H15S POPC/H15S OL  /' tail2 >tail2_tmp
    sed 's/H16R POPC/H16R OL  /' tail2_tmp >tail2
    sed 's/H16S POPC/H16S OL  /' tail2 >tail2_tmp
    sed 's/H17R POPC/H17R OL  /' tail2_tmp >tail2
    sed 's/H17S POPC/H17S OL  /' tail2 >tail2_tmp
    sed 's/H18R POPC/H18R OL  /' tail2_tmp >tail2
    sed 's/H18S POPC/H18S OL  /' tail2 >tail2_tmp
    sed 's/H18T POPC/H18T OL  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   POPC/ N31 PC  /' head > head_tmp
    sed 's/ C13 POPC/ C33 PC  /' head_tmp > head
    sed 's/ C14 POPC/ C34 PC  /' head > head_tmp
    sed 's/ C15 POPC/ C35 PC  /' head_tmp > head
    sed 's/ C12 POPC/ C32 PC  /' head > head_tmp
    sed 's/ C11 POPC/ C31 PC  /' head_tmp > head
    sed 's/ O12 POPC/ O32 PC  /' head > head_tmp
    sed 's/ P   POPC/ P31 PC  /' head_tmp > head
    sed 's/ O11 POPC/ O31 PC  /' head > head_tmp
    sed 's/ O13 POPC/ O33 PC  /' head_tmp > head
    sed 's/ O14 POPC/ O34 PC  /' head > head_tmp
    sed 's/ C1  POPC/ C3  PC  /' head_tmp > head
    sed 's/ C2  POPC/ C2  PC  /' head > head_tmp
    sed 's/ C3  POPC/ C1  PC  /' head_tmp > head
    sed 's/ O31 POPC/ O11 PC  /' head > head_tmp
    sed 's/ C31 POPC/ C11 PC  /' head_tmp > head
    sed 's/ O32 POPC/ O12 PC  /' head > head_tmp
    sed 's/ O21 POPC/ O21 PC  /' head_tmp > head
    sed 's/ C21 POPC/ C21 PC  /' head > head_tmp
    sed 's/ O22 POPC/ O22 PC  /' head_tmp > head
    sed 's/ HX  POPC/ HR  PC  /' head > head_tmp
    sed 's/ HY  POPC/ HS  PC  /' head_tmp > head
    sed 's/ HS  POPC/ HX  PC  /' head > head_tmp
    sed 's/ HA  POPC/ HA  PC  /' head_tmp > head
    sed 's/ HB  POPC/ HB  PC  /' head > head_tmp
    sed 's/H11A POPC/ H1A PC  /' head_tmp > head
    sed 's/H11B POPC/ H1B PC  /' head > head_tmp
    sed 's/H12A POPC/ H2A PC  /' head_tmp > head
    sed 's/H12B POPC/ H2B PC  /' head > head_tmp
    sed 's/H13A POPC/ H3A PC  /' head_tmp > head
    sed 's/H13B POPC/ H3B PC  /' head > head_tmp
    sed 's/H13C POPC/ H3C PC  /' head_tmp > head
    sed 's/H14A POPC/ H4A PC  /' head > head_tmp
    sed 's/H14B POPC/ H4B PC  /' head_tmp > head
    sed 's/H14C POPC/ H4C PC  /' head > head_tmp
    sed 's/H15A POPC/ H5A PC  /' head_tmp > head
    sed 's/H15B POPC/ H5B PC  /' head > head_tmp
    sed 's/H15C POPC/ H5C PC  /' head_tmp > head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#POPE
if ( $pope_found == 1 ) then
  #POPE is 125 atoms per lipid. 
  set atoms_per_lipid = 125

  echo " Processing: POPE"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l POPE_tmp | awk '{print$1}'`
  echo "  POPE Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted POPE Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " POPE Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../POPE_tmp
  rm -f ../POPE_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing POPE Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 POPE" $i > tail1
    grep -A42 "C33 POPE" $i >> tail1

    grep -A2 "C22 POPE" $i > tail2
    grep -A46 "C23 POPE" $i >> tail2

    head -n23 $i > head
    grep -A5 "C3  POPE" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 POPE/ C12 PA  /' tail1 >tail1_tmp
    sed 's/ C33 POPE/ C13 PA  /' tail1_tmp >tail1
    sed 's/ C34 POPE/ C14 PA  /' tail1 >tail1_tmp
    sed 's/ C35 POPE/ C15 PA  /' tail1_tmp >tail1
    sed 's/ C36 POPE/ C16 PA  /' tail1 >tail1_tmp
    sed 's/ C37 POPE/ C17 PA  /' tail1_tmp >tail1
    sed 's/ C38 POPE/ C18 PA  /' tail1 >tail1_tmp
    sed 's/ C39 POPE/ C19 PA  /' tail1_tmp >tail1
    sed 's/C310 POPE/C110 PA  /' tail1 >tail1_tmp
    sed 's/C311 POPE/C111 PA  /' tail1_tmp >tail1
    sed 's/C312 POPE/C112 PA  /' tail1 >tail1_tmp
    sed 's/C313 POPE/C113 PA  /' tail1_tmp >tail1
    sed 's/C314 POPE/C114 PA  /' tail1 >tail1_tmp
    sed 's/C315 POPE/C115 PA  /' tail1_tmp >tail1
    sed 's/C316 POPE/C116 PA  /' tail1 >tail1_tmp
    sed 's/ H2X POPE/ H2R PA  /' tail1_tmp >tail1
    sed 's/ H2Y POPE/ H2S PA  /' tail1 >tail1_tmp
    sed 's/ H3X POPE/ H3R PA  /' tail1_tmp >tail1
    sed 's/ H3Y POPE/ H3S PA  /' tail1 >tail1_tmp
    sed 's/ H4X POPE/ H4R PA  /' tail1_tmp >tail1
    sed 's/ H4Y POPE/ H4S PA  /' tail1 >tail1_tmp
    sed 's/ H5X POPE/ H5R PA  /' tail1_tmp >tail1
    sed 's/ H5Y POPE/ H5S PA  /' tail1 >tail1_tmp
    sed 's/ H6X POPE/ H6R PA  /' tail1_tmp >tail1
    sed 's/ H6Y POPE/ H6S PA  /' tail1 >tail1_tmp
    sed 's/ H7X POPE/ H7R PA  /' tail1_tmp >tail1
    sed 's/ H7Y POPE/ H7S PA  /' tail1 >tail1_tmp
    sed 's/ H8X POPE/ H8R PA  /' tail1_tmp >tail1
    sed 's/ H8Y POPE/ H8S PA  /' tail1 >tail1_tmp
    sed 's/ H9X POPE/ H9R PA  /' tail1_tmp >tail1
    sed 's/ H9Y POPE/ H9S PA  /' tail1 >tail1_tmp
    sed 's/H10X POPE/H10R PA  /' tail1_tmp >tail1
    sed 's/H10Y POPE/H10S PA  /' tail1 >tail1_tmp
    sed 's/H11X POPE/H11R PA  /' tail1_tmp >tail1
    sed 's/H11Y POPE/H11S PA  /' tail1 >tail1_tmp
    sed 's/H12X POPE/H12R PA  /' tail1_tmp >tail1
    sed 's/H12Y POPE/H12S PA  /' tail1 >tail1_tmp
    sed 's/H13X POPE/H13R PA  /' tail1_tmp >tail1
    sed 's/H13Y POPE/H13S PA  /' tail1 >tail1_tmp
    sed 's/H14X POPE/H14R PA  /' tail1_tmp >tail1
    sed 's/H14Y POPE/H14S PA  /' tail1 >tail1_tmp
    sed 's/H15X POPE/H15R PA  /' tail1_tmp >tail1
    sed 's/H15Y POPE/H15S PA  /' tail1 >tail1_tmp
    sed 's/H16X POPE/H16R PA  /' tail1_tmp >tail1
    sed 's/H16Y POPE/H16S PA  /' tail1 >tail1_tmp
    sed 's/H16Z POPE/H16T PA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 POPE/ C12 OL  /' tail2 >tail2_tmp
    sed 's/ C23 POPE/ C13 OL  /' tail2_tmp >tail2
    sed 's/ C24 POPE/ C14 OL  /' tail2 >tail2_tmp
    sed 's/ C25 POPE/ C15 OL  /' tail2_tmp >tail2
    sed 's/ C26 POPE/ C16 OL  /' tail2 >tail2_tmp
    sed 's/ C27 POPE/ C17 OL  /' tail2_tmp >tail2
    sed 's/ C28 POPE/ C18 OL  /' tail2 >tail2_tmp
    sed 's/ C29 POPE/ C19 OL  /' tail2_tmp >tail2
    sed 's/C210 POPE/C110 OL  /' tail2 >tail2_tmp
    sed 's/C211 POPE/C111 OL  /' tail2_tmp >tail2
    sed 's/C212 POPE/C112 OL  /' tail2 >tail2_tmp
    sed 's/C213 POPE/C113 OL  /' tail2_tmp >tail2
    sed 's/C214 POPE/C114 OL  /' tail2 >tail2_tmp
    sed 's/C215 POPE/C115 OL  /' tail2_tmp >tail2
    sed 's/C216 POPE/C116 OL  /' tail2 >tail2_tmp
    sed 's/C217 POPE/C117 OL  /' tail2_tmp >tail2
    sed 's/C218 POPE/C118 OL  /' tail2 >tail2_tmp
    sed 's/ H2R POPE/ H2R OL  /' tail2_tmp >tail2
    sed 's/ H2S POPE/ H2S OL  /' tail2 >tail2_tmp
    sed 's/ H3R POPE/ H3R OL  /' tail2_tmp >tail2
    sed 's/ H3S POPE/ H3S OL  /' tail2 >tail2_tmp
    sed 's/ H4R POPE/ H4R OL  /' tail2_tmp >tail2
    sed 's/ H4S POPE/ H4S OL  /' tail2 >tail2_tmp
    sed 's/ H5R POPE/ H5R OL  /' tail2_tmp >tail2
    sed 's/ H5S POPE/ H5S OL  /' tail2 >tail2_tmp
    sed 's/ H6R POPE/ H6R OL  /' tail2_tmp >tail2
    sed 's/ H6S POPE/ H6S OL  /' tail2 >tail2_tmp
    sed 's/ H7R POPE/ H7R OL  /' tail2_tmp >tail2
    sed 's/ H7S POPE/ H7S OL  /' tail2 >tail2_tmp
    sed 's/ H8R POPE/ H8R OL  /' tail2_tmp >tail2
    sed 's/ H8S POPE/ H8S OL  /' tail2 >tail2_tmp
    sed 's/ H91 POPE/ H9R OL  /' tail2_tmp >tail2
    sed 's/H101 POPE/H10R OL  /' tail2 >tail2_tmp
    sed 's/H11R POPE/H11R OL  /' tail2_tmp >tail2
    sed 's/H11S POPE/H11S OL  /' tail2 >tail2_tmp
    sed 's/H12R POPE/H12R OL  /' tail2_tmp >tail2
    sed 's/H12S POPE/H12S OL  /' tail2 >tail2_tmp
    sed 's/H13R POPE/H13R OL  /' tail2_tmp >tail2
    sed 's/H13S POPE/H13S OL  /' tail2 >tail2_tmp
    sed 's/H14R POPE/H14R OL  /' tail2_tmp >tail2
    sed 's/H14S POPE/H14S OL  /' tail2 >tail2_tmp
    sed 's/H15R POPE/H15R OL  /' tail2_tmp >tail2
    sed 's/H15S POPE/H15S OL  /' tail2 >tail2_tmp
    sed 's/H16R POPE/H16R OL  /' tail2_tmp >tail2
    sed 's/H16S POPE/H16S OL  /' tail2 >tail2_tmp
    sed 's/H17R POPE/H17R OL  /' tail2_tmp >tail2
    sed 's/H17S POPE/H17S OL  /' tail2 >tail2_tmp
    sed 's/H18R POPE/H18R OL  /' tail2_tmp >tail2
    sed 's/H18S POPE/H18S OL  /' tail2 >tail2_tmp
    sed 's/H18T POPE/H18T OL  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   POPE/ N31 PE  /' head > head_tmp
    sed 's/ C12 POPE/ C32 PE  /' head_tmp > head
    sed 's/ C11 POPE/ C31 PE  /' head > head_tmp
    sed 's/ O12 POPE/ O32 PE  /' head_tmp > head
    sed 's/ P   POPE/ P31 PE  /' head > head_tmp
    sed 's/ O11 POPE/ O31 PE  /' head_tmp > head
    sed 's/ O13 POPE/ O33 PE  /' head > head_tmp
    sed 's/ O14 POPE/ O34 PE  /' head_tmp > head
    sed 's/ C1  POPE/ C3  PE  /' head > head_tmp
    sed 's/ C2  POPE/ C2  PE  /' head_tmp > head
    sed 's/ C3  POPE/ C1  PE  /' head > head_tmp
    sed 's/ O31 POPE/ O11 PE  /' head_tmp > head
    sed 's/ C31 POPE/ C11 PE  /' head > head_tmp
    sed 's/ O32 POPE/ O12 PE  /' head_tmp > head
    sed 's/ O21 POPE/ O21 PE  /' head > head_tmp
    sed 's/ C21 POPE/ C21 PE  /' head_tmp > head
    sed 's/ O22 POPE/ O22 PE  /' head > head_tmp
    sed 's/ HX  POPE/ HR  PE  /' head_tmp > head
    sed 's/ HY  POPE/ HS  PE  /' head > head_tmp
    sed 's/ HS  POPE/ HX  PE  /' head_tmp > head
    sed 's/ HA  POPE/ HA  PE  /' head > head_tmp
    sed 's/ HB  POPE/ HB  PE  /' head_tmp > head
    sed 's/H11A POPE/ H1A PE  /' head > head_tmp
    sed 's/H11B POPE/ H1B PE  /' head_tmp > head
    sed 's/H12A POPE/ H2A PE  /' head > head_tmp
    sed 's/H12B POPE/ H2B PE  /' head_tmp > head
    sed 's/ HN1 POPE/HN1A PE  /' head > head_tmp
    sed 's/ HN2 POPE/HN1B PE  /' head_tmp > head
    sed 's/ HN3 POPE/HN1C PE  /' head > head_tmp
    mv head_tmp head  
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif


#POPS
if ( $pops_found == 1 ) then
  #POPS is 127 atoms per lipid. 
  set atoms_per_lipid = 127

  echo " Processing: POPS"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l POPS_tmp | awk '{print$1}'`
  echo "  POPS Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted POPS Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " POPS Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../POPS_tmp
  rm -f ../POPS_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing POPS Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 POPS" $i > tail1
    grep -A42 "C33 POPS" $i >> tail1

    grep -A2 "C22 POPS" $i > tail2
    grep -A46 "C23 POPS" $i >> tail2

    head -n25 $i > head
    grep -A5 "C3  POPS" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 POPS/ C12 PA  /' tail1 >tail1_tmp
    sed 's/ C33 POPS/ C13 PA  /' tail1_tmp >tail1
    sed 's/ C34 POPS/ C14 PA  /' tail1 >tail1_tmp
    sed 's/ C35 POPS/ C15 PA  /' tail1_tmp >tail1
    sed 's/ C36 POPS/ C16 PA  /' tail1 >tail1_tmp
    sed 's/ C37 POPS/ C17 PA  /' tail1_tmp >tail1
    sed 's/ C38 POPS/ C18 PA  /' tail1 >tail1_tmp
    sed 's/ C39 POPS/ C19 PA  /' tail1_tmp >tail1
    sed 's/C310 POPS/C110 PA  /' tail1 >tail1_tmp
    sed 's/C311 POPS/C111 PA  /' tail1_tmp >tail1
    sed 's/C312 POPS/C112 PA  /' tail1 >tail1_tmp
    sed 's/C313 POPS/C113 PA  /' tail1_tmp >tail1
    sed 's/C314 POPS/C114 PA  /' tail1 >tail1_tmp
    sed 's/C315 POPS/C115 PA  /' tail1_tmp >tail1
    sed 's/C316 POPS/C116 PA  /' tail1 >tail1_tmp
    sed 's/C317 POPS/C117 PA  /' tail1_tmp >tail1
    sed 's/C318 POPS/C118 PA  /' tail1 >tail1_tmp
    sed 's/ H2X POPS/ H2R PA  /' tail1_tmp >tail1
    sed 's/ H2Y POPS/ H2S PA  /' tail1 >tail1_tmp
    sed 's/ H3X POPS/ H3R PA  /' tail1_tmp >tail1
    sed 's/ H3Y POPS/ H3S PA  /' tail1 >tail1_tmp
    sed 's/ H4X POPS/ H4R PA  /' tail1_tmp >tail1
    sed 's/ H4Y POPS/ H4S PA  /' tail1 >tail1_tmp
    sed 's/ H5X POPS/ H5R PA  /' tail1_tmp >tail1
    sed 's/ H5Y POPS/ H5S PA  /' tail1 >tail1_tmp
    sed 's/ H6X POPS/ H6R PA  /' tail1_tmp >tail1
    sed 's/ H6Y POPS/ H6S PA  /' tail1 >tail1_tmp
    sed 's/ H7X POPS/ H7R PA  /' tail1_tmp >tail1
    sed 's/ H7Y POPS/ H7S PA  /' tail1 >tail1_tmp
    sed 's/ H8X POPS/ H8R PA  /' tail1_tmp >tail1
    sed 's/ H8Y POPS/ H8S PA  /' tail1 >tail1_tmp
    sed 's/ H9X POPS/ H9R PA  /' tail1_tmp >tail1
    sed 's/ H9Y POPS/ H9S PA  /' tail1 >tail1_tmp
    sed 's/H10X POPS/H10R PA  /' tail1_tmp >tail1
    sed 's/H10Y POPS/H10S PA  /' tail1 >tail1_tmp
    sed 's/H11X POPS/H11R PA  /' tail1_tmp >tail1
    sed 's/H11Y POPS/H11S PA  /' tail1 >tail1_tmp
    sed 's/H12X POPS/H12R PA  /' tail1_tmp >tail1
    sed 's/H12Y POPS/H12S PA  /' tail1 >tail1_tmp
    sed 's/H13X POPS/H13R PA  /' tail1_tmp >tail1
    sed 's/H13Y POPS/H13S PA  /' tail1 >tail1_tmp
    sed 's/H14X POPS/H14R PA  /' tail1_tmp >tail1
    sed 's/H14Y POPS/H14S PA  /' tail1 >tail1_tmp
    sed 's/H15X POPS/H15R PA  /' tail1_tmp >tail1
    sed 's/H15Y POPS/H15S PA  /' tail1 >tail1_tmp
    sed 's/H16X POPS/H16R PA  /' tail1_tmp >tail1
    sed 's/H16Y POPS/H16S PA  /' tail1 >tail1_tmp
    sed 's/H16Z POPS/H16T PA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 POPS/ C12 OL  /' tail2 >tail2_tmp
    sed 's/ C23 POPS/ C13 OL  /' tail2_tmp >tail2
    sed 's/ C24 POPS/ C14 OL  /' tail2 >tail2_tmp
    sed 's/ C25 POPS/ C15 OL  /' tail2_tmp >tail2
    sed 's/ C26 POPS/ C16 OL  /' tail2 >tail2_tmp
    sed 's/ C27 POPS/ C17 OL  /' tail2_tmp >tail2
    sed 's/ C28 POPS/ C18 OL  /' tail2 >tail2_tmp
    sed 's/ C29 POPS/ C19 OL  /' tail2_tmp >tail2
    sed 's/C210 POPS/C110 OL  /' tail2 >tail2_tmp
    sed 's/C211 POPS/C111 OL  /' tail2_tmp >tail2
    sed 's/C212 POPS/C112 OL  /' tail2 >tail2_tmp
    sed 's/C213 POPS/C113 OL  /' tail2_tmp >tail2
    sed 's/C214 POPS/C114 OL  /' tail2 >tail2_tmp
    sed 's/C215 POPS/C115 OL  /' tail2_tmp >tail2
    sed 's/C216 POPS/C116 OL  /' tail2 >tail2_tmp
    sed 's/C217 POPS/C117 OL  /' tail2_tmp >tail2
    sed 's/C218 POPS/C118 OL  /' tail2 >tail2_tmp
    sed 's/ H2R POPS/ H2R OL  /' tail2_tmp >tail2
    sed 's/ H2S POPS/ H2S OL  /' tail2 >tail2_tmp
    sed 's/ H3R POPS/ H3R OL  /' tail2_tmp >tail2
    sed 's/ H3S POPS/ H3S OL  /' tail2 >tail2_tmp
    sed 's/ H4R POPS/ H4R OL  /' tail2_tmp >tail2
    sed 's/ H4S POPS/ H4S OL  /' tail2 >tail2_tmp
    sed 's/ H5R POPS/ H5R OL  /' tail2_tmp >tail2
    sed 's/ H5S POPS/ H5S OL  /' tail2 >tail2_tmp
    sed 's/ H6R POPS/ H6R OL  /' tail2_tmp >tail2
    sed 's/ H6S POPS/ H6S OL  /' tail2 >tail2_tmp
    sed 's/ H7R POPS/ H7R OL  /' tail2_tmp >tail2
    sed 's/ H7S POPS/ H7S OL  /' tail2 >tail2_tmp
    sed 's/ H8R POPS/ H8R OL  /' tail2_tmp >tail2
    sed 's/ H8S POPS/ H8S OL  /' tail2 >tail2_tmp
    sed 's/ H91 POPS/ H9R OL  /' tail2_tmp >tail2
    sed 's/H101 POPS/H10R OL  /' tail2 >tail2_tmp
    sed 's/H11R POPS/H11R OL  /' tail2_tmp >tail2
    sed 's/H11S POPS/H11S OL  /' tail2 >tail2_tmp
    sed 's/H12R POPS/H12R OL  /' tail2_tmp >tail2
    sed 's/H12S POPS/H12S OL  /' tail2 >tail2_tmp
    sed 's/H13R POPS/H13R OL  /' tail2_tmp >tail2
    sed 's/H13S POPS/H13S OL  /' tail2 >tail2_tmp
    sed 's/H14R POPS/H14R OL  /' tail2_tmp >tail2
    sed 's/H14S POPS/H14S OL  /' tail2 >tail2_tmp
    sed 's/H15R POPS/H15R OL  /' tail2_tmp >tail2
    sed 's/H15S POPS/H15S OL  /' tail2 >tail2_tmp
    sed 's/H16R POPS/H16R OL  /' tail2_tmp >tail2
    sed 's/H16S POPS/H16S OL  /' tail2 >tail2_tmp
    sed 's/H17R POPS/H17R OL  /' tail2_tmp >tail2
    sed 's/H17S POPS/H17S OL  /' tail2 >tail2_tmp
    sed 's/H18R POPS/H18R OL  /' tail2_tmp >tail2
    sed 's/H18S POPS/H18S OL  /' tail2 >tail2_tmp
    sed 's/H18T POPS/H18T OL  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   POPS/ N31 PS  /' head > head_tmp
    sed 's/ C12 POPS/ C32 PS  /' head_tmp > head
    sed 's/ C13 POPS/ C33 PS  /' head > head_tmp
    sed 's/O13A POPS/ O35 PS  /' head_tmp > head
    sed 's/O13B POPS/ O36 PS  /' head > head_tmp
    sed 's/ O12 POPS/ O32 PS  /' head_tmp > head
    sed 's/ P   POPS/ P31 PS  /' head > head_tmp
    sed 's/ O11 POPS/ O31 PS  /' head_tmp > head
    sed 's/ O13 POPS/ O33 PS  /' head > head_tmp
    sed 's/ C11 POPS/ C31 PS  /' head_tmp > head
    sed 's/ O14 POPS/ O34 PS  /' head > head_tmp
    sed 's/ C1  POPS/ C3  PS  /' head_tmp > head
    sed 's/ C2  POPS/ C2  PS  /' head > head_tmp
    sed 's/ C3  POPS/ C1  PS  /' head_tmp > head
    sed 's/ O31 POPS/ O11 PS  /' head > head_tmp
    sed 's/ C31 POPS/ C11 PS  /' head_tmp > head
    sed 's/ O32 POPS/ O12 PS  /' head > head_tmp
    sed 's/ O21 POPS/ O21 PS  /' head_tmp > head
    sed 's/ C21 POPS/ C21 PS  /' head > head_tmp
    sed 's/ O22 POPS/ O22 PS  /' head_tmp > head
    sed 's/ HX  POPS/ HR  PS  /' head > head_tmp
    sed 's/ HY  POPS/ HS  PS  /' head_tmp > head
    sed 's/ HS  POPS/ HX  PS  /' head > head_tmp
    sed 's/ HA  POPS/ HA  PS  /' head_tmp > head
    sed 's/ HB  POPS/ HB  PS  /' head > head_tmp
    sed 's/H11A POPS/ H1A PS  /' head_tmp > head
    sed 's/H11B POPS/ H1B PS  /' head > head_tmp
    sed 's/H12A POPS/ H2A PS  /' head_tmp > head
    sed 's/ HN1 POPS/HN1A PS  /' head > head_tmp
    sed 's/ HN2 POPS/HN1B PS  /' head_tmp > head
    sed 's/ HN3 POPS/HN1C PS  /' head > head_tmp
    mv head_tmp head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#POPG
if ( $popg_found == 1 ) then
  #POPG is 127 atoms per lipid. 
  set atoms_per_lipid = 127

  echo " Processing: POPG"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l POPG_tmp | awk '{print$1}'`
  echo "  POPG Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted POPG Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " POPG Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../POPG_tmp
  rm -f ../POPG_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing POPG Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 POPG" $i > tail1
    grep -A42 "C33 POPG" $i >> tail1

    grep -A2 "C22 POPG" $i > tail2
    grep -A46 "C23 POPG" $i >> tail2

    head -n25 $i > head
    grep -A5 "C3  POPG" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 POPG/ C12 PA  /' tail1 >tail1_tmp
    sed 's/ C33 POPG/ C13 PA  /' tail1_tmp >tail1
    sed 's/ C34 POPG/ C14 PA  /' tail1 >tail1_tmp
    sed 's/ C35 POPG/ C15 PA  /' tail1_tmp >tail1
    sed 's/ C36 POPG/ C16 PA  /' tail1 >tail1_tmp
    sed 's/ C37 POPG/ C17 PA  /' tail1_tmp >tail1
    sed 's/ C38 POPG/ C18 PA  /' tail1 >tail1_tmp
    sed 's/ C39 POPG/ C19 PA  /' tail1_tmp >tail1
    sed 's/C310 POPG/C110 PA  /' tail1 >tail1_tmp
    sed 's/C311 POPG/C111 PA  /' tail1_tmp >tail1
    sed 's/C312 POPG/C112 PA  /' tail1 >tail1_tmp
    sed 's/C313 POPG/C113 PA  /' tail1_tmp >tail1
    sed 's/C314 POPG/C114 PA  /' tail1 >tail1_tmp
    sed 's/C315 POPG/C115 PA  /' tail1_tmp >tail1
    sed 's/C316 POPG/C116 PA  /' tail1 >tail1_tmp
    sed 's/C317 POPG/C117 PA  /' tail1_tmp >tail1
    sed 's/C318 POPG/C118 PA  /' tail1 >tail1_tmp
    sed 's/ H2X POPG/ H2R PA  /' tail1_tmp >tail1
    sed 's/ H2Y POPG/ H2S PA  /' tail1 >tail1_tmp
    sed 's/ H3X POPG/ H3R PA  /' tail1_tmp >tail1
    sed 's/ H3Y POPG/ H3S PA  /' tail1 >tail1_tmp
    sed 's/ H4X POPG/ H4R PA  /' tail1_tmp >tail1
    sed 's/ H4Y POPG/ H4S PA  /' tail1 >tail1_tmp
    sed 's/ H5X POPG/ H5R PA  /' tail1_tmp >tail1
    sed 's/ H5Y POPG/ H5S PA  /' tail1 >tail1_tmp
    sed 's/ H6X POPG/ H6R PA  /' tail1_tmp >tail1
    sed 's/ H6Y POPG/ H6S PA  /' tail1 >tail1_tmp
    sed 's/ H7X POPG/ H7R PA  /' tail1_tmp >tail1
    sed 's/ H7Y POPG/ H7S PA  /' tail1 >tail1_tmp
    sed 's/ H8X POPG/ H8R PA  /' tail1_tmp >tail1
    sed 's/ H8Y POPG/ H8S PA  /' tail1 >tail1_tmp
    sed 's/ H9X POPG/ H9R PA  /' tail1_tmp >tail1
    sed 's/ H9Y POPG/ H9S PA  /' tail1 >tail1_tmp
    sed 's/H10X POPG/H10R PA  /' tail1_tmp >tail1
    sed 's/H10Y POPG/H10S PA  /' tail1 >tail1_tmp
    sed 's/H11X POPG/H11R PA  /' tail1_tmp >tail1
    sed 's/H11Y POPG/H11S PA  /' tail1 >tail1_tmp
    sed 's/H12X POPG/H12R PA  /' tail1_tmp >tail1
    sed 's/H12Y POPG/H12S PA  /' tail1 >tail1_tmp
    sed 's/H13X POPG/H13R PA  /' tail1_tmp >tail1
    sed 's/H13Y POPG/H13S PA  /' tail1 >tail1_tmp
    sed 's/H14X POPG/H14R PA  /' tail1_tmp >tail1
    sed 's/H14Y POPG/H14S PA  /' tail1 >tail1_tmp
    sed 's/H15X POPG/H15R PA  /' tail1_tmp >tail1
    sed 's/H15Y POPG/H15S PA  /' tail1 >tail1_tmp
    sed 's/H16X POPG/H16R PA  /' tail1_tmp >tail1
    sed 's/H16Y POPG/H16S PA  /' tail1 >tail1_tmp
    sed 's/H16Z POPG/H16T PA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 POPG/ C12 OL  /' tail2 >tail2_tmp
    sed 's/ C23 POPG/ C13 OL  /' tail2_tmp >tail2
    sed 's/ C24 POPG/ C14 OL  /' tail2 >tail2_tmp
    sed 's/ C25 POPG/ C15 OL  /' tail2_tmp >tail2
    sed 's/ C26 POPG/ C16 OL  /' tail2 >tail2_tmp
    sed 's/ C27 POPG/ C17 OL  /' tail2_tmp >tail2
    sed 's/ C28 POPG/ C18 OL  /' tail2 >tail2_tmp
    sed 's/ C29 POPG/ C19 OL  /' tail2_tmp >tail2
    sed 's/C210 POPG/C110 OL  /' tail2 >tail2_tmp
    sed 's/C211 POPG/C111 OL  /' tail2_tmp >tail2
    sed 's/C212 POPG/C112 OL  /' tail2 >tail2_tmp
    sed 's/C213 POPG/C113 OL  /' tail2_tmp >tail2
    sed 's/C214 POPG/C114 OL  /' tail2 >tail2_tmp
    sed 's/C215 POPG/C115 OL  /' tail2_tmp >tail2
    sed 's/C216 POPG/C116 OL  /' tail2 >tail2_tmp
    sed 's/C217 POPG/C117 OL  /' tail2_tmp >tail2
    sed 's/C218 POPG/C118 OL  /' tail2 >tail2_tmp
    sed 's/ H2R POPG/ H2R OL  /' tail2_tmp >tail2
    sed 's/ H2S POPG/ H2S OL  /' tail2 >tail2_tmp
    sed 's/ H3R POPG/ H3R OL  /' tail2_tmp >tail2
    sed 's/ H3S POPG/ H3S OL  /' tail2 >tail2_tmp
    sed 's/ H4R POPG/ H4R OL  /' tail2_tmp >tail2
    sed 's/ H4S POPG/ H4S OL  /' tail2 >tail2_tmp
    sed 's/ H5R POPG/ H5R OL  /' tail2_tmp >tail2
    sed 's/ H5S POPG/ H5S OL  /' tail2 >tail2_tmp
    sed 's/ H6R POPG/ H6R OL  /' tail2_tmp >tail2
    sed 's/ H6S POPG/ H6S OL  /' tail2 >tail2_tmp
    sed 's/ H7R POPG/ H7R OL  /' tail2_tmp >tail2
    sed 's/ H7S POPG/ H7S OL  /' tail2 >tail2_tmp
    sed 's/ H8R POPG/ H8R OL  /' tail2_tmp >tail2
    sed 's/ H8S POPG/ H8S OL  /' tail2 >tail2_tmp
    sed 's/ H91 POPG/ H9R OL  /' tail2_tmp >tail2
    sed 's/H101 POPG/H10R OL  /' tail2 >tail2_tmp
    sed 's/H11R POPG/H11R OL  /' tail2_tmp >tail2
    sed 's/H11S POPG/H11S OL  /' tail2 >tail2_tmp
    sed 's/H12R POPG/H12R OL  /' tail2_tmp >tail2
    sed 's/H12S POPG/H12S OL  /' tail2 >tail2_tmp
    sed 's/H13R POPG/H13R OL  /' tail2_tmp >tail2
    sed 's/H13S POPG/H13S OL  /' tail2 >tail2_tmp
    sed 's/H14R POPG/H14R OL  /' tail2_tmp >tail2
    sed 's/H14S POPG/H14S OL  /' tail2 >tail2_tmp
    sed 's/H15R POPG/H15R OL  /' tail2_tmp >tail2
    sed 's/H15S POPG/H15S OL  /' tail2 >tail2_tmp
    sed 's/H16R POPG/H16R OL  /' tail2_tmp >tail2
    sed 's/H16S POPG/H16S OL  /' tail2 >tail2_tmp
    sed 's/H17R POPG/H17R OL  /' tail2_tmp >tail2
    sed 's/H17S POPG/H17S OL  /' tail2 >tail2_tmp
    sed 's/H18R POPG/H18R OL  /' tail2_tmp >tail2
    sed 's/H18S POPG/H18S OL  /' tail2 >tail2_tmp
    sed 's/H18T POPG/H18T OL  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ C13 POPG/ C33 PGR /' head > head_tmp
    sed 's/H13A POPG/ H3A PGR /' head_tmp > head
    sed 's/H13B POPG/ H3B PGR /' head > head_tmp
    sed 's/ OC3 POPG/ O36 PGR /' head_tmp > head
    sed 's/ HO3 POPG/HO6A PGR /' head > head_tmp
    sed 's/ C12 POPG/ C32 PGR /' head_tmp > head
    sed 's/H12A POPG/ H2A PGR /' head > head_tmp
    sed 's/ OC2 POPG/ O35 PGR /' head_tmp > head
    sed 's/ HO2 POPG/HO5A PGR /' head > head_tmp
    sed 's/ C11 POPG/ C31 PGR /' head_tmp > head
    sed 's/H11A POPG/ H1A PGR /' head > head_tmp
    sed 's/H11B POPG/ H1B PGR /' head_tmp > head
    sed 's/ P   POPG/ P31 PGR /' head > head_tmp
    sed 's/ O13 POPG/ O33 PGR /' head_tmp > head
    sed 's/ O14 POPG/ O34 PGR /' head > head_tmp
    sed 's/ O12 POPG/ O32 PGR /' head_tmp > head
    sed 's/ O11 POPG/ O31 PGR /' head > head_tmp
    sed 's/ C1  POPG/ C3  PGR /' head_tmp > head
    sed 's/ HA  POPG/ HA  PGR /' head > head_tmp
    sed 's/ HB  POPG/ HB  PGR /' head_tmp > head
    sed 's/ C2  POPG/ C2  PGR /' head > head_tmp
    sed 's/ HS  POPG/ HX  PGR /' head_tmp > head
    sed 's/ O21 POPG/ O21 PGR /' head > head_tmp
    sed 's/ C21 POPG/ C21 PGR /' head_tmp > head
    sed 's/ O22 POPG/ O22 PGR /' head > head_tmp
    sed 's/ C3  POPG/ C1  PGR /' head_tmp > head
    sed 's/ HX  POPG/ HR  PGR /' head > head_tmp
    sed 's/ HY  POPG/ HS  PGR /' head_tmp > head
    sed 's/ O31 POPG/ O11 PGR /' head > head_tmp
    sed 's/ C31 POPG/ C11 PGR /' head_tmp > head
    sed 's/ O32 POPG/ O12 PGR /' head > head_tmp
    mv head_tmp head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#POPA
if ( $popa_found == 1 ) then
  #POPA is 116 atoms per lipid. 
  set atoms_per_lipid = 116

  echo " Processing: POPA"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l POPA_tmp | awk '{print$1}'`
  echo "  POPA Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted POPA Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " POPA Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../POPA_tmp
  rm -f ../POPA_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing POPA Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 POPA" $i > tail1
    grep -A42 "C33 POPA" $i >> tail1

    grep -A2 "C22 POPA" $i > tail2
    grep -A46 "C23 POPA" $i >> tail2

    head -n14 $i > head
    grep -A5 "C3  POPA" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 POPA/ C12 PA  /' tail1 >tail1_tmp
    sed 's/ C33 POPA/ C13 PA  /' tail1_tmp >tail1
    sed 's/ C34 POPA/ C14 PA  /' tail1 >tail1_tmp
    sed 's/ C35 POPA/ C15 PA  /' tail1_tmp >tail1
    sed 's/ C36 POPA/ C16 PA  /' tail1 >tail1_tmp
    sed 's/ C37 POPA/ C17 PA  /' tail1_tmp >tail1
    sed 's/ C38 POPA/ C18 PA  /' tail1 >tail1_tmp
    sed 's/ C39 POPA/ C19 PA  /' tail1_tmp >tail1
    sed 's/C310 POPA/C110 PA  /' tail1 >tail1_tmp
    sed 's/C311 POPA/C111 PA  /' tail1_tmp >tail1
    sed 's/C312 POPA/C112 PA  /' tail1 >tail1_tmp
    sed 's/C313 POPA/C113 PA  /' tail1_tmp >tail1
    sed 's/C314 POPA/C114 PA  /' tail1 >tail1_tmp
    sed 's/C315 POPA/C115 PA  /' tail1_tmp >tail1
    sed 's/C316 POPA/C116 PA  /' tail1 >tail1_tmp
    sed 's/C317 POPA/C117 PA  /' tail1_tmp >tail1
    sed 's/C318 POPA/C118 PA  /' tail1 >tail1_tmp
    sed 's/ H2X POPA/ H2R PA  /' tail1_tmp >tail1
    sed 's/ H2Y POPA/ H2S PA  /' tail1 >tail1_tmp
    sed 's/ H3X POPA/ H3R PA  /' tail1_tmp >tail1
    sed 's/ H3Y POPA/ H3S PA  /' tail1 >tail1_tmp
    sed 's/ H4X POPA/ H4R PA  /' tail1_tmp >tail1
    sed 's/ H4Y POPA/ H4S PA  /' tail1 >tail1_tmp
    sed 's/ H5X POPA/ H5R PA  /' tail1_tmp >tail1
    sed 's/ H5Y POPA/ H5S PA  /' tail1 >tail1_tmp
    sed 's/ H6X POPA/ H6R PA  /' tail1_tmp >tail1
    sed 's/ H6Y POPA/ H6S PA  /' tail1 >tail1_tmp
    sed 's/ H7X POPA/ H7R PA  /' tail1_tmp >tail1
    sed 's/ H7Y POPA/ H7S PA  /' tail1 >tail1_tmp
    sed 's/ H8X POPA/ H8R PA  /' tail1_tmp >tail1
    sed 's/ H8Y POPA/ H8S PA  /' tail1 >tail1_tmp
    sed 's/ H9X POPA/ H9R PA  /' tail1_tmp >tail1
    sed 's/ H9Y POPA/ H9S PA  /' tail1 >tail1_tmp
    sed 's/H10X POPA/H10R PA  /' tail1_tmp >tail1
    sed 's/H10Y POPA/H10S PA  /' tail1 >tail1_tmp
    sed 's/H11X POPA/H11R PA  /' tail1_tmp >tail1
    sed 's/H11Y POPA/H11S PA  /' tail1 >tail1_tmp
    sed 's/H12X POPA/H12R PA  /' tail1_tmp >tail1
    sed 's/H12Y POPA/H12S PA  /' tail1 >tail1_tmp
    sed 's/H13X POPA/H13R PA  /' tail1_tmp >tail1
    sed 's/H13Y POPA/H13S PA  /' tail1 >tail1_tmp
    sed 's/H14X POPA/H14R PA  /' tail1_tmp >tail1
    sed 's/H14Y POPA/H14S PA  /' tail1 >tail1_tmp
    sed 's/H15X POPA/H15R PA  /' tail1_tmp >tail1
    sed 's/H15Y POPA/H15S PA  /' tail1 >tail1_tmp
    sed 's/H16X POPA/H16R PA  /' tail1_tmp >tail1
    sed 's/H16Y POPA/H16S PA  /' tail1 >tail1_tmp
    sed 's/H16Z POPA/H16T PA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 POPA/ C12 OL  /' tail2 >tail2_tmp
    sed 's/ C23 POPA/ C13 OL  /' tail2_tmp >tail2
    sed 's/ C24 POPA/ C14 OL  /' tail2 >tail2_tmp
    sed 's/ C25 POPA/ C15 OL  /' tail2_tmp >tail2
    sed 's/ C26 POPA/ C16 OL  /' tail2 >tail2_tmp
    sed 's/ C27 POPA/ C17 OL  /' tail2_tmp >tail2
    sed 's/ C28 POPA/ C18 OL  /' tail2 >tail2_tmp
    sed 's/ C29 POPA/ C19 OL  /' tail2_tmp >tail2
    sed 's/C210 POPA/C110 OL  /' tail2 >tail2_tmp
    sed 's/C211 POPA/C111 OL  /' tail2_tmp >tail2
    sed 's/C212 POPA/C112 OL  /' tail2 >tail2_tmp
    sed 's/C213 POPA/C113 OL  /' tail2_tmp >tail2
    sed 's/C214 POPA/C114 OL  /' tail2 >tail2_tmp
    sed 's/C215 POPA/C115 OL  /' tail2_tmp >tail2
    sed 's/C216 POPA/C116 OL  /' tail2 >tail2_tmp
    sed 's/C217 POPA/C117 OL  /' tail2_tmp >tail2
    sed 's/C218 POPA/C118 OL  /' tail2 >tail2_tmp
    sed 's/ H2R POPA/ H2R OL  /' tail2_tmp >tail2
    sed 's/ H2S POPA/ H2S OL  /' tail2 >tail2_tmp
    sed 's/ H3R POPA/ H3R OL  /' tail2_tmp >tail2
    sed 's/ H3S POPA/ H3S OL  /' tail2 >tail2_tmp
    sed 's/ H4R POPA/ H4R OL  /' tail2_tmp >tail2
    sed 's/ H4S POPA/ H4S OL  /' tail2 >tail2_tmp
    sed 's/ H5R POPA/ H5R OL  /' tail2_tmp >tail2
    sed 's/ H5S POPA/ H5S OL  /' tail2 >tail2_tmp
    sed 's/ H6R POPA/ H6R OL  /' tail2_tmp >tail2
    sed 's/ H6S POPA/ H6S OL  /' tail2 >tail2_tmp
    sed 's/ H7R POPA/ H7R OL  /' tail2_tmp >tail2
    sed 's/ H7S POPA/ H7S OL  /' tail2 >tail2_tmp
    sed 's/ H8R POPA/ H8R OL  /' tail2_tmp >tail2
    sed 's/ H8S POPA/ H8S OL  /' tail2 >tail2_tmp
    sed 's/ H91 POPA/ H9R OL  /' tail2_tmp >tail2
    sed 's/H101 POPA/H10R OL  /' tail2 >tail2_tmp
    sed 's/H11R POPA/H11R OL  /' tail2_tmp >tail2
    sed 's/H11S POPA/H11S OL  /' tail2 >tail2_tmp
    sed 's/H12R POPA/H12R OL  /' tail2_tmp >tail2
    sed 's/H12S POPA/H12S OL  /' tail2 >tail2_tmp
    sed 's/H13R POPA/H13R OL  /' tail2_tmp >tail2
    sed 's/H13S POPA/H13S OL  /' tail2 >tail2_tmp
    sed 's/H14R POPA/H14R OL  /' tail2_tmp >tail2
    sed 's/H14S POPA/H14S OL  /' tail2 >tail2_tmp
    sed 's/H15R POPA/H15R OL  /' tail2_tmp >tail2
    sed 's/H15S POPA/H15S OL  /' tail2 >tail2_tmp
    sed 's/H16R POPA/H16R OL  /' tail2_tmp >tail2
    sed 's/H16S POPA/H16S OL  /' tail2 >tail2_tmp
    sed 's/H17R POPA/H17R OL  /' tail2_tmp >tail2
    sed 's/H17S POPA/H17S OL  /' tail2 >tail2_tmp
    sed 's/H18R POPA/H18R OL  /' tail2_tmp >tail2
    sed 's/H18S POPA/H18S OL  /' tail2 >tail2_tmp
    sed 's/H18T POPA/H18T OL  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ P   POPA/ P31 PH- /' head > head_tmp
    sed 's/ O13 POPA/ O33 PH- /' head_tmp > head
    sed 's/ O14 POPA/ O34 PH- /' head > head_tmp
    sed 's/ O12 POPA/ O32 PH- /' head_tmp > head
    sed 's/ H12 POPA/HO2A PH- /' head > head_tmp
    sed 's/ O11 POPA/ O31 PH- /' head_tmp > head
    sed 's/ C1  POPA/ C3  PH- /' head > head_tmp
    sed 's/ HA  POPA/ HA  PH- /' head_tmp > head
    sed 's/ HB  POPA/ HB  PH- /' head > head_tmp
    sed 's/ C2  POPA/ C2  PH- /' head_tmp > head
    sed 's/ HS  POPA/ HX  PH- /' head > head_tmp
    sed 's/ O21 POPA/ O21 PH- /' head_tmp > head
    sed 's/ C21 POPA/ C21 PH- /' head > head_tmp
    sed 's/ O22 POPA/ O22 PH- /' head_tmp > head
    sed 's/ C3  POPA/ C1  PH- /' head > head_tmp
    sed 's/ HX  POPA/ HR  PH- /' head_tmp > head
    sed 's/ HY  POPA/ HS  PH- /' head > head_tmp
    sed 's/ O31 POPA/ O11 PH- /' head_tmp > head
    sed 's/ C31 POPA/ C11 PH- /' head > head_tmp
    sed 's/ O32 POPA/ O12 PH- /' head_tmp > head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DLPC
if ( $dlpc_found == 1 ) then
  #DLPC is 106 atoms per lipid. 
  set atoms_per_lipid = 106

  echo " Processing: DLPC"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DLPC_tmp | awk '{print$1}'`
  echo "  DLPC Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DLPC Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DLPC Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DLPC_tmp
  rm -f ../DLPC_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DLPC Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DLPC" $i > tail1
    grep -A30 "C33 DLPC" $i >> tail1

    grep -A2 "C22 DLPC" $i > tail2
    grep -A30 "C23 DLPC" $i >> tail2

    head -n32 $i > head
    grep -A5 "C3  DLPC" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DLPC/ C12 LA  /' tail1 >tail1_tmp
    sed 's/ C33 DLPC/ C13 LA  /' tail1_tmp >tail1
    sed 's/ C34 DLPC/ C14 LA  /' tail1 >tail1_tmp
    sed 's/ C35 DLPC/ C15 LA  /' tail1_tmp >tail1
    sed 's/ C36 DLPC/ C16 LA  /' tail1 >tail1_tmp
    sed 's/ C37 DLPC/ C17 LA  /' tail1_tmp >tail1
    sed 's/ C38 DLPC/ C18 LA  /' tail1 >tail1_tmp
    sed 's/ C39 DLPC/ C19 LA  /' tail1_tmp >tail1
    sed 's/C310 DLPC/C110 LA  /' tail1 >tail1_tmp
    sed 's/C311 DLPC/C111 LA  /' tail1_tmp >tail1
    sed 's/C312 DLPC/C112 LA  /' tail1 >tail1_tmp
    sed 's/ H2X DLPC/ H2R LA  /' tail1_tmp >tail1
    sed 's/ H2Y DLPC/ H2S LA  /' tail1 >tail1_tmp
    sed 's/ H3X DLPC/ H3R LA  /' tail1_tmp >tail1
    sed 's/ H3Y DLPC/ H3S LA  /' tail1 >tail1_tmp
    sed 's/ H4X DLPC/ H4R LA  /' tail1_tmp >tail1
    sed 's/ H4Y DLPC/ H4S LA  /' tail1 >tail1_tmp
    sed 's/ H5X DLPC/ H5R LA  /' tail1_tmp >tail1
    sed 's/ H5Y DLPC/ H5S LA  /' tail1 >tail1_tmp
    sed 's/ H6X DLPC/ H6R LA  /' tail1_tmp >tail1
    sed 's/ H6Y DLPC/ H6S LA  /' tail1 >tail1_tmp
    sed 's/ H7X DLPC/ H7R LA  /' tail1_tmp >tail1
    sed 's/ H7Y DLPC/ H7S LA  /' tail1 >tail1_tmp
    sed 's/ H8X DLPC/ H8R LA  /' tail1_tmp >tail1
    sed 's/ H8Y DLPC/ H8S LA  /' tail1 >tail1_tmp
    sed 's/ H9X DLPC/ H9R LA  /' tail1_tmp >tail1
    sed 's/ H9Y DLPC/ H9S LA  /' tail1 >tail1_tmp
    sed 's/H10X DLPC/H10R LA  /' tail1_tmp >tail1
    sed 's/H10Y DLPC/H10S LA  /' tail1 >tail1_tmp
    sed 's/H11X DLPC/H11R LA  /' tail1_tmp >tail1
    sed 's/H11Y DLPC/H11S LA  /' tail1 >tail1_tmp
    sed 's/H12X DLPC/H12R LA  /' tail1_tmp >tail1
    sed 's/H12Y DLPC/H12S LA  /' tail1 >tail1_tmp
    sed 's/H12Z DLPC/H12T LA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DLPC/ C12 LA  /' tail2 >tail2_tmp
    sed 's/ C23 DLPC/ C13 LA  /' tail2_tmp >tail2
    sed 's/ C24 DLPC/ C14 LA  /' tail2 >tail2_tmp
    sed 's/ C25 DLPC/ C15 LA  /' tail2_tmp >tail2
    sed 's/ C26 DLPC/ C16 LA  /' tail2 >tail2_tmp
    sed 's/ C27 DLPC/ C17 LA  /' tail2_tmp >tail2
    sed 's/ C28 DLPC/ C18 LA  /' tail2 >tail2_tmp
    sed 's/ C29 DLPC/ C19 LA  /' tail2_tmp >tail2
    sed 's/C210 DLPC/C110 LA  /' tail2 >tail2_tmp
    sed 's/C211 DLPC/C111 LA  /' tail2_tmp >tail2
    sed 's/C212 DLPC/C112 LA  /' tail2 >tail2_tmp
    sed 's/ H2R DLPC/ H2R LA  /' tail2_tmp >tail2
    sed 's/ H2S DLPC/ H2S LA  /' tail2 >tail2_tmp
    sed 's/ H3R DLPC/ H3R LA  /' tail2_tmp >tail2
    sed 's/ H3S DLPC/ H3S LA  /' tail2 >tail2_tmp
    sed 's/ H4R DLPC/ H4R LA  /' tail2_tmp >tail2
    sed 's/ H4S DLPC/ H4S LA  /' tail2 >tail2_tmp
    sed 's/ H5R DLPC/ H5R LA  /' tail2_tmp >tail2
    sed 's/ H5S DLPC/ H5S LA  /' tail2 >tail2_tmp
    sed 's/ H6R DLPC/ H6R LA  /' tail2_tmp >tail2
    sed 's/ H6S DLPC/ H6S LA  /' tail2 >tail2_tmp
    sed 's/ H7R DLPC/ H7R LA  /' tail2_tmp >tail2
    sed 's/ H7S DLPC/ H7S LA  /' tail2 >tail2_tmp
    sed 's/ H8R DLPC/ H8R LA  /' tail2_tmp >tail2
    sed 's/ H8S DLPC/ H8S LA  /' tail2 >tail2_tmp
    sed 's/ H9R DLPC/ H9R LA  /' tail2_tmp >tail2
    sed 's/ H9S DLPC/ H9S LA  /' tail2 >tail2_tmp
    sed 's/H10R DLPC/H10R LA  /' tail2_tmp >tail2
    sed 's/H10S DLPC/H10S LA  /' tail2 >tail2_tmp
    sed 's/H11R DLPC/H11R LA  /' tail2_tmp >tail2
    sed 's/H11S DLPC/H11S LA  /' tail2 >tail2_tmp
    sed 's/H12R DLPC/H12R LA  /' tail2_tmp >tail2
    sed 's/H12S DLPC/H12S LA  /' tail2 >tail2_tmp
    sed 's/H12T DLPC/H12T LA  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DLPC/ N31 PC  /' head > head_tmp
    sed 's/ C13 DLPC/ C33 PC  /' head_tmp > head
    sed 's/ C14 DLPC/ C34 PC  /' head > head_tmp
    sed 's/ C15 DLPC/ C35 PC  /' head_tmp > head
    sed 's/ C12 DLPC/ C32 PC  /' head > head_tmp
    sed 's/ C11 DLPC/ C31 PC  /' head_tmp > head
    sed 's/ O12 DLPC/ O32 PC  /' head > head_tmp
    sed 's/ P   DLPC/ P31 PC  /' head_tmp > head
    sed 's/ O11 DLPC/ O31 PC  /' head > head_tmp
    sed 's/ O13 DLPC/ O33 PC  /' head_tmp > head
    sed 's/ O14 DLPC/ O34 PC  /' head > head_tmp
    sed 's/ C1  DLPC/ C3  PC  /' head_tmp > head
    sed 's/ C2  DLPC/ C2  PC  /' head > head_tmp
    sed 's/ C3  DLPC/ C1  PC  /' head_tmp > head
    sed 's/ O31 DLPC/ O11 PC  /' head > head_tmp
    sed 's/ C31 DLPC/ C11 PC  /' head_tmp > head
    sed 's/ O32 DLPC/ O12 PC  /' head > head_tmp
    sed 's/ O21 DLPC/ O21 PC  /' head_tmp > head
    sed 's/ C21 DLPC/ C21 PC  /' head > head_tmp
    sed 's/ O22 DLPC/ O22 PC  /' head_tmp > head
    sed 's/ HX  DLPC/ HR  PC  /' head > head_tmp
    sed 's/ HY  DLPC/ HS  PC  /' head_tmp > head
    sed 's/ HS  DLPC/ HX  PC  /' head > head_tmp
    sed 's/ HA  DLPC/ HA  PC  /' head_tmp > head
    sed 's/ HB  DLPC/ HB  PC  /' head > head_tmp
    sed 's/H11A DLPC/ H1A PC  /' head_tmp > head
    sed 's/H11B DLPC/ H1B PC  /' head > head_tmp
    sed 's/H12A DLPC/ H2A PC  /' head_tmp > head
    sed 's/H12B DLPC/ H2B PC  /' head > head_tmp
    sed 's/H13A DLPC/ H3A PC  /' head_tmp > head
    sed 's/H13B DLPC/ H3B PC  /' head > head_tmp
    sed 's/H13C DLPC/ H3C PC  /' head_tmp > head
    sed 's/H14A DLPC/ H4A PC  /' head > head_tmp
    sed 's/H14B DLPC/ H4B PC  /' head_tmp > head
    sed 's/H14C DLPC/ H4C PC  /' head > head_tmp
    sed 's/H15A DLPC/ H5A PC  /' head_tmp > head
    sed 's/H15B DLPC/ H5B PC  /' head > head_tmp
    sed 's/H15C DLPC/ H5C PC  /' head_tmp > head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DLPE
if ( $dlpe_found == 1 ) then
  #DLPE is 97 atoms per lipid. 
  set atoms_per_lipid = 97

  echo " Processing: DLPE"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DLPE_tmp | awk '{print$1}'`
  echo "  DLPE Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DLPE Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DLPE Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DLPE_tmp
  rm -f ../DLPE_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DLPE Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DLPE" $i > tail1
    grep -A30 "C33 DLPE" $i >> tail1

    grep -A2 "C22 DLPE" $i > tail2
    grep -A30 "C23 DLPE" $i >> tail2

    head -n23 $i > head
    grep -A5 "C3  DLPE" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DLPE/ C12 LA  /' tail1 >tail1_tmp
    sed 's/ C33 DLPE/ C13 LA  /' tail1_tmp >tail1
    sed 's/ C34 DLPE/ C14 LA  /' tail1 >tail1_tmp
    sed 's/ C35 DLPE/ C15 LA  /' tail1_tmp >tail1
    sed 's/ C36 DLPE/ C16 LA  /' tail1 >tail1_tmp
    sed 's/ C37 DLPE/ C17 LA  /' tail1_tmp >tail1
    sed 's/ C38 DLPE/ C18 LA  /' tail1 >tail1_tmp
    sed 's/ C39 DLPE/ C19 LA  /' tail1_tmp >tail1
    sed 's/C310 DLPE/C110 LA  /' tail1 >tail1_tmp
    sed 's/C311 DLPE/C111 LA  /' tail1_tmp >tail1
    sed 's/C312 DLPE/C112 LA  /' tail1 >tail1_tmp
    sed 's/ H2X DLPE/ H2R LA  /' tail1_tmp >tail1
    sed 's/ H2Y DLPE/ H2S LA  /' tail1 >tail1_tmp
    sed 's/ H3X DLPE/ H3R LA  /' tail1_tmp >tail1
    sed 's/ H3Y DLPE/ H3S LA  /' tail1 >tail1_tmp
    sed 's/ H4X DLPE/ H4R LA  /' tail1_tmp >tail1
    sed 's/ H4Y DLPE/ H4S LA  /' tail1 >tail1_tmp
    sed 's/ H5X DLPE/ H5R LA  /' tail1_tmp >tail1
    sed 's/ H5Y DLPE/ H5S LA  /' tail1 >tail1_tmp
    sed 's/ H6X DLPE/ H6R LA  /' tail1_tmp >tail1
    sed 's/ H6Y DLPE/ H6S LA  /' tail1 >tail1_tmp
    sed 's/ H7X DLPE/ H7R LA  /' tail1_tmp >tail1
    sed 's/ H7Y DLPE/ H7S LA  /' tail1 >tail1_tmp
    sed 's/ H8X DLPE/ H8R LA  /' tail1_tmp >tail1
    sed 's/ H8Y DLPE/ H8S LA  /' tail1 >tail1_tmp
    sed 's/ H9X DLPE/ H9R LA  /' tail1_tmp >tail1
    sed 's/ H9Y DLPE/ H9S LA  /' tail1 >tail1_tmp
    sed 's/H10X DLPE/H10R LA  /' tail1_tmp >tail1
    sed 's/H10Y DLPE/H10S LA  /' tail1 >tail1_tmp
    sed 's/H11X DLPE/H11R LA  /' tail1_tmp >tail1
    sed 's/H11Y DLPE/H11S LA  /' tail1 >tail1_tmp
    sed 's/H12X DLPE/H12R LA  /' tail1_tmp >tail1
    sed 's/H12Y DLPE/H12S LA  /' tail1 >tail1_tmp
    sed 's/H12Z DLPE/H12T LA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DLPE/ C12 LA  /' tail2 >tail2_tmp
    sed 's/ C23 DLPE/ C13 LA  /' tail2_tmp >tail2
    sed 's/ C24 DLPE/ C14 LA  /' tail2 >tail2_tmp
    sed 's/ C25 DLPE/ C15 LA  /' tail2_tmp >tail2
    sed 's/ C26 DLPE/ C16 LA  /' tail2 >tail2_tmp
    sed 's/ C27 DLPE/ C17 LA  /' tail2_tmp >tail2
    sed 's/ C28 DLPE/ C18 LA  /' tail2 >tail2_tmp
    sed 's/ C29 DLPE/ C19 LA  /' tail2_tmp >tail2
    sed 's/C210 DLPE/C110 LA  /' tail2 >tail2_tmp
    sed 's/C211 DLPE/C111 LA  /' tail2_tmp >tail2
    sed 's/C212 DLPE/C112 LA  /' tail2 >tail2_tmp
    sed 's/ H2R DLPE/ H2R LA  /' tail2_tmp >tail2
    sed 's/ H2S DLPE/ H2S LA  /' tail2 >tail2_tmp
    sed 's/ H3R DLPE/ H3R LA  /' tail2_tmp >tail2
    sed 's/ H3S DLPE/ H3S LA  /' tail2 >tail2_tmp
    sed 's/ H4R DLPE/ H4R LA  /' tail2_tmp >tail2
    sed 's/ H4S DLPE/ H4S LA  /' tail2 >tail2_tmp
    sed 's/ H5R DLPE/ H5R LA  /' tail2_tmp >tail2
    sed 's/ H5S DLPE/ H5S LA  /' tail2 >tail2_tmp
    sed 's/ H6R DLPE/ H6R LA  /' tail2_tmp >tail2
    sed 's/ H6S DLPE/ H6S LA  /' tail2 >tail2_tmp
    sed 's/ H7R DLPE/ H7R LA  /' tail2_tmp >tail2
    sed 's/ H7S DLPE/ H7S LA  /' tail2 >tail2_tmp
    sed 's/ H8R DLPE/ H8R LA  /' tail2_tmp >tail2
    sed 's/ H8S DLPE/ H8S LA  /' tail2 >tail2_tmp
    sed 's/ H9R DLPE/ H9R LA  /' tail2_tmp >tail2
    sed 's/ H9S DLPE/ H9S LA  /' tail2 >tail2_tmp
    sed 's/H10R DLPE/H10R LA  /' tail2_tmp >tail2
    sed 's/H10S DLPE/H10S LA  /' tail2 >tail2_tmp
    sed 's/H11R DLPE/H11R LA  /' tail2_tmp >tail2
    sed 's/H11S DLPE/H11S LA  /' tail2 >tail2_tmp
    sed 's/H12R DLPE/H12R LA  /' tail2_tmp >tail2
    sed 's/H12S DLPE/H12S LA  /' tail2 >tail2_tmp
    sed 's/H12T DLPE/H12T LA  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DLPE/ N31 PE  /' head > head_tmp
    sed 's/ C12 DLPE/ C32 PE  /' head_tmp > head
    sed 's/ C11 DLPE/ C31 PE  /' head > head_tmp
    sed 's/ O12 DLPE/ O32 PE  /' head_tmp > head
    sed 's/ P   DLPE/ P31 PE  /' head > head_tmp
    sed 's/ O11 DLPE/ O31 PE  /' head_tmp > head
    sed 's/ O13 DLPE/ O33 PE  /' head > head_tmp
    sed 's/ O14 DLPE/ O34 PE  /' head_tmp > head
    sed 's/ C1  DLPE/ C3  PE  /' head > head_tmp
    sed 's/ C2  DLPE/ C2  PE  /' head_tmp > head
    sed 's/ C3  DLPE/ C1  PE  /' head > head_tmp
    sed 's/ O31 DLPE/ O11 PE  /' head_tmp > head
    sed 's/ C31 DLPE/ C11 PE  /' head > head_tmp
    sed 's/ O32 DLPE/ O12 PE  /' head_tmp > head
    sed 's/ O21 DLPE/ O21 PE  /' head > head_tmp
    sed 's/ C21 DLPE/ C21 PE  /' head_tmp > head
    sed 's/ O22 DLPE/ O22 PE  /' head > head_tmp
    sed 's/ HX  DLPE/ HR  PE  /' head_tmp > head
    sed 's/ HY  DLPE/ HS  PE  /' head > head_tmp
    sed 's/ HS  DLPE/ HX  PE  /' head_tmp > head
    sed 's/ HA  DLPE/ HA  PE  /' head > head_tmp
    sed 's/ HB  DLPE/ HB  PE  /' head_tmp > head
    sed 's/H11A DLPE/ H1A PE  /' head > head_tmp
    sed 's/H11B DLPE/ H1B PE  /' head_tmp > head
    sed 's/H12A DLPE/ H2A PE  /' head > head_tmp
    sed 's/H12B DLPE/ H2B PE  /' head_tmp > head
    sed 's/ HN1 DLPE/HN1A PE  /' head > head_tmp
    sed 's/ HN2 DLPE/HN1B PE  /' head_tmp > head
    sed 's/ HN3 DLPE/HN1C PE  /' head > head_tmp
    mv head_tmp head 
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DLPS
if ( $dlps_found == 1 ) then
  #DLPS is 99 atoms per lipid. 
  set atoms_per_lipid = 99

  echo " Processing: DLPS"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DLPS_tmp | awk '{print$1}'`
  echo "  DLPS Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DLPS Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DLPS Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DLPS_tmp
  rm -f ../DLPS_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DLPS Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DLPS" $i > tail1
    grep -A30 "C33 DLPS" $i >> tail1

    grep -A2 "C22 DLPS" $i > tail2
    grep -A30 "C23 DLPS" $i >> tail2

    head -n25 $i > head
    grep -A5 "C3  DLPS" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DLPS/ C12 LA  /' tail1 >tail1_tmp
    sed 's/ C33 DLPS/ C13 LA  /' tail1_tmp >tail1
    sed 's/ C34 DLPS/ C14 LA  /' tail1 >tail1_tmp
    sed 's/ C35 DLPS/ C15 LA  /' tail1_tmp >tail1
    sed 's/ C36 DLPS/ C16 LA  /' tail1 >tail1_tmp
    sed 's/ C37 DLPS/ C17 LA  /' tail1_tmp >tail1
    sed 's/ C38 DLPS/ C18 LA  /' tail1 >tail1_tmp
    sed 's/ C39 DLPS/ C19 LA  /' tail1_tmp >tail1
    sed 's/C310 DLPS/C110 LA  /' tail1 >tail1_tmp
    sed 's/C311 DLPS/C111 LA  /' tail1_tmp >tail1
    sed 's/C312 DLPS/C112 LA  /' tail1 >tail1_tmp
    sed 's/ H2X DLPS/ H2R LA  /' tail1_tmp >tail1
    sed 's/ H2Y DLPS/ H2S LA  /' tail1 >tail1_tmp
    sed 's/ H3X DLPS/ H3R LA  /' tail1_tmp >tail1
    sed 's/ H3Y DLPS/ H3S LA  /' tail1 >tail1_tmp
    sed 's/ H4X DLPS/ H4R LA  /' tail1_tmp >tail1
    sed 's/ H4Y DLPS/ H4S LA  /' tail1 >tail1_tmp
    sed 's/ H5X DLPS/ H5R LA  /' tail1_tmp >tail1
    sed 's/ H5Y DLPS/ H5S LA  /' tail1 >tail1_tmp
    sed 's/ H6X DLPS/ H6R LA  /' tail1_tmp >tail1
    sed 's/ H6Y DLPS/ H6S LA  /' tail1 >tail1_tmp
    sed 's/ H7X DLPS/ H7R LA  /' tail1_tmp >tail1
    sed 's/ H7Y DLPS/ H7S LA  /' tail1 >tail1_tmp
    sed 's/ H8X DLPS/ H8R LA  /' tail1_tmp >tail1
    sed 's/ H8Y DLPS/ H8S LA  /' tail1 >tail1_tmp
    sed 's/ H9X DLPS/ H9R LA  /' tail1_tmp >tail1
    sed 's/ H9Y DLPS/ H9S LA  /' tail1 >tail1_tmp
    sed 's/H10X DLPS/H10R LA  /' tail1_tmp >tail1
    sed 's/H10Y DLPS/H10S LA  /' tail1 >tail1_tmp
    sed 's/H11X DLPS/H11R LA  /' tail1_tmp >tail1
    sed 's/H11Y DLPS/H11S LA  /' tail1 >tail1_tmp
    sed 's/H12X DLPS/H12R LA  /' tail1_tmp >tail1
    sed 's/H12Y DLPS/H12S LA  /' tail1 >tail1_tmp
    sed 's/H12Z DLPS/H12T LA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DLPS/ C12 LA  /' tail2 >tail2_tmp
    sed 's/ C23 DLPS/ C13 LA  /' tail2_tmp >tail2
    sed 's/ C24 DLPS/ C14 LA  /' tail2 >tail2_tmp
    sed 's/ C25 DLPS/ C15 LA  /' tail2_tmp >tail2
    sed 's/ C26 DLPS/ C16 LA  /' tail2 >tail2_tmp
    sed 's/ C27 DLPS/ C17 LA  /' tail2_tmp >tail2
    sed 's/ C28 DLPS/ C18 LA  /' tail2 >tail2_tmp
    sed 's/ C29 DLPS/ C19 LA  /' tail2_tmp >tail2
    sed 's/C210 DLPS/C110 LA  /' tail2 >tail2_tmp
    sed 's/C211 DLPS/C111 LA  /' tail2_tmp >tail2
    sed 's/C212 DLPS/C112 LA  /' tail2 >tail2_tmp
    sed 's/ H2R DLPS/ H2R LA  /' tail2_tmp >tail2
    sed 's/ H2S DLPS/ H2S LA  /' tail2 >tail2_tmp
    sed 's/ H3R DLPS/ H3R LA  /' tail2_tmp >tail2
    sed 's/ H3S DLPS/ H3S LA  /' tail2 >tail2_tmp
    sed 's/ H4R DLPS/ H4R LA  /' tail2_tmp >tail2
    sed 's/ H4S DLPS/ H4S LA  /' tail2 >tail2_tmp
    sed 's/ H5R DLPS/ H5R LA  /' tail2_tmp >tail2
    sed 's/ H5S DLPS/ H5S LA  /' tail2 >tail2_tmp
    sed 's/ H6R DLPS/ H6R LA  /' tail2_tmp >tail2
    sed 's/ H6S DLPS/ H6S LA  /' tail2 >tail2_tmp
    sed 's/ H7R DLPS/ H7R LA  /' tail2_tmp >tail2
    sed 's/ H7S DLPS/ H7S LA  /' tail2 >tail2_tmp
    sed 's/ H8R DLPS/ H8R LA  /' tail2_tmp >tail2
    sed 's/ H8S DLPS/ H8S LA  /' tail2 >tail2_tmp
    sed 's/ H9R DLPS/ H9R LA  /' tail2_tmp >tail2
    sed 's/ H9S DLPS/ H9S LA  /' tail2 >tail2_tmp
    sed 's/H10R DLPS/H10R LA  /' tail2_tmp >tail2
    sed 's/H10S DLPS/H10S LA  /' tail2 >tail2_tmp
    sed 's/H11R DLPS/H11R LA  /' tail2_tmp >tail2
    sed 's/H11S DLPS/H11S LA  /' tail2 >tail2_tmp
    sed 's/H12R DLPS/H12R LA  /' tail2_tmp >tail2
    sed 's/H12S DLPS/H12S LA  /' tail2 >tail2_tmp
    sed 's/H12T DLPS/H12T LA  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DLPS/ N31 PS  /' head > head_tmp
    sed 's/ C12 DLPS/ C32 PS  /' head_tmp > head
    sed 's/ C13 DLPS/ C33 PS  /' head > head_tmp
    sed 's/O13A DLPS/ O35 PS  /' head_tmp > head
    sed 's/O13B DLPS/ O36 PS  /' head > head_tmp
    sed 's/ O12 DLPS/ O32 PS  /' head_tmp > head
    sed 's/ P   DLPS/ P31 PS  /' head > head_tmp
    sed 's/ O11 DLPS/ O31 PS  /' head_tmp > head
    sed 's/ O13 DLPS/ O33 PS  /' head > head_tmp
    sed 's/ C11 DLPS/ C31 PS  /' head_tmp > head
    sed 's/ O14 DLPS/ O34 PS  /' head > head_tmp
    sed 's/ C1  DLPS/ C3  PS  /' head_tmp > head
    sed 's/ C2  DLPS/ C2  PS  /' head > head_tmp
    sed 's/ C3  DLPS/ C1  PS  /' head_tmp > head
    sed 's/ O31 DLPS/ O11 PS  /' head > head_tmp
    sed 's/ C31 DLPS/ C11 PS  /' head_tmp > head
    sed 's/ O32 DLPS/ O12 PS  /' head > head_tmp
    sed 's/ O21 DLPS/ O21 PS  /' head_tmp > head
    sed 's/ C21 DLPS/ C21 PS  /' head > head_tmp
    sed 's/ O22 DLPS/ O22 PS  /' head_tmp > head
    sed 's/ HX  DLPS/ HR  PS  /' head > head_tmp
    sed 's/ HY  DLPS/ HS  PS  /' head_tmp > head
    sed 's/ HS  DLPS/ HX  PS  /' head > head_tmp
    sed 's/ HA  DLPS/ HA  PS  /' head_tmp > head
    sed 's/ HB  DLPS/ HB  PS  /' head > head_tmp
    sed 's/H11A DLPS/ H1A PS  /' head_tmp > head
    sed 's/H11B DLPS/ H1B PS  /' head > head_tmp
    sed 's/H12A DLPS/ H2A PS  /' head_tmp > head
    sed 's/ HN1 DLPS/HN1A PS  /' head > head_tmp
    sed 's/ HN2 DLPS/HN1B PS  /' head_tmp > head
    sed 's/ HN3 DLPS/HN1C PS  /' head > head_tmp
    mv head_tmp head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DLPG
if ( $dlpg_found == 1 ) then
  #DLPG is 99 atoms per lipid. 
  set atoms_per_lipid = 99

  echo " Processing: DLPG"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DLPG_tmp | awk '{print$1}'`
  echo "  DLPG Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DLPG Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DLPG Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DLPG_tmp
  rm -f ../DLPG_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DLPG Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DLPG" $i > tail1
    grep -A30 "C33 DLPG" $i >> tail1

    grep -A2 "C22 DLPG" $i > tail2
    grep -A30 "C23 DLPG" $i >> tail2

    head -n25 $i > head
    grep -A5 "C3  DLPG" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DLPG/ C12 LA  /' tail1 >tail1_tmp
    sed 's/ C33 DLPG/ C13 LA  /' tail1_tmp >tail1
    sed 's/ C34 DLPG/ C14 LA  /' tail1 >tail1_tmp
    sed 's/ C35 DLPG/ C15 LA  /' tail1_tmp >tail1
    sed 's/ C36 DLPG/ C16 LA  /' tail1 >tail1_tmp
    sed 's/ C37 DLPG/ C17 LA  /' tail1_tmp >tail1
    sed 's/ C38 DLPG/ C18 LA  /' tail1 >tail1_tmp
    sed 's/ C39 DLPG/ C19 LA  /' tail1_tmp >tail1
    sed 's/C310 DLPG/C110 LA  /' tail1 >tail1_tmp
    sed 's/C311 DLPG/C111 LA  /' tail1_tmp >tail1
    sed 's/C312 DLPG/C112 LA  /' tail1 >tail1_tmp
    sed 's/ H2X DLPG/ H2R LA  /' tail1_tmp >tail1
    sed 's/ H2Y DLPG/ H2S LA  /' tail1 >tail1_tmp
    sed 's/ H3X DLPG/ H3R LA  /' tail1_tmp >tail1
    sed 's/ H3Y DLPG/ H3S LA  /' tail1 >tail1_tmp
    sed 's/ H4X DLPG/ H4R LA  /' tail1_tmp >tail1
    sed 's/ H4Y DLPG/ H4S LA  /' tail1 >tail1_tmp
    sed 's/ H5X DLPG/ H5R LA  /' tail1_tmp >tail1
    sed 's/ H5Y DLPG/ H5S LA  /' tail1 >tail1_tmp
    sed 's/ H6X DLPG/ H6R LA  /' tail1_tmp >tail1
    sed 's/ H6Y DLPG/ H6S LA  /' tail1 >tail1_tmp
    sed 's/ H7X DLPG/ H7R LA  /' tail1_tmp >tail1
    sed 's/ H7Y DLPG/ H7S LA  /' tail1 >tail1_tmp
    sed 's/ H8X DLPG/ H8R LA  /' tail1_tmp >tail1
    sed 's/ H8Y DLPG/ H8S LA  /' tail1 >tail1_tmp
    sed 's/ H9X DLPG/ H9R LA  /' tail1_tmp >tail1
    sed 's/ H9Y DLPG/ H9S LA  /' tail1 >tail1_tmp
    sed 's/H10X DLPG/H10R LA  /' tail1_tmp >tail1
    sed 's/H10Y DLPG/H10S LA  /' tail1 >tail1_tmp
    sed 's/H11X DLPG/H11R LA  /' tail1_tmp >tail1
    sed 's/H11Y DLPG/H11S LA  /' tail1 >tail1_tmp
    sed 's/H12X DLPG/H12R LA  /' tail1_tmp >tail1
    sed 's/H12Y DLPG/H12S LA  /' tail1 >tail1_tmp
    sed 's/H12Z DLPG/H12T LA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DLPG/ C12 LA  /' tail2 >tail2_tmp
    sed 's/ C23 DLPG/ C13 LA  /' tail2_tmp >tail2
    sed 's/ C24 DLPG/ C14 LA  /' tail2 >tail2_tmp
    sed 's/ C25 DLPG/ C15 LA  /' tail2_tmp >tail2
    sed 's/ C26 DLPG/ C16 LA  /' tail2 >tail2_tmp
    sed 's/ C27 DLPG/ C17 LA  /' tail2_tmp >tail2
    sed 's/ C28 DLPG/ C18 LA  /' tail2 >tail2_tmp
    sed 's/ C29 DLPG/ C19 LA  /' tail2_tmp >tail2
    sed 's/C210 DLPG/C110 LA  /' tail2 >tail2_tmp
    sed 's/C211 DLPG/C111 LA  /' tail2_tmp >tail2
    sed 's/C212 DLPG/C112 LA  /' tail2 >tail2_tmp
    sed 's/ H2R DLPG/ H2R LA  /' tail2_tmp >tail2
    sed 's/ H2S DLPG/ H2S LA  /' tail2 >tail2_tmp
    sed 's/ H3R DLPG/ H3R LA  /' tail2_tmp >tail2
    sed 's/ H3S DLPG/ H3S LA  /' tail2 >tail2_tmp
    sed 's/ H4R DLPG/ H4R LA  /' tail2_tmp >tail2
    sed 's/ H4S DLPG/ H4S LA  /' tail2 >tail2_tmp
    sed 's/ H5R DLPG/ H5R LA  /' tail2_tmp >tail2
    sed 's/ H5S DLPG/ H5S LA  /' tail2 >tail2_tmp
    sed 's/ H6R DLPG/ H6R LA  /' tail2_tmp >tail2
    sed 's/ H6S DLPG/ H6S LA  /' tail2 >tail2_tmp
    sed 's/ H7R DLPG/ H7R LA  /' tail2_tmp >tail2
    sed 's/ H7S DLPG/ H7S LA  /' tail2 >tail2_tmp
    sed 's/ H8R DLPG/ H8R LA  /' tail2_tmp >tail2
    sed 's/ H8S DLPG/ H8S LA  /' tail2 >tail2_tmp
    sed 's/ H9R DLPG/ H9R LA  /' tail2_tmp >tail2
    sed 's/ H9S DLPG/ H9S LA  /' tail2 >tail2_tmp
    sed 's/H10R DLPG/H10R LA  /' tail2_tmp >tail2
    sed 's/H10S DLPG/H10S LA  /' tail2 >tail2_tmp
    sed 's/H11R DLPG/H11R LA  /' tail2_tmp >tail2
    sed 's/H11S DLPG/H11S LA  /' tail2 >tail2_tmp
    sed 's/H12R DLPG/H12R LA  /' tail2_tmp >tail2
    sed 's/H12S DLPG/H12S LA  /' tail2 >tail2_tmp
    sed 's/H12T DLPG/H12T LA  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ C13 DLPG/ C33 PGR /' head > head_tmp
    sed 's/H13A DLPG/ H3A PGR /' head_tmp > head
    sed 's/H13B DLPG/ H3B PGR /' head > head_tmp
    sed 's/ OC3 DLPG/ O36 PGR /' head_tmp > head
    sed 's/ HO3 DLPG/HO6A PGR /' head > head_tmp
    sed 's/ C12 DLPG/ C32 PGR /' head_tmp > head
    sed 's/H12A DLPG/ H2A PGR /' head > head_tmp
    sed 's/ OC2 DLPG/ O35 PGR /' head_tmp > head
    sed 's/ HO2 DLPG/HO5A PGR /' head > head_tmp
    sed 's/ C11 DLPG/ C31 PGR /' head_tmp > head
    sed 's/H11A DLPG/ H1A PGR /' head > head_tmp
    sed 's/H11B DLPG/ H1B PGR /' head_tmp > head
    sed 's/ P   DLPG/ P31 PGR /' head > head_tmp
    sed 's/ O13 DLPG/ O33 PGR /' head_tmp > head
    sed 's/ O14 DLPG/ O34 PGR /' head > head_tmp
    sed 's/ O12 DLPG/ O32 PGR /' head_tmp > head
    sed 's/ O11 DLPG/ O31 PGR /' head > head_tmp
    sed 's/ C1  DLPG/ C3  PGR /' head_tmp > head
    sed 's/ HA  DLPG/ HA  PGR /' head > head_tmp
    sed 's/ HB  DLPG/ HB  PGR /' head_tmp > head
    sed 's/ C2  DLPG/ C2  PGR /' head > head_tmp
    sed 's/ HS  DLPG/ HX  PGR /' head_tmp > head
    sed 's/ O21 DLPG/ O21 PGR /' head > head_tmp
    sed 's/ C21 DLPG/ C21 PGR /' head_tmp > head
    sed 's/ O22 DLPG/ O22 PGR /' head > head_tmp
    sed 's/ C3  DLPG/ C1  PGR /' head_tmp > head
    sed 's/ HX  DLPG/ HR  PGR /' head > head_tmp
    sed 's/ HY  DLPG/ HS  PGR /' head_tmp > head
    sed 's/ O31 DLPG/ O11 PGR /' head > head_tmp
    sed 's/ C31 DLPG/ C11 PGR /' head_tmp > head
    sed 's/ O32 DLPG/ O12 PGR /' head > head_tmp
    mv head_tmp head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DLPA
if ( $dlpa_found == 1 ) then
  #DLPA is 88 atoms per lipid. 
  set atoms_per_lipid = 88

  echo " Processing: DLPA"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DLPA_tmp | awk '{print$1}'`
  echo "  DLPA Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DLPA Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DLPA Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DLPA_tmp
  rm -f ../DLPA_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DLPA Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DLPA" $i > tail1
    grep -A30 "C33 DLPA" $i >> tail1

    grep -A2 "C22 DLPA" $i > tail2
    grep -A30 "C23 DLPA" $i >> tail2

    head -n14 $i > head
    grep -A5 "C3  DLPA" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DLPA/ C12 LA  /' tail1 >tail1_tmp
    sed 's/ C33 DLPA/ C13 LA  /' tail1_tmp >tail1
    sed 's/ C34 DLPA/ C14 LA  /' tail1 >tail1_tmp
    sed 's/ C35 DLPA/ C15 LA  /' tail1_tmp >tail1
    sed 's/ C36 DLPA/ C16 LA  /' tail1 >tail1_tmp
    sed 's/ C37 DLPA/ C17 LA  /' tail1_tmp >tail1
    sed 's/ C38 DLPA/ C18 LA  /' tail1 >tail1_tmp
    sed 's/ C39 DLPA/ C19 LA  /' tail1_tmp >tail1
    sed 's/C310 DLPA/C110 LA  /' tail1 >tail1_tmp
    sed 's/C311 DLPA/C111 LA  /' tail1_tmp >tail1
    sed 's/C312 DLPA/C112 LA  /' tail1 >tail1_tmp
    sed 's/ H2X DLPA/ H2R LA  /' tail1_tmp >tail1
    sed 's/ H2Y DLPA/ H2S LA  /' tail1 >tail1_tmp
    sed 's/ H3X DLPA/ H3R LA  /' tail1_tmp >tail1
    sed 's/ H3Y DLPA/ H3S LA  /' tail1 >tail1_tmp
    sed 's/ H4X DLPA/ H4R LA  /' tail1_tmp >tail1
    sed 's/ H4Y DLPA/ H4S LA  /' tail1 >tail1_tmp
    sed 's/ H5X DLPA/ H5R LA  /' tail1_tmp >tail1
    sed 's/ H5Y DLPA/ H5S LA  /' tail1 >tail1_tmp
    sed 's/ H6X DLPA/ H6R LA  /' tail1_tmp >tail1
    sed 's/ H6Y DLPA/ H6S LA  /' tail1 >tail1_tmp
    sed 's/ H7X DLPA/ H7R LA  /' tail1_tmp >tail1
    sed 's/ H7Y DLPA/ H7S LA  /' tail1 >tail1_tmp
    sed 's/ H8X DLPA/ H8R LA  /' tail1_tmp >tail1
    sed 's/ H8Y DLPA/ H8S LA  /' tail1 >tail1_tmp
    sed 's/ H9X DLPA/ H9R LA  /' tail1_tmp >tail1
    sed 's/ H9Y DLPA/ H9S LA  /' tail1 >tail1_tmp
    sed 's/H10X DLPA/H10R LA  /' tail1_tmp >tail1
    sed 's/H10Y DLPA/H10S LA  /' tail1 >tail1_tmp
    sed 's/H11X DLPA/H11R LA  /' tail1_tmp >tail1
    sed 's/H11Y DLPA/H11S LA  /' tail1 >tail1_tmp
    sed 's/H12X DLPA/H12R LA  /' tail1_tmp >tail1
    sed 's/H12Y DLPA/H12S LA  /' tail1 >tail1_tmp
    sed 's/H12Z DLPA/H12T LA  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DLPA/ C12 LA  /' tail2 >tail2_tmp
    sed 's/ C23 DLPA/ C13 LA  /' tail2_tmp >tail2
    sed 's/ C24 DLPA/ C14 LA  /' tail2 >tail2_tmp
    sed 's/ C25 DLPA/ C15 LA  /' tail2_tmp >tail2
    sed 's/ C26 DLPA/ C16 LA  /' tail2 >tail2_tmp
    sed 's/ C27 DLPA/ C17 LA  /' tail2_tmp >tail2
    sed 's/ C28 DLPA/ C18 LA  /' tail2 >tail2_tmp
    sed 's/ C29 DLPA/ C19 LA  /' tail2_tmp >tail2
    sed 's/C210 DLPA/C110 LA  /' tail2 >tail2_tmp
    sed 's/C211 DLPA/C111 LA  /' tail2_tmp >tail2
    sed 's/C212 DLPA/C112 LA  /' tail2 >tail2_tmp
    sed 's/ H2R DLPA/ H2R LA  /' tail2_tmp >tail2
    sed 's/ H2S DLPA/ H2S LA  /' tail2 >tail2_tmp
    sed 's/ H3R DLPA/ H3R LA  /' tail2_tmp >tail2
    sed 's/ H3S DLPA/ H3S LA  /' tail2 >tail2_tmp
    sed 's/ H4R DLPA/ H4R LA  /' tail2_tmp >tail2
    sed 's/ H4S DLPA/ H4S LA  /' tail2 >tail2_tmp
    sed 's/ H5R DLPA/ H5R LA  /' tail2_tmp >tail2
    sed 's/ H5S DLPA/ H5S LA  /' tail2 >tail2_tmp
    sed 's/ H6R DLPA/ H6R LA  /' tail2_tmp >tail2
    sed 's/ H6S DLPA/ H6S LA  /' tail2 >tail2_tmp
    sed 's/ H7R DLPA/ H7R LA  /' tail2_tmp >tail2
    sed 's/ H7S DLPA/ H7S LA  /' tail2 >tail2_tmp
    sed 's/ H8R DLPA/ H8R LA  /' tail2_tmp >tail2
    sed 's/ H8S DLPA/ H8S LA  /' tail2 >tail2_tmp
    sed 's/ H9R DLPA/ H9R LA  /' tail2_tmp >tail2
    sed 's/ H9S DLPA/ H9S LA  /' tail2 >tail2_tmp
    sed 's/H10R DLPA/H10R LA  /' tail2_tmp >tail2
    sed 's/H10S DLPA/H10S LA  /' tail2 >tail2_tmp
    sed 's/H11R DLPA/H11R LA  /' tail2_tmp >tail2
    sed 's/H11S DLPA/H11S LA  /' tail2 >tail2_tmp
    sed 's/H12R DLPA/H12R LA  /' tail2_tmp >tail2
    sed 's/H12S DLPA/H12S LA  /' tail2 >tail2_tmp
    sed 's/H12T DLPA/H12T LA  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ P   DLPA/ P31 PH- /' head > head_tmp
    sed 's/ O13 DLPA/ O33 PH- /' head_tmp > head
    sed 's/ O14 DLPA/ O34 PH- /' head > head_tmp
    sed 's/ O12 DLPA/ O32 PH- /' head_tmp > head
    sed 's/ H12 DLPA/HO2A PH- /' head > head_tmp
    sed 's/ O11 DLPA/ O31 PH- /' head_tmp > head
    sed 's/ C1  DLPA/ C3  PH- /' head > head_tmp
    sed 's/ HA  DLPA/ HA  PH- /' head_tmp > head
    sed 's/ HB  DLPA/ HB  PH- /' head > head_tmp
    sed 's/ C2  DLPA/ C2  PH- /' head_tmp > head
    sed 's/ HS  DLPA/ HX  PH- /' head > head_tmp
    sed 's/ O21 DLPA/ O21 PH- /' head_tmp > head
    sed 's/ C21 DLPA/ C21 PH- /' head > head_tmp
    sed 's/ O22 DLPA/ O22 PH- /' head_tmp > head
    sed 's/ C3  DLPA/ C1  PH- /' head > head_tmp
    sed 's/ HX  DLPA/ HR  PH- /' head_tmp > head
    sed 's/ HY  DLPA/ HS  PH- /' head > head_tmp
    sed 's/ O31 DLPA/ O11 PH- /' head_tmp > head
    sed 's/ C31 DLPA/ C11 PH- /' head > head_tmp
    sed 's/ O32 DLPA/ O12 PH- /' head_tmp > head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DMPC
if ( $dmpc_found == 1 ) then
  #DMPC is 118 atoms per lipid. 
  set atoms_per_lipid = 118

  echo " Processing: DMPC"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DMPC_tmp | awk '{print$1}'`
  echo "  DMPC Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DMPC Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DMPC Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DMPC_tmp
  rm -f ../DMPC_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DMPC Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DMPC" $i > tail1
    grep -A36 "C33 DMPC" $i >> tail1

    grep -A2 "C22 DMPC" $i > tail2
    grep -A36 "C23 DMPC" $i >> tail2

    head -n32 $i > head
    grep -A5 "C3  DMPC" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DMPC/ C12 MY  /' tail1 >tail1_tmp
    sed 's/ C33 DMPC/ C13 MY  /' tail1_tmp >tail1
    sed 's/ C34 DMPC/ C14 MY  /' tail1 >tail1_tmp
    sed 's/ C35 DMPC/ C15 MY  /' tail1_tmp >tail1
    sed 's/ C36 DMPC/ C16 MY  /' tail1 >tail1_tmp
    sed 's/ C37 DMPC/ C17 MY  /' tail1_tmp >tail1
    sed 's/ C38 DMPC/ C18 MY  /' tail1 >tail1_tmp
    sed 's/ C39 DMPC/ C19 MY  /' tail1_tmp >tail1
    sed 's/C310 DMPC/C110 MY  /' tail1 >tail1_tmp
    sed 's/C311 DMPC/C111 MY  /' tail1_tmp >tail1
    sed 's/C312 DMPC/C112 MY  /' tail1 >tail1_tmp
    sed 's/C313 DMPC/C113 MY  /' tail1_tmp >tail1
    sed 's/C314 DMPC/C114 MY  /' tail1 >tail1_tmp
    sed 's/ H2X DMPC/ H2R MY  /' tail1_tmp >tail1
    sed 's/ H2Y DMPC/ H2S MY  /' tail1 >tail1_tmp
    sed 's/ H3X DMPC/ H3R MY  /' tail1_tmp >tail1
    sed 's/ H3Y DMPC/ H3S MY  /' tail1 >tail1_tmp
    sed 's/ H4X DMPC/ H4R MY  /' tail1_tmp >tail1
    sed 's/ H4Y DMPC/ H4S MY  /' tail1 >tail1_tmp
    sed 's/ H5X DMPC/ H5R MY  /' tail1_tmp >tail1
    sed 's/ H5Y DMPC/ H5S MY  /' tail1 >tail1_tmp
    sed 's/ H6X DMPC/ H6R MY  /' tail1_tmp >tail1
    sed 's/ H6Y DMPC/ H6S MY  /' tail1 >tail1_tmp
    sed 's/ H7X DMPC/ H7R MY  /' tail1_tmp >tail1
    sed 's/ H7Y DMPC/ H7S MY  /' tail1 >tail1_tmp
    sed 's/ H8X DMPC/ H8R MY  /' tail1_tmp >tail1
    sed 's/ H8Y DMPC/ H8S MY  /' tail1 >tail1_tmp
    sed 's/ H9X DMPC/ H9R MY  /' tail1_tmp >tail1
    sed 's/ H9Y DMPC/ H9S MY  /' tail1 >tail1_tmp
    sed 's/H10X DMPC/H10R MY  /' tail1_tmp >tail1
    sed 's/H10Y DMPC/H10S MY  /' tail1 >tail1_tmp
    sed 's/H11X DMPC/H11R MY  /' tail1_tmp >tail1
    sed 's/H11Y DMPC/H11S MY  /' tail1 >tail1_tmp
    sed 's/H12X DMPC/H12R MY  /' tail1_tmp >tail1
    sed 's/H12Y DMPC/H12S MY  /' tail1 >tail1_tmp
    sed 's/H13X DMPC/H13R MY  /' tail1_tmp >tail1
    sed 's/H13Y DMPC/H13S MY  /' tail1 >tail1_tmp
    sed 's/H14X DMPC/H14R MY  /' tail1_tmp >tail1
    sed 's/H14Y DMPC/H14S MY  /' tail1 >tail1_tmp
    sed 's/H14Z DMPC/H14T MY  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DMPC/ C12 MY  /' tail2 >tail2_tmp
    sed 's/ C23 DMPC/ C13 MY  /' tail2_tmp >tail2
    sed 's/ C24 DMPC/ C14 MY  /' tail2 >tail2_tmp
    sed 's/ C25 DMPC/ C15 MY  /' tail2_tmp >tail2
    sed 's/ C26 DMPC/ C16 MY  /' tail2 >tail2_tmp
    sed 's/ C27 DMPC/ C17 MY  /' tail2_tmp >tail2
    sed 's/ C28 DMPC/ C18 MY  /' tail2 >tail2_tmp
    sed 's/ C29 DMPC/ C19 MY  /' tail2_tmp >tail2
    sed 's/C210 DMPC/C110 MY  /' tail2 >tail2_tmp
    sed 's/C211 DMPC/C111 MY  /' tail2_tmp >tail2
    sed 's/C212 DMPC/C112 MY  /' tail2 >tail2_tmp
    sed 's/C213 DMPC/C113 MY  /' tail2_tmp >tail2
    sed 's/C214 DMPC/C114 MY  /' tail2 >tail2_tmp
    sed 's/ H2R DMPC/ H2R MY  /' tail2_tmp >tail2
    sed 's/ H2S DMPC/ H2S MY  /' tail2 >tail2_tmp
    sed 's/ H3R DMPC/ H3R MY  /' tail2_tmp >tail2
    sed 's/ H3S DMPC/ H3S MY  /' tail2 >tail2_tmp
    sed 's/ H4R DMPC/ H4R MY  /' tail2_tmp >tail2
    sed 's/ H4S DMPC/ H4S MY  /' tail2 >tail2_tmp
    sed 's/ H5R DMPC/ H5R MY  /' tail2_tmp >tail2
    sed 's/ H5S DMPC/ H5S MY  /' tail2 >tail2_tmp
    sed 's/ H6R DMPC/ H6R MY  /' tail2_tmp >tail2
    sed 's/ H6S DMPC/ H6S MY  /' tail2 >tail2_tmp
    sed 's/ H7R DMPC/ H7R MY  /' tail2_tmp >tail2
    sed 's/ H7S DMPC/ H7S MY  /' tail2 >tail2_tmp
    sed 's/ H8R DMPC/ H8R MY  /' tail2_tmp >tail2
    sed 's/ H8S DMPC/ H8S MY  /' tail2 >tail2_tmp
    sed 's/ H9R DMPC/ H9R MY  /' tail2_tmp >tail2
    sed 's/ H9S DMPC/ H9S MY  /' tail2 >tail2_tmp
    sed 's/H10R DMPC/H10R MY  /' tail2_tmp >tail2
    sed 's/H10S DMPC/H10S MY  /' tail2 >tail2_tmp
    sed 's/H11R DMPC/H11R MY  /' tail2_tmp >tail2
    sed 's/H11S DMPC/H11S MY  /' tail2 >tail2_tmp
    sed 's/H12R DMPC/H12R MY  /' tail2_tmp >tail2
    sed 's/H12S DMPC/H12S MY  /' tail2 >tail2_tmp
    sed 's/H13R DMPC/H13R MY  /' tail2_tmp >tail2
    sed 's/H13S DMPC/H13S MY  /' tail2 >tail2_tmp
    sed 's/H14R DMPC/H14R MY  /' tail2_tmp >tail2
    sed 's/H14S DMPC/H14S MY  /' tail2 >tail2_tmp
    sed 's/H14T DMPC/H14T MY  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DMPC/ N31 PC  /' head > head_tmp
    sed 's/ C13 DMPC/ C33 PC  /' head_tmp > head
    sed 's/ C14 DMPC/ C34 PC  /' head > head_tmp
    sed 's/ C15 DMPC/ C35 PC  /' head_tmp > head
    sed 's/ C12 DMPC/ C32 PC  /' head > head_tmp
    sed 's/ C11 DMPC/ C31 PC  /' head_tmp > head
    sed 's/ O12 DMPC/ O32 PC  /' head > head_tmp
    sed 's/ P   DMPC/ P31 PC  /' head_tmp > head
    sed 's/ O11 DMPC/ O31 PC  /' head > head_tmp
    sed 's/ O13 DMPC/ O33 PC  /' head_tmp > head
    sed 's/ O14 DMPC/ O34 PC  /' head > head_tmp
    sed 's/ C1  DMPC/ C3  PC  /' head_tmp > head
    sed 's/ C2  DMPC/ C2  PC  /' head > head_tmp
    sed 's/ C3  DMPC/ C1  PC  /' head_tmp > head
    sed 's/ O31 DMPC/ O11 PC  /' head > head_tmp
    sed 's/ C31 DMPC/ C11 PC  /' head_tmp > head
    sed 's/ O32 DMPC/ O12 PC  /' head > head_tmp
    sed 's/ O21 DMPC/ O21 PC  /' head_tmp > head
    sed 's/ C21 DMPC/ C21 PC  /' head > head_tmp
    sed 's/ O22 DMPC/ O22 PC  /' head_tmp > head
    sed 's/ HX  DMPC/ HR  PC  /' head > head_tmp
    sed 's/ HY  DMPC/ HS  PC  /' head_tmp > head
    sed 's/ HS  DMPC/ HX  PC  /' head > head_tmp
    sed 's/ HA  DMPC/ HA  PC  /' head_tmp > head
    sed 's/ HB  DMPC/ HB  PC  /' head > head_tmp
    sed 's/H11A DMPC/ H1A PC  /' head_tmp > head
    sed 's/H11B DMPC/ H1B PC  /' head > head_tmp
    sed 's/H12A DMPC/ H2A PC  /' head_tmp > head
    sed 's/H12B DMPC/ H2B PC  /' head > head_tmp
    sed 's/H13A DMPC/ H3A PC  /' head_tmp > head
    sed 's/H13B DMPC/ H3B PC  /' head > head_tmp
    sed 's/H13C DMPC/ H3C PC  /' head_tmp > head
    sed 's/H14A DMPC/ H4A PC  /' head > head_tmp
    sed 's/H14B DMPC/ H4B PC  /' head_tmp > head
    sed 's/H14C DMPC/ H4C PC  /' head > head_tmp
    sed 's/H15A DMPC/ H5A PC  /' head_tmp > head
    sed 's/H15B DMPC/ H5B PC  /' head > head_tmp
    sed 's/H15C DMPC/ H5C PC  /' head_tmp > head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DMPE
if ( $dmpe_found == 1 ) then
  #DMPE is 109 atoms per lipid. 
  set atoms_per_lipid = 109

  echo " Processing: DMPE"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DMPE_tmp | awk '{print$1}'`
  echo "  DMPE Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DMPE Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DMPE Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DMPE_tmp
  rm -f ../DMPE_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DMPE Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DMPE" $i > tail1
    grep -A36 "C33 DMPE" $i >> tail1

    grep -A2 "C22 DMPE" $i > tail2
    grep -A36 "C23 DMPE" $i >> tail2

    head -n23 $i > head
    grep -A5 "C3  DMPE" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DMPE/ C12 MY  /' tail1 >tail1_tmp
    sed 's/ C33 DMPE/ C13 MY  /' tail1_tmp >tail1
    sed 's/ C34 DMPE/ C14 MY  /' tail1 >tail1_tmp
    sed 's/ C35 DMPE/ C15 MY  /' tail1_tmp >tail1
    sed 's/ C36 DMPE/ C16 MY  /' tail1 >tail1_tmp
    sed 's/ C37 DMPE/ C17 MY  /' tail1_tmp >tail1
    sed 's/ C38 DMPE/ C18 MY  /' tail1 >tail1_tmp
    sed 's/ C39 DMPE/ C19 MY  /' tail1_tmp >tail1
    sed 's/C310 DMPE/C110 MY  /' tail1 >tail1_tmp
    sed 's/C311 DMPE/C111 MY  /' tail1_tmp >tail1
    sed 's/C312 DMPE/C112 MY  /' tail1 >tail1_tmp
    sed 's/C313 DMPE/C113 MY  /' tail1_tmp >tail1
    sed 's/C314 DMPE/C114 MY  /' tail1 >tail1_tmp
    sed 's/ H2X DMPE/ H2R MY  /' tail1_tmp >tail1
    sed 's/ H2Y DMPE/ H2S MY  /' tail1 >tail1_tmp
    sed 's/ H3X DMPE/ H3R MY  /' tail1_tmp >tail1
    sed 's/ H3Y DMPE/ H3S MY  /' tail1 >tail1_tmp
    sed 's/ H4X DMPE/ H4R MY  /' tail1_tmp >tail1
    sed 's/ H4Y DMPE/ H4S MY  /' tail1 >tail1_tmp
    sed 's/ H5X DMPE/ H5R MY  /' tail1_tmp >tail1
    sed 's/ H5Y DMPE/ H5S MY  /' tail1 >tail1_tmp
    sed 's/ H6X DMPE/ H6R MY  /' tail1_tmp >tail1
    sed 's/ H6Y DMPE/ H6S MY  /' tail1 >tail1_tmp
    sed 's/ H7X DMPE/ H7R MY  /' tail1_tmp >tail1
    sed 's/ H7Y DMPE/ H7S MY  /' tail1 >tail1_tmp
    sed 's/ H8X DMPE/ H8R MY  /' tail1_tmp >tail1
    sed 's/ H8Y DMPE/ H8S MY  /' tail1 >tail1_tmp
    sed 's/ H9X DMPE/ H9R MY  /' tail1_tmp >tail1
    sed 's/ H9Y DMPE/ H9S MY  /' tail1 >tail1_tmp
    sed 's/H10X DMPE/H10R MY  /' tail1_tmp >tail1
    sed 's/H10Y DMPE/H10S MY  /' tail1 >tail1_tmp
    sed 's/H11X DMPE/H11R MY  /' tail1_tmp >tail1
    sed 's/H11Y DMPE/H11S MY  /' tail1 >tail1_tmp
    sed 's/H12X DMPE/H12R MY  /' tail1_tmp >tail1
    sed 's/H12Y DMPE/H12S MY  /' tail1 >tail1_tmp
    sed 's/H13X DMPE/H13R MY  /' tail1_tmp >tail1
    sed 's/H13Y DMPE/H13S MY  /' tail1 >tail1_tmp
    sed 's/H14X DMPE/H14R MY  /' tail1_tmp >tail1
    sed 's/H14Y DMPE/H14S MY  /' tail1 >tail1_tmp
    sed 's/H14Z DMPE/H14T MY  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DMPE/ C12 MY  /' tail2 >tail2_tmp
    sed 's/ C23 DMPE/ C13 MY  /' tail2_tmp >tail2
    sed 's/ C24 DMPE/ C14 MY  /' tail2 >tail2_tmp
    sed 's/ C25 DMPE/ C15 MY  /' tail2_tmp >tail2
    sed 's/ C26 DMPE/ C16 MY  /' tail2 >tail2_tmp
    sed 's/ C27 DMPE/ C17 MY  /' tail2_tmp >tail2
    sed 's/ C28 DMPE/ C18 MY  /' tail2 >tail2_tmp
    sed 's/ C29 DMPE/ C19 MY  /' tail2_tmp >tail2
    sed 's/C210 DMPE/C110 MY  /' tail2 >tail2_tmp
    sed 's/C211 DMPE/C111 MY  /' tail2_tmp >tail2
    sed 's/C212 DMPE/C112 MY  /' tail2 >tail2_tmp
    sed 's/C213 DMPE/C113 MY  /' tail2_tmp >tail2
    sed 's/C214 DMPE/C114 MY  /' tail2 >tail2_tmp
    sed 's/ H2R DMPE/ H2R MY  /' tail2_tmp >tail2
    sed 's/ H2S DMPE/ H2S MY  /' tail2 >tail2_tmp
    sed 's/ H3R DMPE/ H3R MY  /' tail2_tmp >tail2
    sed 's/ H3S DMPE/ H3S MY  /' tail2 >tail2_tmp
    sed 's/ H4R DMPE/ H4R MY  /' tail2_tmp >tail2
    sed 's/ H4S DMPE/ H4S MY  /' tail2 >tail2_tmp
    sed 's/ H5R DMPE/ H5R MY  /' tail2_tmp >tail2
    sed 's/ H5S DMPE/ H5S MY  /' tail2 >tail2_tmp
    sed 's/ H6R DMPE/ H6R MY  /' tail2_tmp >tail2
    sed 's/ H6S DMPE/ H6S MY  /' tail2 >tail2_tmp
    sed 's/ H7R DMPE/ H7R MY  /' tail2_tmp >tail2
    sed 's/ H7S DMPE/ H7S MY  /' tail2 >tail2_tmp
    sed 's/ H8R DMPE/ H8R MY  /' tail2_tmp >tail2
    sed 's/ H8S DMPE/ H8S MY  /' tail2 >tail2_tmp
    sed 's/ H9R DMPE/ H9R MY  /' tail2_tmp >tail2
    sed 's/ H9S DMPE/ H9S MY  /' tail2 >tail2_tmp
    sed 's/H10R DMPE/H10R MY  /' tail2_tmp >tail2
    sed 's/H10S DMPE/H10S MY  /' tail2 >tail2_tmp
    sed 's/H11R DMPE/H11R MY  /' tail2_tmp >tail2
    sed 's/H11S DMPE/H11S MY  /' tail2 >tail2_tmp
    sed 's/H12R DMPE/H12R MY  /' tail2_tmp >tail2
    sed 's/H12S DMPE/H12S MY  /' tail2 >tail2_tmp
    sed 's/H13R DMPE/H13R MY  /' tail2_tmp >tail2
    sed 's/H13S DMPE/H13S MY  /' tail2 >tail2_tmp
    sed 's/H14R DMPE/H14R MY  /' tail2_tmp >tail2
    sed 's/H14S DMPE/H14S MY  /' tail2 >tail2_tmp
    sed 's/H14T DMPE/H14T MY  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DMPE/ N31 PE  /' head > head_tmp
    sed 's/ C12 DMPE/ C32 PE  /' head_tmp > head
    sed 's/ C11 DMPE/ C31 PE  /' head > head_tmp
    sed 's/ O12 DMPE/ O32 PE  /' head_tmp > head
    sed 's/ P   DMPE/ P31 PE  /' head > head_tmp
    sed 's/ O11 DMPE/ O31 PE  /' head_tmp > head
    sed 's/ O13 DMPE/ O33 PE  /' head > head_tmp
    sed 's/ O14 DMPE/ O34 PE  /' head_tmp > head
    sed 's/ C1  DMPE/ C3  PE  /' head > head_tmp
    sed 's/ C2  DMPE/ C2  PE  /' head_tmp > head
    sed 's/ C3  DMPE/ C1  PE  /' head > head_tmp
    sed 's/ O31 DMPE/ O11 PE  /' head_tmp > head
    sed 's/ C31 DMPE/ C11 PE  /' head > head_tmp
    sed 's/ O32 DMPE/ O12 PE  /' head_tmp > head
    sed 's/ O21 DMPE/ O21 PE  /' head > head_tmp
    sed 's/ C21 DMPE/ C21 PE  /' head_tmp > head
    sed 's/ O22 DMPE/ O22 PE  /' head > head_tmp
    sed 's/ HX  DMPE/ HR  PE  /' head_tmp > head
    sed 's/ HY  DMPE/ HS  PE  /' head > head_tmp
    sed 's/ HS  DMPE/ HX  PE  /' head_tmp > head
    sed 's/ HA  DMPE/ HA  PE  /' head > head_tmp
    sed 's/ HB  DMPE/ HB  PE  /' head_tmp > head
    sed 's/H11A DMPE/ H1A PE  /' head > head_tmp
    sed 's/H11B DMPE/ H1B PE  /' head_tmp > head
    sed 's/H12A DMPE/ H2A PE  /' head > head_tmp
    sed 's/H12B DMPE/ H2B PE  /' head_tmp > head
    sed 's/ HN1 DMPE/HN1A PE  /' head > head_tmp
    sed 's/ HN2 DMPE/HN1B PE  /' head_tmp > head
    sed 's/ HN3 DMPE/HN1C PE  /' head > head_tmp
    mv head_tmp head  
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DMPS
if ( $dmps_found == 1 ) then
  #DMPS is 111 atoms per lipid. 
  set atoms_per_lipid = 111

  echo " Processing: DMPS"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DMPC_tmp | awk '{print$1}'`
  echo "  DMPS Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DMPS Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DMPS Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DMPS_tmp
  rm -f ../DMPS_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DMPS Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DMPS" $i > tail1
    grep -A36 "C33 DMPS" $i >> tail1

    grep -A2 "C22 DMPS" $i > tail2
    grep -A36 "C23 DMPS" $i >> tail2

    head -n25 $i > head
    grep -A5 "C3  DMPS" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DMPS/ C12 MY  /' tail1 >tail1_tmp
    sed 's/ C33 DMPS/ C13 MY  /' tail1_tmp >tail1
    sed 's/ C34 DMPS/ C14 MY  /' tail1 >tail1_tmp
    sed 's/ C35 DMPS/ C15 MY  /' tail1_tmp >tail1
    sed 's/ C36 DMPS/ C16 MY  /' tail1 >tail1_tmp
    sed 's/ C37 DMPS/ C17 MY  /' tail1_tmp >tail1
    sed 's/ C38 DMPS/ C18 MY  /' tail1 >tail1_tmp
    sed 's/ C39 DMPS/ C19 MY  /' tail1_tmp >tail1
    sed 's/C310 DMPS/C110 MY  /' tail1 >tail1_tmp
    sed 's/C311 DMPS/C111 MY  /' tail1_tmp >tail1
    sed 's/C312 DMPS/C112 MY  /' tail1 >tail1_tmp
    sed 's/C313 DMPS/C113 MY  /' tail1_tmp >tail1
    sed 's/C314 DMPS/C114 MY  /' tail1 >tail1_tmp
    sed 's/ H2X DMPS/ H2R MY  /' tail1_tmp >tail1
    sed 's/ H2Y DMPS/ H2S MY  /' tail1 >tail1_tmp
    sed 's/ H3X DMPS/ H3R MY  /' tail1_tmp >tail1
    sed 's/ H3Y DMPS/ H3S MY  /' tail1 >tail1_tmp
    sed 's/ H4X DMPS/ H4R MY  /' tail1_tmp >tail1
    sed 's/ H4Y DMPS/ H4S MY  /' tail1 >tail1_tmp
    sed 's/ H5X DMPS/ H5R MY  /' tail1_tmp >tail1
    sed 's/ H5Y DMPS/ H5S MY  /' tail1 >tail1_tmp
    sed 's/ H6X DMPS/ H6R MY  /' tail1_tmp >tail1
    sed 's/ H6Y DMPS/ H6S MY  /' tail1 >tail1_tmp
    sed 's/ H7X DMPS/ H7R MY  /' tail1_tmp >tail1
    sed 's/ H7Y DMPS/ H7S MY  /' tail1 >tail1_tmp
    sed 's/ H8X DMPS/ H8R MY  /' tail1_tmp >tail1
    sed 's/ H8Y DMPS/ H8S MY  /' tail1 >tail1_tmp
    sed 's/ H9X DMPS/ H9R MY  /' tail1_tmp >tail1
    sed 's/ H9Y DMPS/ H9S MY  /' tail1 >tail1_tmp
    sed 's/H10X DMPS/H10R MY  /' tail1_tmp >tail1
    sed 's/H10Y DMPS/H10S MY  /' tail1 >tail1_tmp
    sed 's/H11X DMPS/H11R MY  /' tail1_tmp >tail1
    sed 's/H11Y DMPS/H11S MY  /' tail1 >tail1_tmp
    sed 's/H12X DMPS/H12R MY  /' tail1_tmp >tail1
    sed 's/H12Y DMPS/H12S MY  /' tail1 >tail1_tmp
    sed 's/H13X DMPS/H13R MY  /' tail1_tmp >tail1
    sed 's/H13Y DMPS/H13S MY  /' tail1 >tail1_tmp
    sed 's/H14X DMPS/H14R MY  /' tail1_tmp >tail1
    sed 's/H14Y DMPS/H14S MY  /' tail1 >tail1_tmp
    sed 's/H14Z DMPS/H14T MY  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DMPS/ C12 MY  /' tail2 >tail2_tmp
    sed 's/ C23 DMPS/ C13 MY  /' tail2_tmp >tail2
    sed 's/ C24 DMPS/ C14 MY  /' tail2 >tail2_tmp
    sed 's/ C25 DMPS/ C15 MY  /' tail2_tmp >tail2
    sed 's/ C26 DMPS/ C16 MY  /' tail2 >tail2_tmp
    sed 's/ C27 DMPS/ C17 MY  /' tail2_tmp >tail2
    sed 's/ C28 DMPS/ C18 MY  /' tail2 >tail2_tmp
    sed 's/ C29 DMPS/ C19 MY  /' tail2_tmp >tail2
    sed 's/C210 DMPS/C110 MY  /' tail2 >tail2_tmp
    sed 's/C211 DMPS/C111 MY  /' tail2_tmp >tail2
    sed 's/C212 DMPS/C112 MY  /' tail2 >tail2_tmp
    sed 's/C213 DMPS/C113 MY  /' tail2_tmp >tail2
    sed 's/C214 DMPS/C114 MY  /' tail2 >tail2_tmp
    sed 's/ H2R DMPS/ H2R MY  /' tail2_tmp >tail2
    sed 's/ H2S DMPS/ H2S MY  /' tail2 >tail2_tmp
    sed 's/ H3R DMPS/ H3R MY  /' tail2_tmp >tail2
    sed 's/ H3S DMPS/ H3S MY  /' tail2 >tail2_tmp
    sed 's/ H4R DMPS/ H4R MY  /' tail2_tmp >tail2
    sed 's/ H4S DMPS/ H4S MY  /' tail2 >tail2_tmp
    sed 's/ H5R DMPS/ H5R MY  /' tail2_tmp >tail2
    sed 's/ H5S DMPS/ H5S MY  /' tail2 >tail2_tmp
    sed 's/ H6R DMPS/ H6R MY  /' tail2_tmp >tail2
    sed 's/ H6S DMPS/ H6S MY  /' tail2 >tail2_tmp
    sed 's/ H7R DMPS/ H7R MY  /' tail2_tmp >tail2
    sed 's/ H7S DMPS/ H7S MY  /' tail2 >tail2_tmp
    sed 's/ H8R DMPS/ H8R MY  /' tail2_tmp >tail2
    sed 's/ H8S DMPS/ H8S MY  /' tail2 >tail2_tmp
    sed 's/ H9R DMPS/ H9R MY  /' tail2_tmp >tail2
    sed 's/ H9S DMPS/ H9S MY  /' tail2 >tail2_tmp
    sed 's/H10R DMPS/H10R MY  /' tail2_tmp >tail2
    sed 's/H10S DMPS/H10S MY  /' tail2 >tail2_tmp
    sed 's/H11R DMPS/H11R MY  /' tail2_tmp >tail2
    sed 's/H11S DMPS/H11S MY  /' tail2 >tail2_tmp
    sed 's/H12R DMPS/H12R MY  /' tail2_tmp >tail2
    sed 's/H12S DMPS/H12S MY  /' tail2 >tail2_tmp
    sed 's/H13R DMPS/H13R MY  /' tail2_tmp >tail2
    sed 's/H13S DMPS/H13S MY  /' tail2 >tail2_tmp
    sed 's/H14R DMPS/H14R MY  /' tail2_tmp >tail2
    sed 's/H14S DMPS/H14S MY  /' tail2 >tail2_tmp
    sed 's/H14T DMPS/H14T MY  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ N   DMPS/ N31 PS  /' head > head_tmp
    sed 's/ C12 DMPS/ C32 PS  /' head_tmp > head
    sed 's/ C13 DMPS/ C33 PS  /' head > head_tmp
    sed 's/O13A DMPS/ O35 PS  /' head_tmp > head
    sed 's/O13B DMPS/ O36 PS  /' head > head_tmp
    sed 's/ O12 DMPS/ O32 PS  /' head_tmp > head
    sed 's/ P   DMPS/ P31 PS  /' head > head_tmp
    sed 's/ O11 DMPS/ O31 PS  /' head_tmp > head
    sed 's/ O13 DMPS/ O33 PS  /' head > head_tmp
    sed 's/ C11 DMPS/ C31 PS  /' head_tmp > head
    sed 's/ O14 DMPS/ O34 PS  /' head > head_tmp
    sed 's/ C1  DMPS/ C3  PS  /' head_tmp > head
    sed 's/ C2  DMPS/ C2  PS  /' head > head_tmp
    sed 's/ C3  DMPS/ C1  PS  /' head_tmp > head
    sed 's/ O31 DMPS/ O11 PS  /' head > head_tmp
    sed 's/ C31 DMPS/ C11 PS  /' head_tmp > head
    sed 's/ O32 DMPS/ O12 PS  /' head > head_tmp
    sed 's/ O21 DMPS/ O21 PS  /' head_tmp > head
    sed 's/ C21 DMPS/ C21 PS  /' head > head_tmp
    sed 's/ O22 DMPS/ O22 PS  /' head_tmp > head
    sed 's/ HX  DMPS/ HR  PS  /' head > head_tmp
    sed 's/ HY  DMPS/ HS  PS  /' head_tmp > head
    sed 's/ HS  DMPS/ HX  PS  /' head > head_tmp
    sed 's/ HA  DMPS/ HA  PS  /' head_tmp > head
    sed 's/ HB  DMPS/ HB  PS  /' head > head_tmp
    sed 's/H11A DMPS/ H1A PS  /' head_tmp > head
    sed 's/H11B DMPS/ H1B PS  /' head > head_tmp
    sed 's/H12A DMPS/ H2A PS  /' head_tmp > head
    sed 's/ HN1 DMPS/HN1A PS  /' head > head_tmp
    sed 's/ HN2 DMPS/HN1B PS  /' head_tmp > head
    sed 's/ HN3 DMPS/HN1C PS  /' head > head_tmp
    mv head_tmp head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DMPG
if ( $dmpg_found == 1 ) then
  #DMPG is 111 atoms per lipid. 
  set atoms_per_lipid = 111

  echo " Processing: DMPG"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DMPG_tmp | awk '{print$1}'`
  echo "  DMPG Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DMPG Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DMPG Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DMPG_tmp
  rm -f ../DMPG_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DMPG Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DMPG" $i > tail1
    grep -A36 "C33 DMPG" $i >> tail1

    grep -A2 "C22 DMPG" $i > tail2
    grep -A36 "C23 DMPG" $i >> tail2

    head -n25 $i > head
    grep -A5 "C3  DMPG" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DMPG/ C12 MY  /' tail1 >tail1_tmp
    sed 's/ C33 DMPG/ C13 MY  /' tail1_tmp >tail1
    sed 's/ C34 DMPG/ C14 MY  /' tail1 >tail1_tmp
    sed 's/ C35 DMPG/ C15 MY  /' tail1_tmp >tail1
    sed 's/ C36 DMPG/ C16 MY  /' tail1 >tail1_tmp
    sed 's/ C37 DMPG/ C17 MY  /' tail1_tmp >tail1
    sed 's/ C38 DMPG/ C18 MY  /' tail1 >tail1_tmp
    sed 's/ C39 DMPG/ C19 MY  /' tail1_tmp >tail1
    sed 's/C310 DMPG/C110 MY  /' tail1 >tail1_tmp
    sed 's/C311 DMPG/C111 MY  /' tail1_tmp >tail1
    sed 's/C312 DMPG/C112 MY  /' tail1 >tail1_tmp
    sed 's/C313 DMPG/C113 MY  /' tail1_tmp >tail1
    sed 's/C314 DMPG/C114 MY  /' tail1 >tail1_tmp
    sed 's/ H2X DMPG/ H2R MY  /' tail1_tmp >tail1
    sed 's/ H2Y DMPG/ H2S MY  /' tail1 >tail1_tmp
    sed 's/ H3X DMPG/ H3R MY  /' tail1_tmp >tail1
    sed 's/ H3Y DMPG/ H3S MY  /' tail1 >tail1_tmp
    sed 's/ H4X DMPG/ H4R MY  /' tail1_tmp >tail1
    sed 's/ H4Y DMPG/ H4S MY  /' tail1 >tail1_tmp
    sed 's/ H5X DMPG/ H5R MY  /' tail1_tmp >tail1
    sed 's/ H5Y DMPG/ H5S MY  /' tail1 >tail1_tmp
    sed 's/ H6X DMPG/ H6R MY  /' tail1_tmp >tail1
    sed 's/ H6Y DMPG/ H6S MY  /' tail1 >tail1_tmp
    sed 's/ H7X DMPG/ H7R MY  /' tail1_tmp >tail1
    sed 's/ H7Y DMPG/ H7S MY  /' tail1 >tail1_tmp
    sed 's/ H8X DMPG/ H8R MY  /' tail1_tmp >tail1
    sed 's/ H8Y DMPG/ H8S MY  /' tail1 >tail1_tmp
    sed 's/ H9X DMPG/ H9R MY  /' tail1_tmp >tail1
    sed 's/ H9Y DMPG/ H9S MY  /' tail1 >tail1_tmp
    sed 's/H10X DMPG/H10R MY  /' tail1_tmp >tail1
    sed 's/H10Y DMPG/H10S MY  /' tail1 >tail1_tmp
    sed 's/H11X DMPG/H11R MY  /' tail1_tmp >tail1
    sed 's/H11Y DMPG/H11S MY  /' tail1 >tail1_tmp
    sed 's/H12X DMPG/H12R MY  /' tail1_tmp >tail1
    sed 's/H12Y DMPG/H12S MY  /' tail1 >tail1_tmp
    sed 's/H13X DMPG/H13R MY  /' tail1_tmp >tail1
    sed 's/H13Y DMPG/H13S MY  /' tail1 >tail1_tmp
    sed 's/H14X DMPG/H14R MY  /' tail1_tmp >tail1
    sed 's/H14Y DMPG/H14S MY  /' tail1 >tail1_tmp
    sed 's/H14Z DMPG/H14T MY  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DMPG/ C12 MY  /' tail2 >tail2_tmp
    sed 's/ C23 DMPG/ C13 MY  /' tail2_tmp >tail2
    sed 's/ C24 DMPG/ C14 MY  /' tail2 >tail2_tmp
    sed 's/ C25 DMPG/ C15 MY  /' tail2_tmp >tail2
    sed 's/ C26 DMPG/ C16 MY  /' tail2 >tail2_tmp
    sed 's/ C27 DMPG/ C17 MY  /' tail2_tmp >tail2
    sed 's/ C28 DMPG/ C18 MY  /' tail2 >tail2_tmp
    sed 's/ C29 DMPG/ C19 MY  /' tail2_tmp >tail2
    sed 's/C210 DMPG/C110 MY  /' tail2 >tail2_tmp
    sed 's/C211 DMPG/C111 MY  /' tail2_tmp >tail2
    sed 's/C212 DMPG/C112 MY  /' tail2 >tail2_tmp
    sed 's/C213 DMPG/C113 MY  /' tail2_tmp >tail2
    sed 's/C214 DMPG/C114 MY  /' tail2 >tail2_tmp
    sed 's/ H2R DMPG/ H2R MY  /' tail2_tmp >tail2
    sed 's/ H2S DMPG/ H2S MY  /' tail2 >tail2_tmp
    sed 's/ H3R DMPG/ H3R MY  /' tail2_tmp >tail2
    sed 's/ H3S DMPG/ H3S MY  /' tail2 >tail2_tmp
    sed 's/ H4R DMPG/ H4R MY  /' tail2_tmp >tail2
    sed 's/ H4S DMPG/ H4S MY  /' tail2 >tail2_tmp
    sed 's/ H5R DMPG/ H5R MY  /' tail2_tmp >tail2
    sed 's/ H5S DMPG/ H5S MY  /' tail2 >tail2_tmp
    sed 's/ H6R DMPG/ H6R MY  /' tail2_tmp >tail2
    sed 's/ H6S DMPG/ H6S MY  /' tail2 >tail2_tmp
    sed 's/ H7R DMPG/ H7R MY  /' tail2_tmp >tail2
    sed 's/ H7S DMPG/ H7S MY  /' tail2 >tail2_tmp
    sed 's/ H8R DMPG/ H8R MY  /' tail2_tmp >tail2
    sed 's/ H8S DMPG/ H8S MY  /' tail2 >tail2_tmp
    sed 's/ H9R DMPG/ H9R MY  /' tail2_tmp >tail2
    sed 's/ H9S DMPG/ H9S MY  /' tail2 >tail2_tmp
    sed 's/H10R DMPG/H10R MY  /' tail2_tmp >tail2
    sed 's/H10S DMPG/H10S MY  /' tail2 >tail2_tmp
    sed 's/H11R DMPG/H11R MY  /' tail2_tmp >tail2
    sed 's/H11S DMPG/H11S MY  /' tail2 >tail2_tmp
    sed 's/H12R DMPG/H12R MY  /' tail2_tmp >tail2
    sed 's/H12S DMPG/H12S MY  /' tail2 >tail2_tmp
    sed 's/H13R DMPG/H13R MY  /' tail2_tmp >tail2
    sed 's/H13S DMPG/H13S MY  /' tail2 >tail2_tmp
    sed 's/H14R DMPG/H14R MY  /' tail2_tmp >tail2
    sed 's/H14S DMPG/H14S MY  /' tail2 >tail2_tmp
    sed 's/H14T DMPG/H14T MY  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ C13 DMPG/ C33 PGR /' head > head_tmp
    sed 's/H13A DMPG/ H3A PGR /' head_tmp > head
    sed 's/H13B DMPG/ H3B PGR /' head > head_tmp
    sed 's/ OC3 DMPG/ O36 PGR /' head_tmp > head
    sed 's/ HO3 DMPG/HO6A PGR /' head > head_tmp
    sed 's/ C12 DMPG/ C32 PGR /' head_tmp > head
    sed 's/H12A DMPG/ H2A PGR /' head > head_tmp
    sed 's/ OC2 DMPG/ O35 PGR /' head_tmp > head
    sed 's/ HO2 DMPG/HO5A PGR /' head > head_tmp
    sed 's/ C11 DMPG/ C31 PGR /' head_tmp > head
    sed 's/H11A DMPG/ H1A PGR /' head > head_tmp
    sed 's/H11B DMPG/ H1B PGR /' head_tmp > head
    sed 's/ P   DMPG/ P31 PGR /' head > head_tmp
    sed 's/ O13 DMPG/ O33 PGR /' head_tmp > head
    sed 's/ O14 DMPG/ O34 PGR /' head > head_tmp
    sed 's/ O12 DMPG/ O32 PGR /' head_tmp > head
    sed 's/ O11 DMPG/ O31 PGR /' head > head_tmp
    sed 's/ C1  DMPG/ C3  PGR /' head_tmp > head
    sed 's/ HA  DMPG/ HA  PGR /' head > head_tmp
    sed 's/ HB  DMPG/ HB  PGR /' head_tmp > head
    sed 's/ C2  DMPG/ C2  PGR /' head > head_tmp
    sed 's/ HS  DMPG/ HX  PGR /' head_tmp > head
    sed 's/ O21 DMPG/ O21 PGR /' head > head_tmp
    sed 's/ C21 DMPG/ C21 PGR /' head_tmp > head
    sed 's/ O22 DMPG/ O22 PGR /' head > head_tmp
    sed 's/ C3  DMPG/ C1  PGR /' head_tmp > head
    sed 's/ HX  DMPG/ HR  PGR /' head > head_tmp
    sed 's/ HY  DMPG/ HS  PGR /' head_tmp > head
    sed 's/ O31 DMPG/ O11 PGR /' head > head_tmp
    sed 's/ C31 DMPG/ C11 PGR /' head_tmp > head
    sed 's/ O32 DMPG/ O12 PGR /' head > head_tmp
    mv head_tmp head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#DMPA
if ( $dmpa_found == 1 ) then
  #DMPA is 100 atoms per lipid. 
  set atoms_per_lipid = 100

  echo " Processing: DMPA"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l DMPA_tmp | awk '{print$1}'`
  echo "  DMPA Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted DMPA Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " DMPA Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../DMPA_tmp
  rm -f ../DMPA_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing DMPA Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the head and two tails.
    grep -A2 "C32 DMPA" $i > tail1
    grep -A36 "C33 DMPA" $i >> tail1

    grep -A2 "C22 DMPA" $i > tail2
    grep -A36 "C23 DMPA" $i >> tail2

    head -n14 $i > head
    grep -A5 "C3  DMPA" $i >> head

    #Modify residue and atom names
    #Tail 1
    sed 's/ C32 DMPA/ C12 MY  /' tail1 >tail1_tmp
    sed 's/ C33 DMPA/ C13 MY  /' tail1_tmp >tail1
    sed 's/ C34 DMPA/ C14 MY  /' tail1 >tail1_tmp
    sed 's/ C35 DMPA/ C15 MY  /' tail1_tmp >tail1
    sed 's/ C36 DMPA/ C16 MY  /' tail1 >tail1_tmp
    sed 's/ C37 DMPA/ C17 MY  /' tail1_tmp >tail1
    sed 's/ C38 DMPA/ C18 MY  /' tail1 >tail1_tmp
    sed 's/ C39 DMPA/ C19 MY  /' tail1_tmp >tail1
    sed 's/C310 DMPA/C110 MY  /' tail1 >tail1_tmp
    sed 's/C311 DMPA/C111 MY  /' tail1_tmp >tail1
    sed 's/C312 DMPA/C112 MY  /' tail1 >tail1_tmp
    sed 's/C313 DMPA/C113 MY  /' tail1_tmp >tail1
    sed 's/C314 DMPA/C114 MY  /' tail1 >tail1_tmp
    sed 's/ H2X DMPA/ H2R MY  /' tail1_tmp >tail1
    sed 's/ H2Y DMPA/ H2S MY  /' tail1 >tail1_tmp
    sed 's/ H3X DMPA/ H3R MY  /' tail1_tmp >tail1
    sed 's/ H3Y DMPA/ H3S MY  /' tail1 >tail1_tmp
    sed 's/ H4X DMPA/ H4R MY  /' tail1_tmp >tail1
    sed 's/ H4Y DMPA/ H4S MY  /' tail1 >tail1_tmp
    sed 's/ H5X DMPA/ H5R MY  /' tail1_tmp >tail1
    sed 's/ H5Y DMPA/ H5S MY  /' tail1 >tail1_tmp
    sed 's/ H6X DMPA/ H6R MY  /' tail1_tmp >tail1
    sed 's/ H6Y DMPA/ H6S MY  /' tail1 >tail1_tmp
    sed 's/ H7X DMPA/ H7R MY  /' tail1_tmp >tail1
    sed 's/ H7Y DMPA/ H7S MY  /' tail1 >tail1_tmp
    sed 's/ H8X DMPA/ H8R MY  /' tail1_tmp >tail1
    sed 's/ H8Y DMPA/ H8S MY  /' tail1 >tail1_tmp
    sed 's/ H9X DMPA/ H9R MY  /' tail1_tmp >tail1
    sed 's/ H9Y DMPA/ H9S MY  /' tail1 >tail1_tmp
    sed 's/H10X DMPA/H10R MY  /' tail1_tmp >tail1
    sed 's/H10Y DMPA/H10S MY  /' tail1 >tail1_tmp
    sed 's/H11X DMPA/H11R MY  /' tail1_tmp >tail1
    sed 's/H11Y DMPA/H11S MY  /' tail1 >tail1_tmp
    sed 's/H12X DMPA/H12R MY  /' tail1_tmp >tail1
    sed 's/H12Y DMPA/H12S MY  /' tail1 >tail1_tmp
    sed 's/H13X DMPA/H13R MY  /' tail1_tmp >tail1
    sed 's/H13Y DMPA/H13S MY  /' tail1 >tail1_tmp
    sed 's/H14X DMPA/H14R MY  /' tail1_tmp >tail1
    sed 's/H14Y DMPA/H14S MY  /' tail1 >tail1_tmp
    sed 's/H14Z DMPA/H14T MY  /' tail1_tmp >tail1
    rm -f tail1_tmp

    #Tail 2
    sed 's/ C22 DMPA/ C12 MY  /' tail2 >tail2_tmp
    sed 's/ C23 DMPA/ C13 MY  /' tail2_tmp >tail2
    sed 's/ C24 DMPA/ C14 MY  /' tail2 >tail2_tmp
    sed 's/ C25 DMPA/ C15 MY  /' tail2_tmp >tail2
    sed 's/ C26 DMPA/ C16 MY  /' tail2 >tail2_tmp
    sed 's/ C27 DMPA/ C17 MY  /' tail2_tmp >tail2
    sed 's/ C28 DMPA/ C18 MY  /' tail2 >tail2_tmp
    sed 's/ C29 DMPA/ C19 MY  /' tail2_tmp >tail2
    sed 's/C210 DMPA/C110 MY  /' tail2 >tail2_tmp
    sed 's/C211 DMPA/C111 MY  /' tail2_tmp >tail2
    sed 's/C212 DMPA/C112 MY  /' tail2 >tail2_tmp
    sed 's/C213 DMPA/C113 MY  /' tail2_tmp >tail2
    sed 's/C214 DMPA/C114 MY  /' tail2 >tail2_tmp
    sed 's/ H2R DMPA/ H2R MY  /' tail2_tmp >tail2
    sed 's/ H2S DMPA/ H2S MY  /' tail2 >tail2_tmp
    sed 's/ H3R DMPA/ H3R MY  /' tail2_tmp >tail2
    sed 's/ H3S DMPA/ H3S MY  /' tail2 >tail2_tmp
    sed 's/ H4R DMPA/ H4R MY  /' tail2_tmp >tail2
    sed 's/ H4S DMPA/ H4S MY  /' tail2 >tail2_tmp
    sed 's/ H5R DMPA/ H5R MY  /' tail2_tmp >tail2
    sed 's/ H5S DMPA/ H5S MY  /' tail2 >tail2_tmp
    sed 's/ H6R DMPA/ H6R MY  /' tail2_tmp >tail2
    sed 's/ H6S DMPA/ H6S MY  /' tail2 >tail2_tmp
    sed 's/ H7R DMPA/ H7R MY  /' tail2_tmp >tail2
    sed 's/ H7S DMPA/ H7S MY  /' tail2 >tail2_tmp
    sed 's/ H8R DMPA/ H8R MY  /' tail2_tmp >tail2
    sed 's/ H8S DMPA/ H8S MY  /' tail2 >tail2_tmp
    sed 's/ H9R DMPA/ H9R MY  /' tail2_tmp >tail2
    sed 's/ H9S DMPA/ H9S MY  /' tail2 >tail2_tmp
    sed 's/H10R DMPA/H10R MY  /' tail2_tmp >tail2
    sed 's/H10S DMPA/H10S MY  /' tail2 >tail2_tmp
    sed 's/H11R DMPA/H11R MY  /' tail2_tmp >tail2
    sed 's/H11S DMPA/H11S MY  /' tail2 >tail2_tmp
    sed 's/H12R DMPA/H12R MY  /' tail2_tmp >tail2
    sed 's/H12S DMPA/H12S MY  /' tail2 >tail2_tmp
    sed 's/H13R DMPA/H13R MY  /' tail2_tmp >tail2
    sed 's/H13S DMPA/H13S MY  /' tail2 >tail2_tmp
    sed 's/H14R DMPA/H14R MY  /' tail2_tmp >tail2
    sed 's/H14S DMPA/H14S MY  /' tail2 >tail2_tmp
    sed 's/H14T DMPA/H14T MY  /' tail2_tmp >tail2
    rm -f tail2_tmp

    #Head
    sed 's/ P   DMPA/ P31 PH- /' head > head_tmp
    sed 's/ O13 DMPA/ O33 PH- /' head_tmp > head
    sed 's/ O14 DMPA/ O34 PH- /' head > head_tmp
    sed 's/ O12 DMPA/ O32 PH- /' head_tmp > head
    sed 's/ H12 DMPA/HO2A PH- /' head > head_tmp
    sed 's/ O11 DMPA/ O31 PH- /' head_tmp > head
    sed 's/ C1  DMPA/ C3  PH- /' head > head_tmp
    sed 's/ HA  DMPA/ HA  PH- /' head_tmp > head
    sed 's/ HB  DMPA/ HB  PH- /' head > head_tmp
    sed 's/ C2  DMPA/ C2  PH- /' head_tmp > head
    sed 's/ HS  DMPA/ HX  PH- /' head > head_tmp
    sed 's/ O21 DMPA/ O21 PH- /' head_tmp > head
    sed 's/ C21 DMPA/ C21 PH- /' head > head_tmp
    sed 's/ O22 DMPA/ O22 PH- /' head_tmp > head
    sed 's/ C3  DMPA/ C1  PH- /' head > head_tmp
    sed 's/ HX  DMPA/ HR  PH- /' head_tmp > head
    sed 's/ HY  DMPA/ HS  PH- /' head > head_tmp
    sed 's/ O31 DMPA/ O11 PH- /' head_tmp > head
    sed 's/ C31 DMPA/ C11 PH- /' head > head_tmp
    sed 's/ O32 DMPA/ O12 PH- /' head_tmp > head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat tail1 >> ../lipid_tmp_final
    cat head >> ../lipid_tmp_final
    cat tail2 >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

#CHL
if ( $chl_found == 1 ) then
  #CHL is 74 atoms per lipid. 
  set atoms_per_lipid = 74

  echo " Processing: CHL"

  #Check how many we have and see if it ads up correctly
  set line_count=`wc -l CHL_tmp | awk '{print$1}'`
  echo "  CHL Line count = $line_count"
  @ lipid_remainder = $line_count % $atoms_per_lipid
  if ( $lipid_remainder != 0 ) then
    echo " ERROR: Extracted CHL Lipids do NOT give a multiple of $atoms_per_lipid atoms."
    echo "                       Line Count = $line_count"
    echo "        Division by $atoms_per_lipid remainder = $lipid_remainder"
  endif
  @ lipid_count = $line_count / $atoms_per_lipid
  echo " CHL Lipid count = $lipid_count" 

  mkdir split_tmp
  cd split_tmp
  split -l $atoms_per_lipid ../CHL_tmp
  rm -f ../CHL_tmp 

  #Step 1 - rearrange atoms to be in the correct order for the Lipid11 force field.
  #SN1 tail goes first.
  set current_processing = 0
  echo -n " Processing CHL Lipid: "
  foreach i ( * )
    @ current_processing++
    echo -n " $current_processing "

    #Extract the CHL
    head -n74 $i > head

    #Modify residue and atom names
    sed 's/ C3  CHL1/ C3  CHL /' head >head_tmp
    sed 's/ O3  CHL1/ O1  CHL /' head_tmp >head
    sed s/" H3' CHL1"/" HO1 CHL "/ head >head_tmp
    sed 's/ H3  CHL1/ H31 CHL /' head_tmp >head
    sed 's/ C4  CHL1/ C4  CHL /' head >head_tmp
    sed 's/ H4A CHL1/ H41 CHL /' head_tmp >head
    sed 's/ H4B CHL1/ H42 CHL /' head >head_tmp
    sed 's/ C5  CHL1/ C5  CHL /' head_tmp >head
    sed 's/ C6  CHL1/ C6  CHL /' head >head_tmp
    sed 's/ H6  CHL1/ H61 CHL /' head_tmp >head
    sed 's/ C7  CHL1/ C7  CHL /' head >head_tmp
    sed 's/ H7A CHL1/ H71 CHL /' head_tmp >head
    sed 's/ H7B CHL1/ H72 CHL /' head >head_tmp
    sed 's/ C8  CHL1/ C8  CHL /' head_tmp >head
    sed 's/ H8  CHL1/ H81 CHL /' head >head_tmp
    sed 's/ C14 CHL1/ C14 CHL /' head_tmp >head
    sed 's/ H14 CHL1/H141 CHL /' head >head_tmp
    sed 's/ C15 CHL1/ C15 CHL /' head_tmp >head
    sed 's/H15A CHL1/H151 CHL /' head >head_tmp
    sed 's/H15B CHL1/H152 CHL /' head_tmp >head
    sed 's/ C16 CHL1/ C16 CHL /' head >head_tmp
    sed 's/H16A CHL1/H161 CHL /' head_tmp >head
    sed 's/H16B CHL1/H162 CHL /' head >head_tmp
    sed 's/ C17 CHL1/ C17 CHL /' head_tmp >head
    sed 's/ H17 CHL1/H171 CHL /' head >head_tmp
    sed 's/ C13 CHL1/ C13 CHL /' head_tmp >head
    sed 's/ C18 CHL1/ C18 CHL /' head >head_tmp
    sed 's/H18A CHL1/H181 CHL /' head_tmp >head
    sed 's/H18B CHL1/H182 CHL /' head >head_tmp
    sed 's/H18C CHL1/H183 CHL /' head_tmp >head
    sed 's/ C12 CHL1/ C12 CHL /' head >head_tmp
    sed 's/H12A CHL1/H121 CHL /' head_tmp >head
    sed 's/H12B CHL1/H122 CHL /' head >head_tmp
    sed 's/ C11 CHL1/ C11 CHL /' head_tmp >head
    sed 's/H11A CHL1/H111 CHL /' head >head_tmp
    sed 's/H11B CHL1/H112 CHL /' head_tmp >head
    sed 's/ C9  CHL1/ C9  CHL /' head >head_tmp
    sed 's/ H9  CHL1/ H91 CHL /' head_tmp >head
    sed 's/ C10 CHL1/ C10 CHL /' head >head_tmp
    sed 's/ C19 CHL1/ C19 CHL /' head_tmp >head
    sed 's/H19A CHL1/H191 CHL /' head >head_tmp
    sed 's/H19B CHL1/H192 CHL /' head_tmp >head
    sed 's/H19C CHL1/H193 CHL /' head >head_tmp
    sed 's/ C1  CHL1/ C1  CHL /' head_tmp >head
    sed 's/ H1A CHL1/ H11 CHL /' head >head_tmp
    sed 's/ H1B CHL1/ H12 CHL /' head_tmp >head
    sed 's/ C2  CHL1/ C2  CHL /' head >head_tmp
    sed 's/ H2A CHL1/ H21 CHL /' head_tmp >head
    sed 's/ H2B CHL1/ H22 CHL /' head >head_tmp
    sed 's/ C20 CHL1/ C20 CHL /' head_tmp >head
    sed 's/ H20 CHL1/H201 CHL /' head >head_tmp
    sed 's/ C21 CHL1/ C21 CHL /' head_tmp >head
    sed 's/H21A CHL1/H211 CHL /' head >head_tmp
    sed 's/H21B CHL1/H212 CHL /' head_tmp >head
    sed 's/H21C CHL1/H213 CHL /' head >head_tmp
    sed 's/ C22 CHL1/ C22 CHL /' head_tmp >head
    sed 's/H22A CHL1/H221 CHL /' head >head_tmp
    sed 's/H22B CHL1/H222 CHL /' head_tmp >head
    sed 's/ C23 CHL1/ C23 CHL /' head >head_tmp
    sed 's/H23A CHL1/H231 CHL /' head_tmp >head
    sed 's/H23B CHL1/H232 CHL /' head >head_tmp
    sed 's/ C24 CHL1/ C24 CHL /' head_tmp >head
    sed 's/H24A CHL1/H241 CHL /' head >head_tmp
    sed 's/H24B CHL1/H242 CHL /' head_tmp >head
    sed 's/ C25 CHL1/ C25 CHL /' head >head_tmp
    sed 's/ H25 CHL1/H251 CHL /' head_tmp >head
    sed 's/ C26 CHL1/ C26 CHL /' head >head_tmp
    sed 's/H26A CHL1/H261 CHL /' head_tmp >head
    sed 's/H26B CHL1/H262 CHL /' head >head_tmp
    sed 's/H26C CHL1/H263 CHL /' head_tmp >head
    sed 's/ C27 CHL1/ C27 CHL /' head >head_tmp
    sed 's/H27A CHL1/H271 CHL /' head_tmp >head
    sed 's/H27B CHL1/H272 CHL /' head >head_tmp
    sed 's/H27C CHL1/H273 CHL /' head_tmp >head
    rm -f head_tmp

    #Append the modified lipid to the lipid_tmp_final file.
    cat head >> ../lipid_tmp_final
    echo "TER " >> ../lipid_tmp_final 

  end  

  cd ../
  rm -rf split_tmp
  echo " "
endif

echo " *** STAGE 6 : COMBINE PROTEIN, LIPID, ION AND WATER FILES INTO OUTPUT PDB ***"
echo " "

if ( $proteins_found == 1 ) then
  mv protein_tmp $output_filename
  echo "TER " >>$output_filename
else
  touch $output_filename
endif

cat lipid_tmp_final >> $output_filename
rm -f lipid_tmp_final

if ( $ions_found == 1 ) then
  cat ions_tmp >> $output_filename
  rm -f ions_tmp
  echo "TER " >> $output_filename
endif

if ( $water_found == 1 ) then
  cat water_tmp >> $output_filename
  rm -f water_tmp
  echo "TER " >> $output_filename
endif

echo "END " >> $output_filename

echo " ******* CONVERSION COMPLETE ******* "
echo " "

exit(0)

