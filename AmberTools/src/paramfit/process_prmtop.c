/*****************************************************
 * AMBER Bond Angle and Dihedral Parameter Optimiser *
 *                                                   *
 *           Written by: Robin Betz  (2011)          *
 *                       Ross Walker (2004)          *
 *                   UC San Diego                    *
 *           San Diego Supercomputer Center          *
 *            La Jolla, California, 92092            *
 *                       USA                         *
 *****************************************************/


/*process_prmtop.c*/

/*This routine is responsible for processing the prmtop file into seperate parameters for
  each bond angle and dihedral type.

  This is necessary because the prmtop file and thus prmtop storage arrays do not keep seperate
  parameters for different bond / angle or dihedral setups that share the same values for these parameters

  E.g O-C-N-H and O-C-N-CH3 share the same parameters and so this is stored in the same spot in memory

  We need to split this in order to optimise them individually
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "function_def.h"
#include "constants.h"

int process_prmtop(global_options_struct *global_options, parm_struct *parm_data)
{
  int i,j;
  int atom1, atom2, atom3, atom4;
  double rk, req, tk, teq, pk, pn, phase;
  char type1[NAME_SIZE];
  char type2[NAME_SIZE];
  char type3[NAME_SIZE];
  char type4[NAME_SIZE];    
  short int is_unique;
  short int is_improper;
  char bond_make_up1[(NAME_SIZE*2)-1];
  char bond_make_up2[(NAME_SIZE*2)-1];
  char angle_make_up1[(NAME_SIZE*3)-2];
  char angle_make_up2[(NAME_SIZE*3)-2];
  char dihedral_make_up1[(NAME_SIZE*4)-3];
  char dihedral_make_up2[(NAME_SIZE*4)-3];    
  int check1,check2;
  char response = '\0';

  /* Ask the user some general categories so as to overwhelm them with less output */
  short int prompt_bonds=YES, prompt_angles=YES, prompt_dihedrals=YES;
  short int b_req=NO, b_kr=NO;
  short int a_kt=NO, a_theq=NO;
  short int d_kp=NO, d_np=NO, d_phase=NO;
  
  /* Make sure a file name has been specified to save in if necessary */
  if (global_options->RUNTYPE==SET_PARAMS && global_options->PARAMETER_FILE_NAME==NULL)
  {
    printf("\n\n!  ERROR: No location specified to save parameter file. Please\n");
    printf("          read the documentation on how to use this runtype and   \n");
    printf("          specify a PARAMETER_FILE_NAME.\n");
    return FAILURE;
  }
  /* If a filename has been specified and it's not save, warn the user of ambiguity */
  if ( global_options->PARAMETERS_TO_FIT==DEFAULT && global_options->PARAMETER_FILE_NAME!=NULL && 
       global_options->RUNTYPE != SET_PARAMS && global_options->VERBOSITY >=MEDIUM )
  {
      printf("! Ambiguous parameters to fit options- parameter file and DEFAULT parameters specified.\n");
      return FAILURE;
  }
  
  /* If a filename has been specified and it's K only, warn the user */
  if (global_options->PARAMETERS_TO_FIT==K_ONLY && global_options->PARAMETER_FILE_NAME!=NULL)
  {
    printf("! Ambiguous parameters to fit options- parameter file and K_ONLY parameters specified.\n");
    return FAILURE;
  }
  
  if (global_options->RUNTYPE==SET_PARAMS)
  {
    global_options->PARAMETERS_TO_FIT = SAVE;
    
    printf("  Would you like to optimise any bond parameters? (y/n):\n");
    while( response!='y' && response!='n' )
      scanf("%c",&response);
    prompt_bonds = (response == 'y') ? YES : NO;
    response = '\0';
    fflush(stdout);

    if (prompt_bonds == YES)
    {
      printf("    - Any bond REQ? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      b_req = (response == 'y') ? YES : NO;
      response = '\0';
      fflush(stdout);
      
      printf("    - Any bond KR? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      b_kr = (response == 'y') ? YES : NO;
      response = '\0';
      fflush(stdout);
    }
    
    printf("  Would you like to optimise any angle parameters? (y/n): \n");
    fflush(stdout);
    while( response!='y' && response!='n' )
      scanf("%c",&response);
    prompt_angles = (response == 'y') ? YES : NO;
    response = '\0';
    if (prompt_angles == YES)
    {
      printf("    - Any angle KT? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      a_kt = (response == 'y') ? YES : NO;
      response = '\0';
      fflush(stdout);
      
      printf("    - Any angle THEQ? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      a_theq = (response == 'y') ? YES : NO;
      response = '\0';
      fflush(stdout);
    }
   
   printf("  Would you like to optimise any dihedral parameters? (y/n): \n");
    while( response!='y' && response!='n' )
      scanf("%c",&response);
    prompt_dihedrals = (response == 'y') ? YES : NO;    
    fflush(stdout);
    
    response = '\0';
    if (prompt_dihedrals == YES)
    {
      printf("    - Any dihedral KP? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      d_kp = (response == 'y') ? YES : NO;
      fflush(stdout);
      
      response = '\0';
      printf("    - Any dihedral NP? (y/n): \n" );
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      d_np = (response == 'y') ? YES : NO;
      fflush(stdout);
      
      response = '\0';
      printf("    - Any dihedral PHASE? (y/n): \n");
      while(response != 'y' && response != 'n')
        scanf("%c", &response);
      d_phase = (response == 'y') ? YES : NO;
      response = '\0';
      fflush(stdout);
    }
  }

  /*We should allocate at least 1 block to the bond, angle and dihedral data structs to avoid problems
    with reallocation*/
  if (global_options->VERBOSITY>=HIGH)
     printf("   Allocating %d bytes for parm_data->*bond_data\n",(int) sizeof(bond_data_struct));
  parm_data->bond_data = (bond_data_struct *)malloc(sizeof(bond_data_struct));
  if (parm_data->bond_data == NULL)
  {
    malloc_failure_char("process_prmtop", "parm_data->bond_data", sizeof(bond_data_struct));
    return ALLOC_FAIL;
  }
  parm_data->mem_allocated+=sizeof(bond_data_struct);
  if (global_options->VERBOSITY>=HIGH)
     printf("   Allocating %d bytes for parm_data->*angle_data\n",(int) sizeof(angle_data_struct));
  parm_data->angle_data = (angle_data_struct *)malloc(sizeof(angle_data_struct));
  if (parm_data->angle_data == NULL)
  {
    malloc_failure_char("process_prmtop", "parm_data->angle_data", sizeof(angle_data_struct));
    return ALLOC_FAIL;
  }
  parm_data->mem_allocated+=sizeof(angle_data_struct);
  if (global_options->VERBOSITY>=HIGH)
     printf("   Allocating %d bytes for parm_data->*dihedral_data\n",(int) sizeof(dihedral_data_struct));
  parm_data->dihedral_data = (dihedral_data_struct *)malloc(sizeof(dihedral_data_struct));
  if (parm_data->dihedral_data == NULL)
  {
    malloc_failure_char("process_prmtop", "parm_data->dihedral_data", sizeof(dihedral_data_struct));
    return ALLOC_FAIL;
  }
  parm_data->mem_allocated+=sizeof(dihedral_data_struct);
  
  if (global_options->VERBOSITY>=HIGH)
    printf("   Prmtop   (unique): Processing prmtop to find unique bonds, angles & dihedrals.\n");
  
  parm_data->unique_bonds_found=0;
  parm_data->unique_angles_found=0;
  parm_data->unique_dihedrals_found=0;    
  /*First of all we will loop over all the bonds and split them into their unique types*/
  /*We will go through each bond in turn and assign it to a location in bond_data. The problem
    we have is that in theory two bond terms could exist that have the same types but different
    parameters, this is not an accute problem for bonds and angles but is for dihedrals*/
  
  for (i=0;i<parm_data->NBONH;++i)
  {
     /*Loop over all bonds involving hydrogen*/
     /*Get the atoms involved in the bond we are looking at*/
     atom1=unObfuscateAtom(parm_data->pbondH[i].ib);
     atom2=unObfuscateAtom(parm_data->pbondH[i].jb);
     /*Get the parameters involved*/
     rk=parm_data->rk[parm_data->pbondH[i].icb-1];
     req=parm_data->req[parm_data->pbondH[i].icb-1];     
     /*Get the atom types*/
     strcpy(type1,parm_data->atom[atom1-1].isymbl);
     strcpy(type2,parm_data->atom[atom2-1].isymbl);
     
     /*Ok, now we have the data for this bond*/
     /*Now we need to determine if it is unique*/
     /*loop over all unique bonds to date and see if we get a match*/
     is_unique=YES;
     for (j=0;j<parm_data->unique_bonds_found;++j)
     {
       /*A bond is unique if both atoms types don't match, and both parameters don't match*/
       /*Note, we can do the parameters in both directions*/
       /*First try the bond one way, then the reverse since directionality is irrelevant*/
       strcpy(bond_make_up1,type1);
       strcat(bond_make_up1,type2);
       strcpy(bond_make_up2,parm_data->bond_data[j].atom_type1);
       strcat(bond_make_up2,parm_data->bond_data[j].atom_type2);
       check1=!strcmp(bond_make_up1,bond_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
       strcpy(bond_make_up1,type2);
       strcat(bond_make_up1,type1);
       strcpy(bond_make_up2,parm_data->bond_data[j].atom_type1);
       strcat(bond_make_up2,parm_data->bond_data[j].atom_type2);
       check2=!strcmp(bond_make_up1,bond_make_up2);             /*strcmp returns zero on a match - so check2 will be 1 if reverse matched*/
       /*if the atom type is not unique and matches here see if the parameters also match
         if they do then we have found the matching bond type for this bond so add the data and break
         else we keep looping*/
       if (check1 || check2)
       {
         if (rk == parm_data->bond_data[j].rk && req == parm_data->bond_data[j].req)
         {
           /*bond type is not unique*/
           is_unique=NO;
           /*fill in the data*/
           ++(parm_data->bond_data[j].number);
           /*WE NEED TO CHECK HERE THAT WE HAVEN'T OVERFLOWED OUR ARRAY*/
           if (parm_data->bond_data[j].number>=MAX_BONDS_PER_TYPE)
           {
             printf("*** ERROR IN PROCESS_PRMTOP - MAX_BONDS_PER_TYPE OF %d\n",MAX_BONDS_PER_TYPE);
             printf("*** EXCEEDED FOR parm_data->bond_data[j], j = %d\n",j);
             printf("*** EDIT MAX_BONDS_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
             return DATA_OVERFLOW;
           }
           parm_data->bond_data[j].atom1[(parm_data->bond_data[j].number)-1]=atom1;
           parm_data->bond_data[j].atom2[(parm_data->bond_data[j].number)-1]=atom2;
           
           break; /*Quit the loop, no point doing the rest if we have matched at this point*/
         }
       }
     } /*essentially if we finish this loop without a hit it must be unique*/
     if (is_unique==YES)
     {
       /*Bond is unique, have to extend number of bond_data structures by 1*/
       ++(parm_data->unique_bonds_found);
       
       /*reallocate memory for the bond_data*/
       if (parm_data->unique_bonds_found>1)
       {
         if (global_options->VERBOSITY>=HIGH)
             printf("   Reallocating %d bytes for parm_data->*bond_data\n",(int) sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found);
          parm_data->bond_data = (struct _bond_data_struct *) realloc(parm_data->bond_data,sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found);
          if (parm_data->bond_data == NULL)
          {
            malloc_failure_char("process_prmtop", "parm_data->bond_data", sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found);
            return ALLOC_FAIL;
          }

         parm_data->mem_allocated+=sizeof(struct _bond_data_struct);
       }
       /*now put the data in the bond_structure*/
       parm_data->bond_data[parm_data->unique_bonds_found-1].number=1;
       strcpy(parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,type1);
       strcpy(parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2,type2);
       parm_data->bond_data[parm_data->unique_bonds_found-1].rk=rk;
       parm_data->bond_data[parm_data->unique_bonds_found-1].req=req;
       parm_data->bond_data[parm_data->unique_bonds_found-1].atom1[0]=atom1;
       parm_data->bond_data[parm_data->unique_bonds_found-1].atom2[0]=atom2;
       
       /*Need to check if we are to fit this parameter or not*/
       parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=NO;
       parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=NO;
       if ( global_options->RUNTYPE==SET_PARAMS && prompt_bonds )
        {
          /*prompt the user whether to fit or not*/
          if (b_kr==YES)
          {
            response='\0';
            printf("Fit Parameter: BOND (%s-%s) KR? (y/n): ",
            parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,
            parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=YES;
          }

          if (b_req==YES)
          {
            response='\0';
            printf("Fit Parameter: BOND (%s-%s) REQ? (y/n): ",
                  parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,
                  parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=YES;
          }
        }
        else if (global_options->PARAMETERS_TO_FIT==DEFAULT)
        {
          parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=YES;  /*DEFAULT IS TO FIT EVERYTHING, this will be adjusted as necessary by other routines*/
          parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=YES;
        }
        else if (global_options->PARAMETERS_TO_FIT==K_ONLY)
        {
          parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=NO;
          parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=NO;
        }
     }
  } /*End of bonds with Hydrogen*/

  /*Now do bonds without hydrogen*/
  for (i=0;i<parm_data->NBONA;++i)
  {
     /*Get the atoms involved in the bond we are looking at*/
     atom1=unObfuscateAtom(parm_data->pbond[i].ib);
     atom2=unObfuscateAtom(parm_data->pbond[i].jb);
     /*Get the parameters involved*/
     rk=parm_data->rk[parm_data->pbond[i].icb-1];
     req=parm_data->req[parm_data->pbond[i].icb-1];
     /*Get the atom types*/
     strcpy(type1,parm_data->atom[atom1-1].isymbl);
     strcpy(type2,parm_data->atom[atom2-1].isymbl);
     
     /*Ok, now we have the data for this bond*/
     /*Now we need to determine if it is unique*/
     /*loop over all unique bonds to date and see if we get a match*/
     is_unique=YES;
     for (j=0;j<parm_data->unique_bonds_found;++j)
     {
       /*A bond is unique if both atoms types don't match, and both parameters don't match*/
       /*Note, we can do the parameters in both directions*/
       /*First try the bond one way, then the reverse since directionality is irrelevant*/
       strcpy(bond_make_up1,type1);
       strcat(bond_make_up1,type2);
       strcpy(bond_make_up2,parm_data->bond_data[j].atom_type1);
       strcat(bond_make_up2,parm_data->bond_data[j].atom_type2);
       check1=!strcmp(bond_make_up1,bond_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
       strcpy(bond_make_up1,type2);
       strcat(bond_make_up1,type1);
       strcpy(bond_make_up2,parm_data->bond_data[j].atom_type1);
       strcat(bond_make_up2,parm_data->bond_data[j].atom_type2);
       check2=!strcmp(bond_make_up1,bond_make_up2);             /*strcmp returns zero on a match - so check2 will be 1 if reverse matched*/
       /*if the atom type is not unique and matches here see if the parameters also match
         if they do then we have found the matching bond type for this bond so add the data and break
         else we keep looping*/
       if (check1 || check2)
       {
         if (rk == parm_data->bond_data[j].rk && req == parm_data->bond_data[j].req)
         {
           /*bond type is not unique*/
           is_unique=NO;
           /*fill in the data*/
           ++(parm_data->bond_data[j].number);
           /*WE NEED TO CHECK HERE THAT WE HAVEN'T OVERFLOWED OUR ARRAY*/
           if (parm_data->bond_data[j].number>=MAX_BONDS_PER_TYPE)
           {
             printf("*** ERROR IN PROCESS_PRMTOP - MAX_BONDS_PER_TYPE OF %d\n",MAX_BONDS_PER_TYPE);
             printf("*** EXCEEDED FOR parm_data->bond_data[j], j = %d\n",j);
             printf("*** EDIT MAX_BONDS_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
             return DATA_OVERFLOW;
           }
           parm_data->bond_data[j].atom1[(parm_data->bond_data[j].number)-1]=atom1;
           parm_data->bond_data[j].atom2[(parm_data->bond_data[j].number)-1]=atom2;
           
           break; /*Quit the loop, no point doing the rest if we have matched at this point*/
         }
       }
     } /*essentially if we finish this loop without a hit it must be unique*/
     if (is_unique==YES)
     {
       /*Bond is unique, have to extend number of bond_data structures by 1*/
       ++(parm_data->unique_bonds_found);
       /*reallocate memory for the bond_data*/
       if (parm_data->unique_bonds_found>1)
       {
         if (global_options->VERBOSITY>=HIGH)
             printf("   Reallocating %d bytes for parm_data->*bond_data\n",(int) (sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found));
          parm_data->bond_data = (struct _bond_data_struct *) realloc(parm_data->bond_data,sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found);
          if (parm_data->bond_data == NULL)
          {
            malloc_failure_char("process_prmtop", "parm_data->bond_data", sizeof(struct _bond_data_struct)*parm_data->unique_bonds_found);
            return ALLOC_FAIL;
          }
          parm_data->mem_allocated+=sizeof(struct _bond_data_struct);
       }
       /*now put the data in the bond_structure*/
       parm_data->bond_data[parm_data->unique_bonds_found-1].number=1;
       strcpy(parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,type1);
       strcpy(parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2,type2);
       parm_data->bond_data[parm_data->unique_bonds_found-1].rk=rk;
       parm_data->bond_data[parm_data->unique_bonds_found-1].req=req;
       parm_data->bond_data[parm_data->unique_bonds_found-1].atom1[0]=atom1;
       parm_data->bond_data[parm_data->unique_bonds_found-1].atom2[0]=atom2;
       
       /*Need to check if we are to fit this parameter or not*/
       parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=NO;
       parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=NO;
       if (global_options->RUNTYPE==SET_PARAMS)
        {
          /*prompt the user whether to fit or not*/
          if (b_kr == YES)
          {
            response='\0';
            printf("Fit Parameter: BOND (%s-%s) KR? (y/n): ",
            parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,
            parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=YES;
          }
          
          if (b_req==YES)
          {
            response='\0';
            printf("Fit Parameter: BOND (%s-%s) REQ? (y/n): ",
                  parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type1,
                  parm_data->bond_data[parm_data->unique_bonds_found-1].atom_type2);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=YES;
          }
        }
        else if (global_options->PARAMETERS_TO_FIT==DEFAULT)
        {
          parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=YES;  /*DEFAULT IS TO FIT EVERYTHING, this will be adjusted as necessary by other routines*/
          parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=YES;
        }
        else if (global_options->PARAMETERS_TO_FIT==K_ONLY)
        {
          parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_KR=NO;
          parm_data->bond_data[parm_data->unique_bonds_found-1].DO_FIT_REQ=NO;
        }
     }
  } /*End of bonds withOUT Hydrogen*/
  
  /*Now we do angles with hydrogen*/
  for (i=0;i<parm_data->NTHETH;++i)
  {
     /*Loop over all angles involving hydrogen*/
     /*Get the atoms involved in the angle we are looking at*/
     atom1=unObfuscateAtom(parm_data->pangleH[i].it);
     atom2=unObfuscateAtom(parm_data->pangleH[i].jt);
     atom3=unObfuscateAtom(parm_data->pangleH[i].kt);
     /*Get the parameters involved*/
     tk=parm_data->tk[parm_data->pangleH[i].ict-1];
     teq=parm_data->teq[parm_data->pangleH[i].ict-1];
     /*Get the atom types*/
     strcpy(type1,parm_data->atom[atom1-1].isymbl);
     strcpy(type2,parm_data->atom[atom2-1].isymbl);
     strcpy(type3,parm_data->atom[atom3-1].isymbl);
     

     /*Ok, now we have the data for this angle*/
     /*Now we need to determine if it is unique*/
     /*loop over all unique angles to date and see if we get a match*/
     is_unique=YES;
     for (j=0;j<parm_data->unique_angles_found;++j)
     {
       /*An angle is unique if all atoms types don't match, and both parameters don't match*/
       /*Note, we can do the atom types in both directions*/
       /*First try the angle one way, then the reverse since directionality is irrelevant*/
       strcpy(angle_make_up1,type1);
       strcat(angle_make_up1,type2);
       strcat(angle_make_up1,type3);       
       strcpy(angle_make_up2,parm_data->angle_data[j].atom_type1);
       strcat(angle_make_up2,parm_data->angle_data[j].atom_type2);
       strcat(angle_make_up2,parm_data->angle_data[j].atom_type3);       
       check1=!strcmp(angle_make_up1,angle_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
       strcpy(angle_make_up1,type3);
       strcat(angle_make_up1,type2);
       strcat(angle_make_up1,type1);
       strcpy(angle_make_up2,parm_data->angle_data[j].atom_type1);
       strcat(angle_make_up2,parm_data->angle_data[j].atom_type2);
       strcat(angle_make_up2,parm_data->angle_data[j].atom_type3);
       check2=!strcmp(angle_make_up1,angle_make_up2);             /*strcmp returns zero on a match - so check2 will be 1 if reverse matched*/
       /*if the atom type is not unique and matches here see if the parameters also match
         if they do then we have found the matching angle type for this angle so add the data and break
         else we keep looping*/
       if (check1 || check2)
       {
         if (tk == parm_data->angle_data[j].tk && teq == parm_data->angle_data[j].teq)
         {
           /*angle type is not unique*/
           is_unique=NO;
           /*fill in the data*/
           ++(parm_data->angle_data[j].number);
           /*WE NEED TO CHECK HERE THAT WE HAVEN'T OVERFLOWED OUR ARRAY*/
           if (parm_data->angle_data[j].number>=MAX_ANGLES_PER_TYPE)
           {
             printf("*** ERROR IN PROCESS_PRMTOP - MAX_ANGLES_PER_TYPE OF %d\n",MAX_ANGLES_PER_TYPE);
             printf("*** EXCEEDED FOR parm_data->angle_data[j], j = %d\n",j);
             printf("*** EDIT MAX_ANGLES_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
             return DATA_OVERFLOW;
           }
           parm_data->angle_data[j].atom1[(parm_data->angle_data[j].number)-1]=atom1;
           parm_data->angle_data[j].atom2[(parm_data->angle_data[j].number)-1]=atom2;
           parm_data->angle_data[j].atom3[(parm_data->angle_data[j].number)-1]=atom3;
           
           break; /*Quit the loop, no point doing the rest if we have matched at this point*/
         }
       }
     } /*essentially if we finish this loop without a hit it must be unique*/
     if (is_unique==YES)
     {
       /*Angle is unique, have to extend number of angle_data structures by 1*/
       ++(parm_data->unique_angles_found);
       /*reallocate memory for the angle_data*/
       if (parm_data->unique_angles_found>1)
       {
         if (global_options->VERBOSITY>=HIGH)
             printf("   Reallocating %d bytes for parm_data->*angle_data\n",(int) (sizeof(struct _angle_data_struct)*parm_data->unique_angles_found));
          parm_data->angle_data = (struct _angle_data_struct *) realloc(parm_data->angle_data,sizeof(struct _angle_data_struct)*parm_data->unique_angles_found);
          if (parm_data->angle_data == NULL)
          {
            malloc_failure_char("process_prmtop", "parm_data->angle_data", sizeof(struct _angle_data_struct)*parm_data->unique_angles_found);
            return ALLOC_FAIL;
          }
          parm_data->mem_allocated+=sizeof(struct _angle_data_struct);
       }
       /*now put the data in the angle_structure*/
       parm_data->angle_data[parm_data->unique_angles_found-1].number=1;
       strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,type1);
       strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,type2);
       strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3,type3);
       parm_data->angle_data[parm_data->unique_angles_found-1].tk=tk;
       parm_data->angle_data[parm_data->unique_angles_found-1].teq=teq;
       parm_data->angle_data[parm_data->unique_angles_found-1].atom1[0]=atom1;
       parm_data->angle_data[parm_data->unique_angles_found-1].atom2[0]=atom2;
       parm_data->angle_data[parm_data->unique_angles_found-1].atom3[0]=atom3;
       
       /*Need to check if we are to fit this parameter or not*/
       parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=NO;
       parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=NO;
       if (global_options->RUNTYPE==SET_PARAMS)
        {
          /*prompt the user whether to fit or not*/
          if (a_kt == YES)
          {
            response='\0';
            printf("Fit Parameter: ANGLE (%s-%s-%s) KT? (y/n): ",
            parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,
            parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,
            parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=YES;
          }
          
        if (a_theq == YES)
        {
          response='\0';
          printf("Fit Parameter: ANGLE (%s-%s-%s) THEQ? (y/n): ",
                 parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,
                 parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,
                 parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3);
          fflush(stdout);
          while( response!='y' && response!='n' )
          {
            scanf("%c",&response);
          }
          if (response=='y')
            parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=YES;
        }
        }
        else if (global_options->PARAMETERS_TO_FIT==DEFAULT)
        {
          parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=YES; 
          parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=YES;
        }
        else if (global_options->PARAMETERS_TO_FIT==K_ONLY)
        {
          parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=NO; 
          parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=NO;
        }
     }
  } /*End of angles with Hydrogen*/

  /*Now we do angles withOUT hydrogen*/
  for (i=0;i<parm_data->NTHETA;++i)
  {
     /*Loop over all angles not involving hydrogen*/
     /*Get the atoms involved in the angle we are looking at*/
     atom1=unObfuscateAtom(parm_data->pangle[i].it);
     atom2=unObfuscateAtom(parm_data->pangle[i].jt);
     atom3=unObfuscateAtom(parm_data->pangle[i].kt);
     /*Get the parameters involved*/
     tk=parm_data->tk[parm_data->pangle[i].ict-1];
     teq=parm_data->teq[parm_data->pangle[i].ict-1];
     /*Get the atom types*/
     strcpy(type1,parm_data->atom[atom1-1].isymbl);
     strcpy(type2,parm_data->atom[atom2-1].isymbl);
     strcpy(type3,parm_data->atom[atom3-1].isymbl);


     /*Ok, now we have the data for this angle*/
     /*Now we need to determine if it is unique*/
     /*loop over all unique angles to date and see if we get a match*/
     is_unique=YES;
     for (j=0;j<parm_data->unique_angles_found;++j)
     {
       /*An angle is unique if all atoms types don't match, and both parameters don't match*/
       /*Note, we can do the atom types in both directions*/
       /*First try the angle one way, then the reverse since directionality is irrelevant*/
       strcpy(angle_make_up1,type1);
       strcat(angle_make_up1,type2);
       strcat(angle_make_up1,type3);
       strcpy(angle_make_up2,parm_data->angle_data[j].atom_type1);
       strcat(angle_make_up2,parm_data->angle_data[j].atom_type2);
       strcat(angle_make_up2,parm_data->angle_data[j].atom_type3);
       check1=!strcmp(angle_make_up1,angle_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
       strcpy(angle_make_up1,type3);
       strcat(angle_make_up1,type2);
       strcat(angle_make_up1,type1);
       strcpy(angle_make_up2,parm_data->angle_data[j].atom_type1);
       strcat(angle_make_up2,parm_data->angle_data[j].atom_type2);
       strcat(angle_make_up2,parm_data->angle_data[j].atom_type3);
       check2=!strcmp(angle_make_up1,angle_make_up2);             /*strcmp returns zero on a match - so check2 will be 1 if reverse matched*/
       /*if the atom type is not unique and matches here see if the parameters also match
         if they do then we have found the matching angle type for this angle so add the data and break
         else we keep looping*/
       if (check1 || check2)
       {
         if (tk == parm_data->angle_data[j].tk && teq == parm_data->angle_data[j].teq)
         {
           /*angle type is not unique*/
           is_unique=NO;
           /*fill in the data*/
           ++(parm_data->angle_data[j].number);
           /*WE NEED TO CHECK HERE THAT WE HAVEN'T OVERFLOWED OUR ARRAY*/
           if (parm_data->angle_data[j].number>=MAX_ANGLES_PER_TYPE)
           {
             printf("*** ERROR IN PROCESS_PRMTOP - MAX_ANGLES_PER_TYPE OF %d\n",MAX_ANGLES_PER_TYPE);
             printf("*** EXCEEDED FOR parm_data->angle_data[j], j = %d\n",j);
             printf("*** EDIT MAX_ANGLES_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
             return DATA_OVERFLOW;
           }
           parm_data->angle_data[j].atom1[(parm_data->angle_data[j].number)-1]=atom1;
           parm_data->angle_data[j].atom2[(parm_data->angle_data[j].number)-1]=atom2;
           parm_data->angle_data[j].atom3[(parm_data->angle_data[j].number)-1]=atom3;
           
           break; /*Quit the loop, no point doing the rest if we have matched at this point*/
         }
       }
     } /*essentially if we finish this loop without a hit it must be unique*/
     if (is_unique==YES)
     {
       /*Angle is unique, have to extend number of angle_data structures by 1*/
       ++(parm_data->unique_angles_found);
       /*reallocate memory for the angle_data*/
       if (parm_data->unique_angles_found>1)
       {
         if (global_options->VERBOSITY>=HIGH)
             printf("   Reallocating %d bytes for parm_data->*angle_data\n",(int) (sizeof(struct _angle_data_struct)*parm_data->unique_angles_found));
          parm_data->angle_data = (struct _angle_data_struct *) realloc(parm_data->angle_data,sizeof(struct _angle_data_struct)*parm_data->unique_angles_found);
          if (parm_data->angle_data == NULL)
          {
            malloc_failure_char("process_prmtop", "parm_data->angle_data", sizeof(struct _angle_data_struct)*parm_data->unique_angles_found);
            return ALLOC_FAIL;
          }
          parm_data->mem_allocated+=sizeof(struct _angle_data_struct);
       }
       /*now put the data in the angle_structure*/
       parm_data->angle_data[parm_data->unique_angles_found-1].number=1;
       strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,type1);
       strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,type2);
       strcpy(parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3,type3);
       parm_data->angle_data[parm_data->unique_angles_found-1].tk=tk;
       parm_data->angle_data[parm_data->unique_angles_found-1].teq=teq;
       parm_data->angle_data[parm_data->unique_angles_found-1].atom1[0]=atom1;
       parm_data->angle_data[parm_data->unique_angles_found-1].atom2[0]=atom2;
       parm_data->angle_data[parm_data->unique_angles_found-1].atom3[0]=atom3;
       
       /*Need to check if we are to fit this parameter or not*/
       parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=NO;
       parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=NO;
       if (global_options->RUNTYPE==SET_PARAMS)
        {
          /*prompt the user whether to fit or not*/
          if (a_kt==YES)
          {
            response='\0';
            printf("Fit Parameter: ANGLE (%s-%s-%s) KT? (y/n): ",
            parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,
            parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,
            parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=YES;
          }
          
          if (a_theq==YES)
          {
            response='\0';
            printf("Fit Parameter: ANGLE (%s-%s-%s) THEQ? (y/n): ",
                   parm_data->angle_data[parm_data->unique_angles_found-1].atom_type1,
                   parm_data->angle_data[parm_data->unique_angles_found-1].atom_type2,
                   parm_data->angle_data[parm_data->unique_angles_found-1].atom_type3);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=YES;
          }
        }
        else if (global_options->PARAMETERS_TO_FIT == DEFAULT)
        {
          parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=YES;
          parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=YES;
        }
        else if (global_options->PARAMETERS_TO_FIT == K_ONLY)
        {
          parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=NO;
          parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=NO;
        }
     }
  } /*End of angles withOUT Hydrogen*/
  
  /*Now we do dihedrals with hydrogen*/
  for (i=0;i<parm_data->NPHIH;++i)
  {
     /*Loop over all dihedrals involving hydrogen*/
     // Check if it's an improper dihedral
     if (parm_data->pdihedralH[i].lp < 0)
       is_improper=YES;
     else
       is_improper=NO;
     
     /*Get the atoms involved in the dihedral we are looking at*/
     atom1=unObfuscateAtom(parm_data->pdihedralH[i].ip);
     atom2=unObfuscateAtom(parm_data->pdihedralH[i].jp);
     atom3=unObfuscateAtom(parm_data->pdihedralH[i].kp);
     atom4=unObfuscateAtom(parm_data->pdihedralH[i].lp);     
     /*Get the parameters involved*/
     pk=parm_data->pk[parm_data->pdihedralH[i].icp-1];
     pn=parm_data->pn[parm_data->pdihedralH[i].icp-1];
     phase=parm_data->phase[parm_data->pdihedralH[i].icp-1];
     
     /*Get the atom types*/
     strcpy(type1,parm_data->atom[atom1-1].isymbl);
     strcpy(type2,parm_data->atom[atom2-1].isymbl);
     strcpy(type3,parm_data->atom[atom3-1].isymbl);
     strcpy(type4,parm_data->atom[atom4-1].isymbl);

     /*Ok, now we have the data for this dihedral*/
     /*Now we need to determine if it is unique*/
     /*loop over all unique dihedrals to date and see if we get a match*/
     is_unique=YES;
     for (j=0;j<parm_data->unique_dihedrals_found;++j)
     {
       /*A dihedral is unique if all atom types don't match, and all parameters don't match*/
       /*Note, we can do the atom types in both directions*/
       /*First try the dihedral one way, then the reverse since directionality is irrelevant*/
       strcpy(dihedral_make_up1,type1);
       strcat(dihedral_make_up1,type2);
       strcat(dihedral_make_up1,type3);
       strcat(dihedral_make_up1,type4);       
       strcpy(dihedral_make_up2,parm_data->dihedral_data[j].atom_type1);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type2);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type3);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type4);       
       check1=!strcmp(dihedral_make_up1,dihedral_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
       strcpy(dihedral_make_up1,type4);
       strcat(dihedral_make_up1,type3);
       strcat(dihedral_make_up1,type2);
       strcat(dihedral_make_up1,type1);
       strcpy(dihedral_make_up2,parm_data->dihedral_data[j].atom_type1);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type2);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type3);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type4);
       check2=!strcmp(dihedral_make_up1,dihedral_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
       /*if the atom type is not unique and matches here see if the parameters also match
         if they do then we have found the matching dihedral type for this dihedral so add the data and break
         else we keep looping*/
       if (check1 || check2)
       {
         if (pk == parm_data->dihedral_data[j].pk && pn == parm_data->dihedral_data[j].pn && phase == parm_data->dihedral_data[j].phase )
         {
           /*dihedral type is not unique*/
           is_unique=NO;
           /*fill in the data*/
           ++(parm_data->dihedral_data[j].number);
           /*WE NEED TO CHECK HERE THAT WE HAVEN'T OVERFLOWED OUR ARRAY*/
           if (parm_data->dihedral_data[j].number>=MAX_DIHEDRALS_PER_TYPE)
           {
             printf("*** ERROR IN PROCESS_PRMTOP - MAX_DIHEDRALS_PER_TYPE OF %d\n",MAX_DIHEDRALS_PER_TYPE);
             printf("*** EXCEEDED FOR parm_data->dihedral_data[j], j = %d\n",j);
             printf("*** EDIT MAX_DIHEDRALS_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
             return DATA_OVERFLOW;
           }
           parm_data->dihedral_data[j].atom1[(parm_data->dihedral_data[j].number)-1]=atom1;
           parm_data->dihedral_data[j].atom2[(parm_data->dihedral_data[j].number)-1]=atom2;
           parm_data->dihedral_data[j].atom3[(parm_data->dihedral_data[j].number)-1]=atom3;
           parm_data->dihedral_data[j].atom4[(parm_data->dihedral_data[j].number)-1]=atom4;
           
           break; /*Quit the loop, no point doing the rest if we have matched at this point*/
         }
       }
     } /*essentially if we finish this loop without a hit it must be unique*/
     if (is_unique==YES)
     {
       /*dihedral is unique, have to extend number of dihedral_data structures by 1*/
       ++(parm_data->unique_dihedrals_found);
       /*reallocate memory for the dihedral_data*/
       if (parm_data->unique_dihedrals_found>1)
       {
         if (global_options->VERBOSITY>=HIGH)
             printf("   Reallocating %d bytes for parm_data->*dihedral_data\n",(int ) (sizeof(struct _dihedral_data_struct)*parm_data->unique_dihedrals_found));
          parm_data->dihedral_data = (struct _dihedral_data_struct *) realloc(parm_data->dihedral_data,sizeof(struct _dihedral_data_struct)*parm_data->unique_dihedrals_found);
          if (parm_data->dihedral_data == NULL)
          {
            malloc_failure_char("process_prmtop", "parm_data->dihedral_data", sizeof(struct _dihedral_data_struct)*parm_data->unique_dihedrals_found);
            return ALLOC_FAIL;
          }
          parm_data->mem_allocated+=sizeof(struct _dihedral_data_struct);
       }
       /*now put the data in the dihedral_structure*/
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].number=1;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_KP=YES;  /*DEFAULT FOR DIHEDRALS IS TO FIT ONLY KP, this will be adjusted as necessary by other routines*/
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_NP=NO;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_PHASE=NO;       
       strcpy(parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type1,type1);
       strcpy(parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type2,type2);
       strcpy(parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type3,type3);
       strcpy(parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type4,type4);       
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].pk=pk;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].pn=pn;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].phase=phase;       
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom1[0]=atom1;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom2[0]=atom2;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom3[0]=atom3;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom4[0]=atom4;   
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].improper=is_improper;
       /*Need to check if we are to fit this parameter or not*/
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_NP=NO;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_KP=NO;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_PHASE=NO;
       
       if (global_options->RUNTYPE==SET_PARAMS)
        {
          /*prompt the user whether to fit or not*/
          if (d_kp == YES)
          {
            response='\0';
            printf("Fit Parameter: ");
            if (parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].improper==YES)
              printf("IMPROPER ");
            printf("DIHEDRAL (%s-%s-%s-%s) KP? (y/n): ",
            parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type1,
            parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type2,
            parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type3,
            parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type4);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_KP=YES;
          }
          
          if (d_np == YES)
          {
            response='\0';
            printf("Fit Parameter: ");
            if (parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].improper==YES)
              printf("IMPROPER ");
            printf("DIHEDRAL (%s-%s-%s-%s) NP? (y/n): ",
                   parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type1,
                   parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type2,
                   parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type3,
                   parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type4);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_NP=YES;
          }
          
          if(d_phase==YES)
          {
            response='\0';
            printf("Fit Parameter: ");
            if (parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].improper==YES)
            {
              printf("IMPROPER -- skipping phase\n");
              parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_PHASE=NO;
            }
            else
            {
              printf("DIHEDRAL (%s-%s-%s-%s) PHASE? (y/n): ",
                    parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type1,
                    parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type2,
                    parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type3,
                    parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type4);
              fflush(stdout);
              while( response!='y' && response!='n' )
              {
                scanf("%c",&response);
              }
              if (response=='y')
                parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_PHASE=YES;
            }
          }
        }
        else if (global_options->PARAMETERS_TO_FIT==DEFAULT)
        {
          parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_KP=YES;  /*DEFAULT FOR DIHEDRALS IS TO FIT ONLY KP, this will be adjusted as necessary by other routines*/
          parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_NP=NO;
          parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_PHASE=NO;
        }
        else if (global_options->PARAMETERS_TO_FIT == K_ONLY)
        {
          parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_KT=NO;
          parm_data->angle_data[parm_data->unique_angles_found-1].DO_FIT_THEQ=NO;
        }
     }
  } /*End of dihedrals with Hydrogen*/

  /*Now we do dihedrals withOUT hydrogen*/
  for (i=0;i<parm_data->NPHIA;++i)
  {
     /*Loop over all dihedrals not involving hydrogen*/
     /*Check if it's an improper dihedral*/
     if (parm_data->pdihedral[i].lp < 0)
       is_improper=YES;
     else
       is_improper=NO;
     
     /*Get the atoms involved in the dihedral we are looking at*/
     atom1=unObfuscateAtom(parm_data->pdihedral[i].ip);
     atom2=unObfuscateAtom(parm_data->pdihedral[i].jp);
     atom3=unObfuscateAtom(parm_data->pdihedral[i].kp);
     atom4=unObfuscateAtom(parm_data->pdihedral[i].lp);
     /*Get the parameters involved*/
     pk=parm_data->pk[parm_data->pdihedral[i].icp-1];
     pn=parm_data->pn[parm_data->pdihedral[i].icp-1];
     phase=parm_data->phase[parm_data->pdihedral[i].icp-1];
     
     /*Get the atom types*/
     strcpy(type1,parm_data->atom[atom1-1].isymbl);
     strcpy(type2,parm_data->atom[atom2-1].isymbl);
     strcpy(type3,parm_data->atom[atom3-1].isymbl);
     strcpy(type4,parm_data->atom[atom4-1].isymbl);

     /*Ok, now we have the data for this dihedral*/
     /*Now we need to determine if it is unique*/
     /*loop over all unique dihedrals to date and see if we get a match*/
     is_unique=YES;
     for (j=0;j<parm_data->unique_dihedrals_found;++j)
     {
       /*A dihedral is unique if all atom types don't match, and all parameters don't match*/
       /*Note, we can do the atom types in both directions*/
       /*First try the dihedral one way, then the reverse since directionality is irrelevant*/
       strcpy(dihedral_make_up1,type1);
       strcat(dihedral_make_up1,type2);
       strcat(dihedral_make_up1,type3);
       strcat(dihedral_make_up1,type4);
       strcpy(dihedral_make_up2,parm_data->dihedral_data[j].atom_type1);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type2);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type3);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type4);
       check1=!strcmp(dihedral_make_up1,dihedral_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
       strcpy(dihedral_make_up1,type4);
       strcat(dihedral_make_up1,type3);
       strcat(dihedral_make_up1,type2);
       strcat(dihedral_make_up1,type1);
       strcpy(dihedral_make_up2,parm_data->dihedral_data[j].atom_type1);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type2);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type3);
       strcat(dihedral_make_up2,parm_data->dihedral_data[j].atom_type4);
       check2=!strcmp(dihedral_make_up1,dihedral_make_up2);             /*strcmp returns zero on a match - so check1 will be 1 if matched*/
       /*if the atom type is not unique and matches here see if the parameters also match
         if they do then we have found the matching dihedral type for this dihedral so add the data and break
         else we keep looping*/
       if (check1 || check2)
       {
         if (pk == parm_data->dihedral_data[j].pk && pn == parm_data->dihedral_data[j].pn && phase == parm_data->dihedral_data[j].phase )
         {
           /*dihedral type is not unique*/
           is_unique=NO;
           /*fill in the data*/
           ++(parm_data->dihedral_data[j].number);
           
           /*WE NEED TO CHECK HERE THAT WE HAVEN'T OVERFLOWED OUR ARRAY*/
           if (parm_data->dihedral_data[j].number>=MAX_DIHEDRALS_PER_TYPE)
           {
             printf("*** ERROR IN PROCESS_PRMTOP - MAX_DIHEDRALS_PER_TYPE OF %d\n",MAX_DIHEDRALS_PER_TYPE);
             printf("*** EXCEEDED FOR parm_data->dihedral_data[j], j = %d\n",j);
             printf("*** EDIT MAX_DIHEDRALS_PER_TYPE IN prmtop_params.h AND RECOMPILE.\n");
             return DATA_OVERFLOW;
           }
           parm_data->dihedral_data[j].atom1[(parm_data->dihedral_data[j].number)-1]=atom1;
           parm_data->dihedral_data[j].atom2[(parm_data->dihedral_data[j].number)-1]=atom2;
           parm_data->dihedral_data[j].atom3[(parm_data->dihedral_data[j].number)-1]=atom3;
           parm_data->dihedral_data[j].atom4[(parm_data->dihedral_data[j].number)-1]=atom4;
           
           break; /*Quit the loop, no point doing the rest if we have matched at this point*/
         }
       }
     } /*essentially if we finish this loop without a hit it must be unique*/
     if (is_unique==YES)
     {
       /*dihedral is unique, have to extend number of dihedral_data structures by 1*/
       ++(parm_data->unique_dihedrals_found);
       /*reallocate memory for the dihedral_data*/
       if (parm_data->unique_dihedrals_found>1)
       {
         if (global_options->VERBOSITY>=HIGH)
             printf("   Reallocating %d bytes for parm_data->*dihedral_data\n",(int ) (sizeof(struct _dihedral_data_struct)*parm_data->unique_dihedrals_found));
          parm_data->dihedral_data = (struct _dihedral_data_struct *) realloc(parm_data->dihedral_data,sizeof(struct _dihedral_data_struct)*parm_data->unique_dihedrals_found);
          if (parm_data->dihedral_data == NULL)
          {
            malloc_failure_char("process_prmtop", "parm_data->dihedral_data", sizeof(struct _dihedral_data_struct)*parm_data->unique_dihedrals_found);
            return ALLOC_FAIL;
          }
          parm_data->mem_allocated+=sizeof(struct _dihedral_data_struct);
       }
       /*now put the data in the dihedral_structure*/
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].number=1;
       strcpy(parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type1,type1);
       strcpy(parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type2,type2);
       strcpy(parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type3,type3);
       strcpy(parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type4,type4);
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].pk=pk;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].pn=pn;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].phase=phase;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom1[0]=atom1;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom2[0]=atom2;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom3[0]=atom3;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom4[0]=atom4;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].improper=is_improper;
       
       /*Need to check if we are to fit this parameter or not*/
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_NP=NO;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_KP=NO;
       parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_PHASE=NO;

       if (global_options->RUNTYPE==SET_PARAMS)
        {
          /*prompt the user whether to fit or not*/
          if (d_kp == YES)
          {
            response='\0';
            printf("Fit Parameter: ");
						if (parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].improper==YES)
							printf("IMPROPER ");
						printf("DIHEDRAL (%s-%s-%s-%s) KP? (y/n): ",
            parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type1,
            parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type2,
            parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type3,
            parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type4);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_KP=YES;
          }
          
          if (d_np == YES)
          {
            response='\0';
            printf("Fit Parameter: ");
						if (parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].improper==YES)
							printf("IMPROPER ");
						printf("DIHEDRAL (%s-%s-%s-%s) NP? (y/n): ",
                   parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type1,
                   parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type2,
                   parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type3,
                   parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type4);
            fflush(stdout);
            while( response!='y' && response!='n' )
            {
              scanf("%c",&response);
            }
            if (response=='y')
              parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_NP=YES;
          }
          
          if (d_phase==YES)
          {
            response='\0';
            printf("Fit Parameter: ");
						if (parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].improper==YES)
            {
							printf("IMPROPER -- skipping phase\n");
              parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_PHASE=NO;
            }
            else
            {
              printf("DIHEDRAL (%s-%s-%s-%s) PHASE? (y/n): ",
                    parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type1,
                    parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type2,
                    parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type3,
                    parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].atom_type4);;
              fflush(stdout);
              while( response!='y' && response!='n' )
              {
                scanf("%c",&response);
              }
              if (response=='y')
                parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_PHASE=YES;
            }
          }

        }
        else if (global_options->PARAMETERS_TO_FIT == DEFAULT)
        {
          /*DEFAULT FOR DIHEDRALS IS TO FIT ONLY KP, this will be adjusted as necessary by other routines*/
          parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_KP=YES;
          parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_NP=NO;
          parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_PHASE=NO;
        }
        else if (global_options->PARAMETERS_TO_FIT == DEFAULT)
        {
          parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_KP=NO;
          parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_NP=NO;
          parm_data->dihedral_data[parm_data->unique_dihedrals_found-1].DO_FIT_PHASE=NO;
        }
     }
  } /*End of dihedrals withOUT Hydrogen*/
  
  // Save the inputted options to a file if desired
  if (global_options->RUNTYPE==SET_PARAMS)
  { 
    if (write_input_parameters(global_options, parm_data) != SUCCESS)
      return FAILURE;
  }
  
  // Load options if indicated
  if (global_options->PARAMETERS_TO_FIT == LOAD)
  {
    if (read_parameter_file(global_options, parm_data) != SUCCESS)
      return FAILURE;
  }
  
  if (global_options->VERBOSITY>=HIGH)
  {
    printf("Contents of parm_data->bond_data array:\n");
    for(i=0;i<parm_data->unique_bonds_found;++i)
    {
      printf("%d: number=%d, %s-%s, rk=%.4f, req=%.4f\n",i,parm_data->bond_data[i].number,parm_data->bond_data[i].atom_type1
      ,parm_data->bond_data[i].atom_type2,parm_data->bond_data[i].rk,parm_data->bond_data[i].req);
      printf("   ATOMS: ");
      for (j=0;j<parm_data->bond_data[i].number;++j)
      {
        printf("%d-%d ",parm_data->bond_data[i].atom1[j],parm_data->bond_data[i].atom2[j]);
      }
      printf("\n");
    }
    
    printf("Contents of parm_data->angle_data array:\n");
    for(i=0;i<parm_data->unique_angles_found;++i)
    {
      printf("%d: number=%d, %s-%s-%s, tk=%.4f, teq=%.4f\n",i,parm_data->angle_data[i].number,parm_data->angle_data[i].atom_type1
      ,parm_data->angle_data[i].atom_type2,parm_data->angle_data[i].atom_type3,parm_data->angle_data[i].tk,parm_data->angle_data[i].teq*RADIAN_TO_DEGREE);
      printf("   ATOMS: ");
      for (j=0;j<parm_data->angle_data[i].number;++j)
      {
        printf("%d-%d-%d ",parm_data->angle_data[i].atom1[j],parm_data->angle_data[i].atom2[j],parm_data->angle_data[i].atom3[j]);
      }
      printf("\n");
    }
    
    printf("Contents of parm_data->dihedral_data array:\n");
    for(i=0;i<parm_data->unique_dihedrals_found;++i)
    {
      printf("%d: number=%d, %s-%s-%s-%s, pk=%.4f, phase=%.4f, pn=%.4f\n",i,parm_data->dihedral_data[i].number,parm_data->dihedral_data[i].atom_type1
      ,parm_data->dihedral_data[i].atom_type2,parm_data->dihedral_data[i].atom_type3,parm_data->dihedral_data[i].atom_type4
      ,parm_data->dihedral_data[i].pk,parm_data->dihedral_data[i].phase*RADIAN_TO_DEGREE,parm_data->dihedral_data[i].pn);
      printf("   ATOMS: ");
      for (j=0;j<parm_data->dihedral_data[i].number;++j)
      {
        printf("%d-%d-%d-%d ",parm_data->dihedral_data[i].atom1[j],parm_data->dihedral_data[i].atom2[j],parm_data->dihedral_data[i].atom3[j],parm_data->dihedral_data[i].atom4[j]);
      }
      printf("\n");
    }
  }
  
  // Count how many terms there are for each dihedral
  for (i=0; i<parm_data->unique_dihedrals_found; ++i)
  {
    parm_data->dihedral_data[i].num_terms = 1;
    for (j=0; j<parm_data->unique_dihedrals_found; ++j)
    {
      if (i!=j && dihedral_types_equal(&parm_data->dihedral_data[i], &parm_data->dihedral_data[j])==YES )
        ++parm_data->dihedral_data[i].num_terms;
    }
  }
  if (global_options->VERBOSITY>=MEDIUM)
  {
    printf("   Prmtop   (unique): Found %d unique bonds.\n",parm_data->unique_bonds_found);
    printf("   Prmtop   (unique): Found %d unique angles.\n",parm_data->unique_angles_found);
    printf("   Prmtop   (unique): Found %d unique dihedrals.\n",parm_data->unique_dihedrals_found);
  }
    
  return SUCCESS;
}

