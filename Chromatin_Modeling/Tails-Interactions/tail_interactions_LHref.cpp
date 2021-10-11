/////
/// tail_interactions_LHref.cpp
/// Adapted by Antoni Luque on 02/25/2014 
/// The number of LH beads are generalized in this code
/// All the frames are now explored (before only the first frame was stored)
///
///
/// tail_interactions.cpp
/// Adapted by Antoni Luque based on Rosana Collepardo's code
/*
 *  tail_interactons_1.cpp
 *  
 *
 *  Created by Ognjen Persic on 4/22/09. Corrected by Rosana Collepardo.
 *  Copyright 2009 Tamar Schlick Lab, New York University. All rights reserved.
 *
 */
 
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <iomanip>
#include "tail_interactions_debug.h"
/////////
/// Modification: TONI
/// This library was missing and in some machines c++ does not load it by default
#include <cstring>
/////////


using namespace std;

const int number_of_core_particles = 300;  // number of DiSCO charges
const int MAX_NUMBER_OF_CORES = 300;        // number of cores
const int MAX_NUMBER_OF_LINKERS = 1000;       // number of linker beads total
const int NUMBER_OF_TAILS = 50;            // number of tail beads per core
const int NUMBER_OF_LH = 28;                // number of linker histone beads 

int n_c;                    // number of cores;
int n_l;                    // number of linker DNA beads (total);
int n_b;                    // number of linker DNA beads (ave per linker); 
int n_t;                    // number of tail beads;
int l_h_p;                  // linker histone presence;
int n_l_h;                  // number of linker histone beads;
int lines_per_frame;        // n_c*4 + n_l*4 + n_t + n_l_h;
int trajectory_number;      // simulation trajectory;
int omit_frames;            // number of frames to ommit
int max_number_of_frames;   // max number of frames
bool read_o;
bool read_m;
int step;


float frame_line[3];
float core_particles[number_of_core_particles][3];
float core_charges[number_of_core_particles];
float core_coordinates[MAX_NUMBER_OF_CORES][3];
float core_orientations[MAX_NUMBER_OF_CORES*3][3];
float dna_coordinates[MAX_NUMBER_OF_CORES*MAX_NUMBER_OF_LINKERS][3];
float dna_orientations[MAX_NUMBER_OF_CORES*MAX_NUMBER_OF_LINKERS*3][3];
float tail_coordinates[NUMBER_OF_TAILS*MAX_NUMBER_OF_CORES][3];
//////////
/// Modification: TONI (LHref)
float linker_histone_coordinates[NUMBER_OF_LH][3];
// float linker_histone_coordinates[MAX_NUMBER_OF_CORES*3][3];
//////////

//////////
/// Modification: TONI (LHref)
bool ifi,whilei;
//////////

int free_inter_H3_1;          // free tails H3
int free_inter_H3_2;          // free tails H3

int free_inter_H4_1;          // free tails H4
int free_inter_H4_2;          // free tails H4

int free_inter_H2A1_1;        // free tails H2A1
int free_inter_H2A1_2;        // free tails H2A1

int free_inter_H2B_1;         // free tails H2B
int free_inter_H2B_2;         // free tails H2B

int free_inter_H2A2_1;        // free tails H2A2
int free_inter_H2A2_2;        // free tails H2A2



int tails_inter_H3_1;          // interactions between tails H3 and all other tails not belonging to its parent nucleosome
int tails_inter_H3_2;          // interactions between tails H3 and all other tails not belonging to its parent nucleosome

int tails_inter_H4_1;          // interactions between tails H4 and all other tails not belonging to its parent nucleosome
int tails_inter_H4_2;          // interactions between tails H4 and all other tails not belonging to its parent nucleosome

int tails_inter_H2A1_1;        // interactions between tails H2A1 and all other tails not belonging to its parent nucleosome 
int tails_inter_H2A1_2;        // interactions between tails H2A1 and all other tails not belonging to its parent nucleosome 

int tails_inter_H2B_1;         // interactions between tails H2B and all other tails not belonging to its parent nucleosome 
int tails_inter_H2B_2;         // interactions between tails H2B and all other tails not belonging to its parent nucleosome 

int tails_inter_H2A2_1;        // interactions between tails H2A2 and all other tails not belonging to its parent nucleosome
int tails_inter_H2A2_2;        // interactions between tails H2A2 and all other tails not belonging to its parent nucleosome


int cores_inter_H3_1;         // interactions between tails H3 and all cores except its parent core
int temp_cores_inter_H3_1;    // temp interactions between tails H3 and all cores except its parent core
int cores_inter_H3_2;         // interactions between tails H3 and all cores except its parent core
int temp_cores_inter_H3_2;    // temp interactions between tails H3 and all cores except its parent core

int cores_inter_H4_1;         // interactions between tails H4 and all cores except its parent core
int cores_inter_H4_2;         // interactions between tails H4 and all cores except its parent core

int cores_inter_H2A1_1;       // interactions between tails H2A1 and all cores except its parent core
int cores_inter_H2A1_2;       // interactions between tails H2A1 and all cores except its parent core

int cores_inter_H2B_1;        // interactions between tails H2B and all cores except its parent core
int cores_inter_H2B_2;        // interactions between tails H2B and all cores except its parent core

int cores_inter_H2A2_1;       // interactions between tails H2A2 and all cores except its parent core
int cores_inter_H2A2_2;       // interactions between tails H2A2 and all cores except its parent core


int DNA_inter_H3_1;         // interactions between tails H3 and all DNA except its parent DNA
int DNA_inter_H3_2;         // interactions between tails H3 and all DNA except its parent DNA

int DNA_inter_H4_1;         // interactions between tails H4 and all DNA except its parent DNA
int DNA_inter_H4_2;         // interactions between tails H4 and all DNA except its parent DNA

int DNA_inter_H2A1_1;       // interactions between tails H2A1 and all DNA except its parent DNA
int DNA_inter_H2A1_2;       // interactions between tails H2A1 and all DNA except its parent DNA

int DNA_inter_H2B_1;        // interactions between tails H2B and all DNA except its parent DNA
int DNA_inter_H2B_2;        // interactions between tails H2B and all DNA except its parent DNA

int DNA_inter_H2A2_1;       // interactions between tails H2A2 and all DNA except its parent DNA
int DNA_inter_H2A2_2;       // interactions between tails H2A2 and all DNA except its parent DNA


int parent_DNA_inter_H3_1;         // interactions between tails H3 and its parent DNA
int parent_DNA_inter_H3_2;         // interactions between tails H3 and its parent DNA

int parent_DNA_inter_H4_1;         // interactions between tails H4 and its parent DNA
int parent_DNA_inter_H4_2;         // interactions between tails H4 and its parent DNA

int parent_DNA_inter_H2A1_1;       // interactions between tails H2A1 and its parent DNA
int parent_DNA_inter_H2A1_2;       // interactions between tails H2A1 and its parent DNA

int parent_DNA_inter_H2B_1;        // interactions between tails H2B and its parent DNA
int parent_DNA_inter_H2B_2;        // interactions between tails H2B and its parent DNA

int parent_DNA_inter_H2A2_1;       // interactions between tails H2A2 and its parent DNA
int parent_DNA_inter_H2A2_2;       // interactions between tails H2A2 and its parent DNA


int parent_core_inter_H3_1;         // interactions between tails H3 and their parent core
int parent_core_inter_H3_2;         // interactions between tails H3 and their parent core

int parent_core_inter_H4_1;         // interactions between tails H4 and their parent core
int parent_core_inter_H4_2;         // interactions between tails H4 and their parent core

int parent_core_inter_H2A1_1;       // interactions between tails H2A1 and their parent core
int parent_core_inter_H2A1_2;       // interactions between tails H2A1 and their parent core

int parent_core_inter_H2B_1;        // interactions between tails H2B and their parent core
int parent_core_inter_H2B_2;        // interactions between tails H2B and their parent core

int parent_core_inter_H2A2_1;       // interactions between tails H2A2 and their parent core
int parent_core_inter_H2A2_2;       // interactions between tails H2A2 and their parent core



int scopyn(char *t, char *s, int n)
{
    int i;	
	
    if(n<=strlen(s))
    {
        for(i=0; i<n; i++)
        {				
            *(t+i) = *(s+i);			
        }
        t[i++]='\0';
        return 0;
    }
    return -1;
}


int scopyn2n(char *t, char *s, int n1, int n2)
{
    int i,j;	
    if((n2<=strlen(s))&(n1>=0)&(n1<=n2))
    {
        j = 0;
        for(i=n1-1; i<n2; i++)
        {				
            *(t+j++) = *(s+i);			
        }
        t[++j]='\0';
        return 0;
    }
    return -1;
}


int command_line_interpreter(int argc, char* argv[])
{
    int i, temp_number;	
    char tmp[80];
    bool read_c = false, read_l = false, read_h = false, read_t = false, read_step = false;
    
    read_o = false;
	read_m = false;
	step = 1;
	
	max_number_of_frames = 10000;

    if(argc>2)
    {
        for (i=1; i<argc; i++)	
        {
            if(!scopyn(tmp,argv[i],2))
            {
                
                if(strcmp(tmp,"-c")==0)
                {				
                    scopyn2n(tmp,argv[i],3,strlen(argv[i]));
                    temp_number = atoi(tmp);					
                    if((temp_number<=0)|(temp_number>MAX_NUMBER_OF_CORES))
                    {
                        printf("Wrong number of cores (max %d)!\n",MAX_NUMBER_OF_CORES);
                        return 1;
                    }
                    else
                    {
                        n_c = temp_number;				
                        n_t = n_c*50;
                        read_c = true;
                    }
                }
				
                if(strcmp(tmp,"-l")==0)
                {				
                    scopyn2n(tmp,argv[i],3,strlen(argv[i]));
                    temp_number = atoi(tmp);
                    
                    if((temp_number<=0)|(temp_number>MAX_NUMBER_OF_LINKERS))
                    {
                        printf("Wrong number of linker beads (max %d)!\n",MAX_NUMBER_OF_LINKERS);
                        return 1;
                    }
                    else
                    {
                        n_l = temp_number;
                        read_l = true;
                    }
                }
			    			

                if(strcmp(tmp,"-h")==0)
                {				
                    scopyn2n(tmp,argv[i],3,strlen(argv[i]));
                    temp_number = atoi(tmp);
                    if((temp_number<0)|(temp_number>1))
                    {
                        printf("Wrong linker histone definition (max = 1)!\n");
                        return -1;
                    }
                    else
		    {
                        l_h_p = temp_number;
			//////////
			/// Modification: TONI (LHref)
			n_l_h = NUMBER_OF_LH;				
			//n_l_h = 3;				
                        read_h = true;
			//////////
                    }
                }
				

                if(strcmp(tmp,"-t")==0)
                {				
                    scopyn2n(tmp,argv[i],3,strlen(argv[i]));
                    temp_number = atoi(tmp);
                    if((temp_number<0)|(temp_number>8))
                    {
                        printf("Wrong trajectory number (max = 8)!\n");
		        return -1;
		    }
		    else
		    {
                         trajectory_number = temp_number;	
			 printf("Trajectory number : %d\n",trajectory_number);								
			 read_t = true;
		    }
                }
				
				
                if(strcmp(tmp,"-o")==0)
                {
                     scopyn2n(tmp,argv[i],3,strlen(argv[i]));
                     temp_number = atoi(tmp);
                     if(temp_number<=0)
                     {
                         printf("Wrong number of frames to be omited!\n");
                         return 1;
                     }
                     else
                     {
                         omit_frames = temp_number;
                         printf("omit number : %d\n",temp_number);
                         read_o = true;
                     }
                        
                 }
				 
				 
                if(strcmp(tmp,"-e")==0)
                {				
                    scopyn2n(tmp,argv[i],3,strlen(argv[i]));
                    temp_number = atoi(tmp);					
                    if((temp_number<=0)|(temp_number>1000))
                    {
                        printf("Wrong step number (0<step<=1000)!\n");
                        return 1;
                    }
                    else
                    {					    
                        step = temp_number;				
						printf("Step = %d\n",step);
                        //n_t = n_c*50;
                        read_step = true;
                    }
                }

					
					
                if(strcmp(tmp,"-m")==0)
                {
                         scopyn2n(tmp,argv[i],3,strlen(argv[i]));
                         temp_number = atoi(tmp);
                         if(temp_number<=0)
                         {
                            printf("Wrong maximum number of frames!\n");
                            return 1;
                         }
                         else
                         {
                            max_number_of_frames = temp_number;
                            printf("max. number : %d\n",temp_number);
                            read_o = true;
                         }
                    }
                }
                else
                    return 1;
        }
    }
    else
    {
        printf("\n Tail interactions \n\n");
        printf(" tail_interactions_2 -[t s l h m e o][number]\n\n");
        printf(" -c[number] = number of chromatin cores\n");
        printf(" -l[number] = total number of linker beads\n");
        printf(" -t[number] = trajectory number 1 = inter_1.dat, etc.\n");
        printf(" -h[number] = linker histone presence; 0 = -LH; 1 = +LH;\n");
        printf("              h SHOULD ALWAYS BE SET TO 1 NOW!\n");
        printf(" -m[number] = maximum number of frames to be analyzed; (default = 10,000)\n");
	printf(" -e[number] = frame step value; (default = 1)\n");
	printf(" -o[number] = number of initial steps to be omitted!\n");
		
        printf("\n Example: > ");
        printf("tail_interactions_2 -c24 -t4 -l125 -t1\n");
        printf(" Analyses 4th trajectory (fort.80) with 24 cores, 6 linker beads and linker histone (LH) per simulation frame.\n\n");
        printf("\n (c), Ogi, 2009.\n");
		printf("    additional changes by gBascom 2016.\n");
        printf("    now works with refined LH and nonuniform NRLs\n\n");

        return 1;
    }
	
    
    if(read_c&read_l&read_t&read_h)
    {	    
	//////////
	/// Modification: TONI (LHref)
	lines_per_frame = n_c*4 + n_l*4 + n_t + l_h_p*n_c*n_l_h;
        // lines_per_frame = n_c*4 + n_l*4 + n_t + l_h_p*n_c*3;
	//////////		
        printf("Lines per frame %d\n",lines_per_frame);
		
        printf("\n");
        return 0;
    }
    else
    {
        printf("Wrong input!\n");
        return 1;
    }
}


int fill_core_particles(char* file_name)
{
    ifstream core_file(file_name, ios::in );
    int i = 0;

    if( core_file == NULL )
    {
        printf( "The file %s was not opened\n",file_name );
        return 0;
    }
    else
    {
        while ((core_file >> core_particles[i][0] >> core_particles[i][1] >> core_particles[i][2])&&(i<300))
        {
	    core_particles[i][0] = core_particles[i][0] /10; 
	    core_particles[i][1] = core_particles[i][1] /10; 
            core_particles[i][2] = core_particles[i][2] /10; 
            i++;
        }

        i=0;
        while (core_file >> core_charges[i])
        {
            i++;
        }
        core_file.close();
        return 0;
    }
}


int min_distance(double drs[5], int types[5])
{
     int i, j, temp_type; // 100 = parent core, 101 = other cores, 102 = other tails 
     double temp_dr;
	 
	 for (i=0;i< 5-1; i++)
	 {
	    for (j=i+1;j<5;j++)
		{
		   if(drs[j]<=drs[i])
		   {
		      temp_dr = drs[i];
	              drs[i]  = drs[j];
		      drs[j]  = temp_dr;
			  
		      temp_type = types[i];
		      types[i]  = types[j];
		      types[j]  = temp_type; 
		   }
		}
	 }

    return types[0];
}

int main(int argc, char *argv[])
{
    char simulation_file_name[256];
    char file_name[100];
    char file_directory[512];
    int current_line;
    int core_orientation_position = 0;
    int core_coordinate_position = 0;
    long int total_line_number = 0;
    int dna_coordinate_index = 0;
    int dna_orientation_index = 0;
    int tail_coordinate_index = 0;
    int simulation_frame = 0;
    int i, j, k, l, m, it, local_core_orientation_position, la, lb, lc, icp, t, it2;
    int i_dna, j_dna, k_dna, start_counter;	
    long int totdata = 0;
    static int tailE [number_of_tails];
    float cut_off_value = cut_off;
    float dx, dy, dz, dr;
	
    
    bool core_finished = false;
    
    double dr_core = DR_MAX;
    double core_tail_distance = 10000;
	
    double dr_parent_core = DR_MAX;
    double parent_core_tail_distance = 10000;
	
    double dr_tail = DR_MAX;
    double tail_tail_distance = 10000;
	
    double dr_DNA = DR_MAX;
    double DNA_tail_distance = 10000;
	
    double dr_parent_DNA = DR_MAX;	
    double parent_DNA_tail_distance = 10000;
		
    double drs  [5];
    int    types[5];
	
    int interaction_type;
    int totdata_H3 = 0, totdata_H4 = 0, totdata_H2A1 = 0, totdata_H2A2 = 0, totdata_H2B = 0;


    if(number_of_tails==10)
    {
        tailE[0] = 0;  // H3
        tailE[1] = 8;  // H3
 	tailE[2] = 16; // H4
	tailE[3] = 21; // H4 
	tailE[4] = 26; // H2A1 
	tailE[5] = 30; // H2A1
	tailE[6] = 34; // H2B
	tailE[7] = 39; // H2B
	tailE[8] = 44; // H2A2
	tailE[9] = 47; // H2A2
    }

    if(number_of_tails==20)
    {
        tailE[0]  = 0;  // H3
        tailE[1]  = 1;  // H3
        tailE[2]  = 8;  // H3		
        tailE[3]  = 9;  // H3
		
        tailE[4]  = 16; // H4
        tailE[5]  = 17; // H4 
        tailE[6]  = 21; // H4
        tailE[7]  = 22; // H4 
		
        tailE[8]  = 26; // H2A1 
        tailE[9]  = 27; // H2A1
        tailE[10] = 30; // H2A1 
        tailE[11] = 31; // H2A1
		

        tailE[12] = 34; // H2B
        tailE[13] = 35; // H2B
        tailE[14] = 39; // H2B
        tailE[15] = 40; // H2B
		
        tailE[16] = 44; // H2A2
        tailE[17] = 45; // H2A2
        tailE[18] = 47; // H2A2
        tailE[19] = 48; // H2A2

    }
	
	
	
    if(number_of_tails==50)
    {
        tailE[0]  = 0;  // H3
        tailE[1]  = 1;  // H3
	tailE[2]  = 2;  // H3
	tailE[3]  = 3;  // H3
	tailE[4]  = 4;  // H3
	tailE[5]  = 5;  // H3
	tailE[6]  = 6;  // H3
	tailE[7]  = 7;  // H3
		
        tailE[8]  =  8;  // H3		
        tailE[9]  =  9;  // H3
        tailE[10] = 10;  // H3
        tailE[11] = 11;  // H3
        tailE[12] = 12;  // H3
        tailE[13] = 13;  // H3
        tailE[14] = 14;  // H3
        tailE[15] = 15;  // H3								
		
																		
        tailE[16] = 16; // H4
        tailE[17] = 17; // H4 
        tailE[18] = 18; // H4 
        tailE[19] = 19; // H4 
        tailE[20] = 20; // H4 						
		
        tailE[21] = 21; // H4
        tailE[22] = 22; // H4 
        tailE[23] = 23; // H4 
        tailE[24] = 24; // H4 				
        tailE[25] = 25; // H4 						
		
		
        tailE[26]  = 26; // H2A1 
        tailE[27]  = 27; // H2A1
        tailE[28]  = 28; // H2A1
        tailE[29]  = 29; // H2A1				
		
        tailE[30] = 30; // H2A1 
        tailE[31] = 31; // H2A1
        tailE[32] = 32; // H2A1
        tailE[33] = 33; // H2A1				
		

        tailE[34] = 34; // H2B
        tailE[35] = 35; // H2B
        tailE[36] = 36; // H2B
        tailE[37] = 37; // H2B
        tailE[38] = 38; // H2B
		
        tailE[39] = 39; // H2B
        tailE[40] = 40; // H2B
        tailE[41] = 41; // H2B
        tailE[42] = 42; // H2B
        tailE[43] = 43; // H2B
		
		
        tailE[44] = 44; // H2A2
        tailE[45] = 45; // H2A2
        tailE[46] = 46; // H2A2
		
        tailE[47] = 47; // H2A2
        tailE[48] = 48; // H2A2
        tailE[49] = 49; // H2A2

    }

	
	
	
	
    cores_inter_H3_1 = 0;         
    cores_inter_H3_2 = 0;         

    cores_inter_H4_1 = 0;         
    cores_inter_H4_2 = 0;         

    cores_inter_H2A1_1 = 0;       
    cores_inter_H2A1_2 = 0;       

    cores_inter_H2B_1 = 0;        
    cores_inter_H2B_2 = 0;        

    cores_inter_H2A2_1 = 0;       
    cores_inter_H2A2_2 = 0;       

    parent_core_inter_H3_1 = 0;         
    parent_core_inter_H3_2 = 0;         

    parent_core_inter_H4_1 = 0;         
    parent_core_inter_H4_2 = 0;         

    parent_core_inter_H2A1_1 = 0;       
    parent_core_inter_H2A1_2 = 0;       

    parent_core_inter_H2B_1 = 0;        
    parent_core_inter_H2B_2 = 0;        

    parent_core_inter_H2A2_1 = 0;       
    parent_core_inter_H2A2_2 = 0;       
	   	   
	   
    DNA_inter_H3_1 = 0;
    DNA_inter_H3_2 = 0;

    DNA_inter_H4_1 = 0;
    DNA_inter_H4_2 = 0;

    DNA_inter_H2A1_1 = 0;
    DNA_inter_H2A1_2 = 0;

    DNA_inter_H2B_1 = 0;
    DNA_inter_H2B_2 = 0;

    DNA_inter_H2A2_1 = 0;
    DNA_inter_H2A2_2 = 0;


    parent_DNA_inter_H3_1 = 0;
    parent_DNA_inter_H3_2 = 0;

    parent_DNA_inter_H4_1 = 0;
    parent_DNA_inter_H4_2 = 0;

    parent_DNA_inter_H2A1_1 = 0;
    parent_DNA_inter_H2A1_2 = 0;

    parent_DNA_inter_H2B_1 = 0;
    parent_DNA_inter_H2B_2 = 0;

    parent_DNA_inter_H2A2_1 = 0;
    parent_DNA_inter_H2A2_2 = 0;


    tails_inter_H3_1 = 0;          
    tails_inter_H3_2 = 0;          


    tails_inter_H4_1 = 0;          
    tails_inter_H4_2 = 0;          

    tails_inter_H2A1_1 = 0;        
    tails_inter_H2A1_2 = 0;        

    tails_inter_H2B_1 = 0;         
    tails_inter_H2B_2 = 0;         

    tails_inter_H2A2_1 = 0;        
    tails_inter_H2A2_2 = 0;        


    free_inter_H3_1 = 0;          
    free_inter_H3_2 = 0;          

    free_inter_H4_1 = 0;          
    free_inter_H4_2 = 0;          

    free_inter_H2A1_1 = 0;        
    free_inter_H2A1_2 = 0;        

    free_inter_H2B_1 = 0;         
    free_inter_H2B_2 = 0;         

    free_inter_H2A2_1 = 0;        
    free_inter_H2A2_2 = 0;        
  		   
	
    printf("\n\n\nSTART!\n");
	
    if(!command_line_interpreter(argc,argv))
    {

      //////////////////////
      /// Modification: TONI
      /// The file inter.dat is in the running directory
      sprintf(simulation_file_name,"inter_%d.dat",trajectory_number);
      /// Old
      //sprintf(simulation_file_name,"../inter.dat",trajectory_number);
      //////////////////////
       fill_core_particles("core_data.reg.150mM");
		
       printf("Here we go!\n");

       ifstream simulation_file(simulation_file_name, ios::in );
				
       printf("File %s is opened.\n",simulation_file_name);
		
       printf("Max. number of frames %d in analysis\n",max_number_of_frames);

        current_line = 0;
	  #ifdef CHECK_ROTATIONS
		ofstream core_orientations_file("core_orientations_1.txt");
	  #endif



#ifdef WRITE_TAIL_DISTANCES

/*    ------------------------ H3 distances -----------------------------------------------    */	  
      sprintf(file_name,"parent_core_tail_H3_distances_t%d.txt",trajectory_number);		  
	  ofstream parent_core_H3_distances(file_name);
	  
	  sprintf(file_name,"core_tail_H3_distances_t%d.txt",trajectory_number);		  
	  ofstream core_H3_distances(file_name);
	  
	  sprintf(file_name,"parent_DNA_tail_H3_distances_t%d.txt",trajectory_number);		  		  	  
	  ofstream parent_DNA_H3_distances(file_name);

	  sprintf(file_name,"DNA_tail_H3_distances_t%d.txt",trajectory_number);		  		  	  	  
	  ofstream DNA_H3_distances(file_name);


/*    ------------------------ H4 distances -----------------------------------------------    */	  
      sprintf(file_name,"parent_core_tail_H4_distances_t%d.txt",trajectory_number);		  
	  ofstream parent_core_H4_distances(file_name);
	  
	  sprintf(file_name,"core_tail_H4_distances_t%d.txt",trajectory_number);		  
	  ofstream core_H4_distances(file_name);
	  
	  sprintf(file_name,"parent_DNA_tail_H4_distances_t%d.txt",trajectory_number);		  		  	  
	  ofstream parent_DNA_H4_distances(file_name);

	  sprintf(file_name,"DNA_tail_H4_distances_t%d.txt",trajectory_number);		  		  	  	  
	  ofstream DNA_H4_distances(file_name);


/*    ------------------------ H2A1 distances ---------------------------------------------    */	  
      sprintf(file_name,"parent_core_tail_H2A1_distances_t%d.txt",trajectory_number);		  
	  ofstream parent_core_H2A1_distances(file_name);
	  
	  sprintf(file_name,"core_tail_H2A1_distances_t%d.txt",trajectory_number);		  
	  ofstream core_H2A1_distances(file_name);
	  
	  sprintf(file_name,"parent_DNA_tail_H2A1_distances_t%d.txt",trajectory_number);		  		  	  
	  ofstream parent_DNA_H2A1_distances(file_name);

	  sprintf(file_name,"DNA_tail_H2A1_distances_t%d.txt",trajectory_number);		  		  	  	  
	  ofstream DNA_H2A1_distances(file_name);


/*    ------------------------ H2B distances -----------------------------------------------    */	  
      sprintf(file_name,"parent_core_tail_H2B_distances_t%d.txt",trajectory_number);		  
	  ofstream parent_core_H2B_distances(file_name);
	  
	  sprintf(file_name,"core_tail_H2B_distances_t%d.txt",trajectory_number);		  
	  ofstream core_H2B_distances(file_name);
	  
	  sprintf(file_name,"parent_DNA_tail_H2B_distances_t%d.txt",trajectory_number);		  		  	  
	  ofstream parent_DNA_H2B_distances(file_name);

	  sprintf(file_name,"DNA_tail_H2B_distances_t%d.txt",trajectory_number);		  		  	  	  
	  ofstream DNA_H2B_distances(file_name);

/*    ------------------------ H2A2 distances -----------------------------------------------    */	  
      sprintf(file_name,"parent_core_tail_H2A2_distances_t%d.txt",trajectory_number);		  
	  ofstream parent_core_H2A2_distances(file_name);
	  
	  sprintf(file_name,"core_tail_H2A2_distances_t%d.txt",trajectory_number);		  
	  ofstream core_H2A2_distances(file_name);
	  
	  sprintf(file_name,"parent_DNA_tail_H2A2_distances_t%d.txt",trajectory_number);		  		  	  
	  ofstream parent_DNA_H2A2_distances(file_name);

	  sprintf(file_name,"DNA_tail_H2A2_distances_t%d.txt",trajectory_number);		  		  	  	  
	  ofstream DNA_H2A2_distances(file_name);
#endif
		
	  //////////////
	  /// Modification: TONI (LHref)
	  /// The macro ALL_FRAMES was giving problems when running over all frames.
	  /// This while condition works better
	  while(!simulation_file.eof()&&(simulation_frame<max_number_of_frames))
	  //      #ifdef ALL_FRAMES
	  //	  while ((simulation_file >> frame_line[0] >> frame_line[1] >> frame_line[2]))
	  //      #else
	  //	  while ((simulation_file >> frame_line[0] >> frame_line[1] >> frame_line[2])&&(simulation_frame<max_number_of_frames))
	  //      #endif	
	  //////////////
	  {
	    //////////////
	    /// Modification: TONI (LHref)
	    /// We store the coordinates of the line
	    /// We double check if we're at the end of the file (the eof flag is activated once the stream reads after the last line)
	    simulation_file >> frame_line[0] >> frame_line[1] >> frame_line[2];
	    if(!simulation_file.eof()) /// This ensures that the last line is not reread
	    {
	    /// printf("check %d %d %f %f %f\n",simulation_frame,current_line,frame_line[0],frame_line[1],frame_line[2]);	  	  
	    //////////////
            // --------------------------- NUCLEOSOME cores --------------------------------------------------
            if(current_line<n_c*4)
            {
                if((current_line % 4)==0)
                {
                    core_coordinates[core_coordinate_position][0]=frame_line[0];
                    core_coordinates[core_coordinate_position][1]=frame_line[1];
                    core_coordinates[core_coordinate_position][2]=frame_line[2];

		    core_coordinate_position++;
                }

                if(((current_line-1) % 4)==0)
                {
                    core_orientations[core_orientation_position][0]=frame_line[0];
                    core_orientations[core_orientation_position][1]=frame_line[1];
                    core_orientations[core_orientation_position][2]=frame_line[2];
                    core_orientation_position++;
                }

                if(((current_line-2) % 4)==0)
                {
                    core_orientations[core_orientation_position][0]=frame_line[0];
                    core_orientations[core_orientation_position][1]=frame_line[1];
                    core_orientations[core_orientation_position][2]=frame_line[2];
                    core_orientation_position++;
                }

                if(((current_line-3) % 4)==0)
                {
                    core_orientations[core_orientation_position][0]=frame_line[0];
                    core_orientations[core_orientation_position][1]=frame_line[1];
                    core_orientations[core_orientation_position][2]=frame_line[2];
                    core_orientation_position++;
                }
            }
           
            // ------------------ DNA beads ---------------------------------------------------
            if((current_line>=n_c*4)&&(current_line<(n_c*4+n_l*4)))
            {
                if((current_line % 4)==0)
                {
                    dna_coordinates[dna_coordinate_index][0]=frame_line[0];
                    dna_coordinates[dna_coordinate_index][1]=frame_line[1];
                    dna_coordinates[dna_coordinate_index][2]=frame_line[2];
                    dna_coordinate_index++;
                }

                if(((current_line-1) % 4)==0)
                {
                    dna_orientations[dna_orientation_index][0]=frame_line[0];
                    dna_orientations[dna_orientation_index][1]=frame_line[1];
                    dna_orientations[dna_orientation_index][2]=frame_line[2];
                    dna_orientation_index++;
                }

                if(((current_line-2) % 4)==0)
                {
                    dna_orientations[dna_orientation_index][0]=frame_line[0];
                    dna_orientations[dna_orientation_index][1]=frame_line[1];
                    dna_orientations[dna_orientation_index][2]=frame_line[2];
                    dna_orientation_index++;
                }

                if(((current_line-3) % 4)==0)
                {
                    dna_orientations[dna_orientation_index][0]=frame_line[0];
                    dna_orientations[dna_orientation_index][1]=frame_line[1];
                    dna_orientations[dna_orientation_index][2]=frame_line[2];
                    dna_orientation_index++;
                }
                
            }
			
            // ------------------ HISTONE tails ---------------------------------------------------
            if((current_line>=(n_c*4+n_l*4))&&(current_line<(n_c*4+n_l*4+n_t)))
            {
                tail_coordinates[tail_coordinate_index][0]=frame_line[0];
                tail_coordinates[tail_coordinate_index][1]=frame_line[1];
                tail_coordinates[tail_coordinate_index][2]=frame_line[2];
                tail_coordinate_index++;
            }
						            
            current_line++;
            total_line_number++;
	    //////////////
	    /// Modification: TONI (LHref)
	    //printf("Frame and line %d %d\n",simulation_frame,current_line);
	    //printf("%d %f \n",core_coordinate_position-1,core_coordinates[core_coordinate_position-1][0]);
	    //////////////

            if((current_line % (lines_per_frame))==0)
            {
                			
                simulation_frame++;
		
		//////////////
		/// Modification: TONI (LHref)
		/// We reset current_line and the main indexes 
		/// In the old code current_line kept incrementing and only the first frame we properly stored
		current_line = 0;
		core_coordinate_position = 0;
		core_orientation_position = 0;
		dna_coordinate_index = 0;
		dna_orientation_index = 0;
		tail_coordinate_index = 0;
		//////////////
		
                if((simulation_frame % 2) == 0)
                    printf("sim_frame=%d omit_frames=%d  max_number_of_frames=%d\n",simulation_frame,omit_frames,max_number_of_frames);
							
               if((((read_o&(simulation_frame>omit_frames))|(!read_o))&(simulation_frame<=max_number_of_frames))&&(((simulation_frame+1) % step)==0))
               {
			   
				
                    if(((simulation_frame+1) % 2)==0)
                            printf("Analyzing frame %d\n",simulation_frame);
				
				
				    /* ------ pass through all cores ------- */
				    
                    for(j = 0; j<n_c; j++)
                    {
					
					
                      /* ------------ pass through tails belongint to the core j --------- */
					    
                        for (k = 0; k<number_of_tails; k++)
                        {	
						
						
                             dr_tail = 40*evd_tt;  
                             DNA_tail_distance = 40*evd_tt;
						
                             dr_core = 40*evd_tc; 
                             core_tail_distance = 40*evd_tc; 
						
                             dr_parent_core = 40*evd_tc; 
                             parent_core_tail_distance = 40*evd_tc;
						
                             dr_DNA = 40*evd_tDNA;    
                             parent_DNA_tail_distance = 40*evd_tDNA; 
						
                             dr_parent_DNA = 40*evd_tDNA;  
                             DNA_tail_distance = 40*evd_tDNA;			
						
						
                             //totdata++;
							 
							     
                             #ifndef ALL_TAIL_PARTICLES
                                it = j*50 + tailE[k];
                             #else 	 
                                it = j*50 + k;
                             #endif	 
                             	 							 
				
				/* --------------------------------------------- Analyze DNA ------------------------------------------- */
				/*NOTE NEED TO DOUBLE CHECK THAT THIS IS DOING WHAT IT IS SUPPOSED TO*/
				n_b=n_l/n_c;
				
				for (m = 0; m<n_l;m++)
				  {
				    
				    dx = tail_coordinates[it][0] - dna_coordinates[m][0];
				    dy = tail_coordinates[it][1] - dna_coordinates[m][1];
				    dz = tail_coordinates[it][2] - dna_coordinates[m][2];		
                                
				    dr = sqrt(dx*dx + dy*dy + dz*dz);
								
				    if((j==int(floor(m/n_b))||((j-1)==int(floor(m/n_b)))))  // neighbouring DNA linkers
				      {
					{
					  if(dr<dr_parent_DNA)	
					    {
					      dr_parent_DNA = dr;
					      
					      parent_DNA_tail_distance = dr - evd_tDNA;
					      
					      if(parent_DNA_tail_distance < 0)
						parent_DNA_tail_distance = 0;
					    }		
					}							
					
				      }
				    else  // other linkers
				      {
					{
					  if(dr<dr_DNA)
					    {
					      dr_DNA = dr;
					      
					      DNA_tail_distance = dr - evd_tDNA;
					      
					      if (DNA_tail_distance < 0)
						DNA_tail_distance = 0;
					    }
					}		   
				      }//if(j== neighbouting DNa or not!
				  }//for(m=0; m<n_l ; m++)
						


	 	/* ----------- Analyze other cores and their interaction with tail 'it'  ------------- */
							 
                            for(m = 0 ; m<n_c; m++)
                            {		

			      //Tail-tail interactions
				if(j!=m)
				  {
				    for(t=0;(t<number_of_tails);t++)
				      {
#ifdef ALL_TAIL_PARTICLES
					it2 = m*50 + t;
#else
					it2 = m*50 + tailE[t];
#endif
					
					dx = tail_coordinates[it][0] - tail_coordinates[it2][0];
					dy = tail_coordinates[it][1] - tail_coordinates[it2][1];
					dz = tail_coordinates[it][2] - tail_coordinates[it2][2];
					
					dr = sqrt(dx*dx + dy*dy + dz*dz);
											
					//if(dr>multiplication_constant*evd_tt)
					{
					  if(dr<dr_tail)
					    {
					      dr_tail = dr;
					      tail_tail_distance = dr - evd_tt;
					      if(tail_tail_distance < 0)
						tail_tail_distance = 0;
					    }	  
					}
				      }//end for t=0 t<number_of_tails;t++
	
				    //INTERACTION WITH OTHER CORES
							     			   					       
                                    la = m*3;
                                    lb = m*3 + 1;
                                    lc = m*3 + 2;
                                    icp = m;
                                    
                                    core_finished = false;
                                          
				    for(t = 0; (t<300); t++)
                                    {
				      dx = tail_coordinates[it][0] - core_coordinates[icp][0] - core_orientations[la][0]*core_particles[t][0] - core_orientations[lb][0]*core_particles[t][1] - core_orientations[lc][0]*core_particles[t][2];
				      dy = tail_coordinates[it][1] - core_coordinates[icp][1] - core_orientations[la][1]*core_particles[t][0] - core_orientations[lb][1]*core_particles[t][1] - core_orientations[lc][1]*core_particles[t][2];
				      dz = tail_coordinates[it][2] - core_coordinates[icp][2] - core_orientations[la][2]*core_particles[t][0] - core_orientations[lb][2]*core_particles[t][1] - core_orientations[lc][2]*core_particles[t][2];
				      
#ifdef CHECK_ROTATIONS
				      core_orientations_file << core_orientations[la][0]<< " " << core_orientations[lb][0] << " " << core_orientations[lc][0] << endl;
				      core_orientations_file << core_orientations[la][1]<< " " << core_orientations[lb][1] << " " << core_orientations[lc][1] << endl;
				      core_orientations_file << core_orientations[la][2]<< " " << core_orientations[lb][2] << " " << core_orientations[lc][2] << endl;
#endif
					  
				      dr = sqrt(dx*dx + dy*dy + dz*dz);	
										
				      {
					if(dr<dr_core)
					  {
					    dr_core = dr;
					    
					    core_tail_distance = dr - evd_tc;
					    if(core_tail_distance < 0)
					      core_tail_distance = 0;
					  }	  
				      }
										
                                    }     // end - for(t = 0; t<300; t++) 
				  }      // end - if(j!=m)
                                else
                                {   
				  // here j=m
				  // --------- interaction with parent nucleosome --------------------------------
				  
				  				  if(((k!=7)||(k!=6)||(k!=5))||((k!=15)||(k!=14)||(k!=13))||((k!=20)||(k!=19)||(k!=18))||((k!=25)||(k!=24)||(k!=23))||((k!=29)||(k!=28))||((k!=33)||(k!=32))||((k!=38)||(k!=37)||(k!=36))||((k!=43)||(k!=42)||(k!=41))||(k==46)||(k==49))
				// if((k!=7)||(k!=15)||(k!=20)||(k!=25)||(k!=29)||(k!=33)||(k!=38)||(k!=43)||(k==46)||(k==49))
				    {
				      core_finished = false;
				      la = m*3;
				      lb = m*3 + 1;
				      lc = m*3 + 2;
				      icp = m;
                                   
				      
				      for(t = 0; (t<300);t++)
					{
					  dx = tail_coordinates[it][0] - core_coordinates[icp][0] - core_orientations[la][0]*core_particles[t][0] - core_orientations[lb][0]*core_particles[t][1] - core_orientations[lc][0]*core_particles[t][2];
					  dy = tail_coordinates[it][1] - core_coordinates[icp][1] - core_orientations[la][1]*core_particles[t][0] - core_orientations[lb][1]*core_particles[t][1] - core_orientations[lc][1]*core_particles[t][2];
					  dz = tail_coordinates[it][2] - core_coordinates[icp][2] - core_orientations[la][2]*core_particles[t][0] - core_orientations[lb][2]*core_particles[t][1] - core_orientations[lc][2]*core_particles[t][2];
					  
					  
					  dr = sqrt(dx*dx + dy*dy + dz*dz);	
					  //if(dr<multiplication_constant*evd_tc)
					  { 
					    if(dr<dr_parent_core)
					      {
						dr_parent_core = dr;
						parent_core_tail_distance = dr - evd_tc;
						
						if(parent_core_tail_distance < 0)
						  parent_core_tail_distance = 0;
													
												
					      }
					  }
					}     // end - for(t = 0; t<300; tl++) 
				    }         // if(((k!=7)||(k!=6)||(k!=5))||((k!=15)||(k!=14)||
                                }            // else

                            }     // end - for(m = 0 ; m<n_c; m++)   / - other cores analysis
							
							
							
			    drs[0] = core_tail_distance;
			    drs[1] = parent_core_tail_distance;
			    drs[2] = DNA_tail_distance;
			    drs[3] = parent_DNA_tail_distance;
			    drs[4] = tail_tail_distance;
							
			    
			    types[0] = 100;  // dr_core;
			    types[1] = 101;  // dr_parent_core;
			    types[2] = 102;  // dr_DNA;
			    types[3] = 103;  // dr_parent_DNA;
			    types[4] = 105;  // dr_tail;
							
			    
			    if((core_tail_distance<=multiplication_constant*evd_tc)||(parent_core_tail_distance<=multiplication_constant*evd_tc)||(DNA_tail_distance<=(multiplication_constant*evd_tDNA))||(parent_DNA_tail_distance<=(multiplication_constant*evd_tDNA))||(tail_tail_distance<multiplication_constant*evd_tt))
			      {
				interaction_type = min_distance(drs, types);								 
			      }
			    else
			      {
				interaction_type = 200;  // free tail							   							   
			      }
							
			    
                            if((tailE[k]>=0)&&(tailE[k]<=7))
			      {
							
				totdata_H3++;
				
				switch(interaction_type)
				  {
				  case 100:cores_inter_H3_1++;
				    break;		   
				  case 101:parent_core_inter_H3_1++;
				    break;					
                                  case 102:DNA_inter_H3_1++;
				    break;
                                  case 103:parent_DNA_inter_H3_1++;
				    break;
                                  case 105:tails_inter_H3_1++;
				    break;
                                  case 200:free_inter_H3_1++;
				    break;										   
				  }
#ifdef WRITE_TAIL_DISTANCES
				parent_core_H3_distances << setw(10) << setprecision(3) << (parent_core_tail_distance + evd_tc) << endl;
				core_H3_distances        << setw(10) << setprecision(3) << (       core_tail_distance + evd_tc) << endl;
				parent_DNA_H3_distances  << setw(10) << setprecision(3) << (parent_DNA_tail_distance  + evd_tDNA + 0.9) << endl;
				DNA_H3_distances         << setw(10) << setprecision(3) << (       DNA_tail_distance  + evd_tDNA + 0.9)  << endl;
#endif
			      }
							
							
                           if((tailE[k]>=8)&&(tailE[k]<=15))
                            {
			      
			      totdata_H3++;
			      
			      switch(interaction_type)
				{
				case 100:cores_inter_H3_2++;
				  break;		   
				case 101:parent_core_inter_H3_2++;
				  break;					
				case 102:DNA_inter_H3_2++;
				  break;
				case 103:parent_DNA_inter_H3_2++;
				  break;
				case 105:tails_inter_H3_2++;
				  break;
				case 200:free_inter_H3_2++;
				  break;										   
                                }
#ifdef WRITE_TAIL_DISTANCES
			      parent_core_H3_distances << setw(10) << setprecision(3) << (parent_core_tail_distance + evd_tc) << endl;
			      core_H3_distances        << setw(10) << setprecision(3) << (       core_tail_distance + evd_tc) << endl;
			      parent_DNA_H3_distances  << setw(10) << setprecision(3) << (parent_DNA_tail_distance  + evd_tDNA + 0.9) << endl;
			      DNA_H3_distances         << setw(10) << setprecision(3) << (       DNA_tail_distance  + evd_tDNA + 0.9)  << endl;
#endif
			      
			    }

			   
                             if((tailE[k]>=16)&&(tailE[k]<=20))
			       {
				 
				 totdata_H4++;
				 
				 switch(interaction_type)
				   {
				   case 100:cores_inter_H4_1++;
				     break;		   
				   case 101:parent_core_inter_H4_1++;
				     break;				
				   case 102:DNA_inter_H4_1++;
				     break;
				   case 103:parent_DNA_inter_H4_1++;
				     break;										   						   
				   case 105:tails_inter_H4_1++;
				     break;
				   case 200:free_inter_H4_1++;
				     break;										   
				   }
								
#ifdef WRITE_TAIL_DISTANCES
				 parent_core_H4_distances << setw(10) << setprecision(3) << (parent_core_tail_distance + evd_tc) << endl;
				 core_H4_distances        << setw(10) << setprecision(3) << (       core_tail_distance + evd_tc) << endl;
				 parent_DNA_H4_distances  << setw(10) << setprecision(3) << (parent_DNA_tail_distance  + evd_tDNA + 0.9) << endl;
				 DNA_H4_distances         << setw(10) << setprecision(3) << (       DNA_tail_distance  + evd_tDNA + 0.9)  << endl;
#endif
								
				 
			       }
							 
                             if((tailE[k]>=21)&&(tailE[k]<=25))
                             {
			       
			       totdata_H4++;
			       
			       switch(interaction_type)
				 {
				 case 100:cores_inter_H4_2++;
				   break;		   
				 case 101:parent_core_inter_H4_2++;
				   break;				
				 case 102:DNA_inter_H4_2++;
				   break;
				 case 103:parent_DNA_inter_H4_2++;
				   break;										   						   
				 case 105:tails_inter_H4_2++;
				   break;
				 case 200:free_inter_H4_2++;
				   break;										   
				 }
			       
#ifdef WRITE_TAIL_DISTANCES
			       parent_core_H4_distances << setw(10) << setprecision(3) << (parent_core_tail_distance + evd_tc) << endl;
			       core_H4_distances        << setw(10) << setprecision(3) << (       core_tail_distance + evd_tc) << endl;
			       parent_DNA_H4_distances  << setw(10) << setprecision(3) << (parent_DNA_tail_distance  + evd_tDNA + 0.9) << endl;
			       DNA_H4_distances         << setw(10) << setprecision(3) << (       DNA_tail_distance  + evd_tDNA + 0.9)  << endl;
#endif
                             }

			     
                             if((tailE[k]>=26)&&(tailE[k]<=29))
			       {
							 
				 totdata_H2A1++;
				 
				 switch(interaction_type)
				   {
				   case 100:cores_inter_H2A1_1++;
				     break;		   
				   case 101:parent_core_inter_H2A1_1++;
				     break;
				   case 102:DNA_inter_H2A1_1++;
				     break;
				   case 103:parent_DNA_inter_H2A1_1++;
				     break;										   										   
				   case 105:tails_inter_H2A1_1++;
				     break;
				   case 200:free_inter_H2A1_1++;
				     break;
				   }
								
#ifdef WRITE_TAIL_DISTANCES
				 parent_core_H2A1_distances << setw(10) << setprecision(3) << (parent_core_tail_distance + evd_tc) << endl;
				 core_H2A1_distances        << setw(10) << setprecision(3) << (       core_tail_distance + evd_tc) << endl;
				 parent_DNA_H2A1_distances  << setw(10) << setprecision(3) << (parent_DNA_tail_distance  + evd_tDNA + 0.9) << endl;
				 DNA_H2A1_distances         << setw(10) << setprecision(3) << (       DNA_tail_distance  + evd_tDNA + 0.9)  << endl;
#endif
			       }											
			     
                             if((tailE[k]>=30)&&(tailE[k]<=33))
			       {
							 
				 totdata_H2A1++;
				 
				 switch(interaction_type)
				   {
				   case 100:cores_inter_H2A1_2++;
				     break;		   
				   case 101:parent_core_inter_H2A1_2++;
				     break;
				   case 102:DNA_inter_H2A1_2++;
				     break;
				   case 103:parent_DNA_inter_H2A1_2++;
				     break;										   										   
				   case 105:tails_inter_H2A1_2++;
				     break;
				   case 200:free_inter_H2A1_2++;
				     break;
				   }
#ifdef WRITE_TAIL_DISTANCES
				 parent_core_H2A1_distances << setw(10) << setprecision(3) << (parent_core_tail_distance + evd_tc) << endl;
				 core_H2A1_distances        << setw(10) << setprecision(3) << (       core_tail_distance + evd_tc) << endl;
				 parent_DNA_H2A1_distances  << setw(10) << setprecision(3) << (parent_DNA_tail_distance  + evd_tDNA + 0.9) << endl;
				 DNA_H2A1_distances         << setw(10) << setprecision(3) << (       DNA_tail_distance  + evd_tDNA + 0.9)  << endl;
#endif
			       }											
							 
							 
							 
										
			     if((tailE[k]>=34)&&(tailE[k]<=38))
			       {
				 
				 totdata_H2B++;
				 
				 switch(interaction_type)
				   {
				   case 100:cores_inter_H2B_1++;
				     break;		   
				   case 101:parent_core_inter_H2B_1++;
				     break;
				   case 102:DNA_inter_H2B_1++;
				     break;
				   case 103:parent_DNA_inter_H2B_1++;
				     break;										   										   
				   case 105:tails_inter_H2B_1++;
				     break;
				   case 200:free_inter_H2B_1++;
				     break;										   
				   }
#ifdef WRITE_TAIL_DISTANCES
				 parent_core_H2B_distances << setw(10) << setprecision(3) << (parent_core_tail_distance + evd_tc) << endl;
				 core_H2B_distances        << setw(10) << setprecision(3) << (       core_tail_distance + evd_tc) << endl;
				 parent_DNA_H2B_distances  << setw(10) << setprecision(3) << (parent_DNA_tail_distance  + evd_tDNA + 0.9) << endl;
				 DNA_H2B_distances         << setw(10) << setprecision(3) << (       DNA_tail_distance  + evd_tDNA + 0.9) << endl;
#endif
                             }
			     
			     if((tailE[k]>=39)&&(tailE[k]<=43))
			       {
				 
				 totdata_H2B++;
				 
				 switch(interaction_type)
				   {
				   case 100:cores_inter_H2B_2++;
				     break;		   
				   case 101:parent_core_inter_H2B_2++;
				     break;
				   case 102:DNA_inter_H2B_2++;
				     break;
				   case 103:parent_DNA_inter_H2B_2++;
				     break;										   										   
				   case 105:tails_inter_H2B_2++;
				     break;
				   case 200:free_inter_H2B_2++;
				     break;										   
				   }
#ifdef WRITE_TAIL_DISTANCES
				 parent_core_H2B_distances << setw(10) << setprecision(3) << (parent_core_tail_distance + evd_tc) << endl;
				 core_H2B_distances        << setw(10) << setprecision(3) << (       core_tail_distance + evd_tc) << endl;
				 parent_DNA_H2B_distances  << setw(10) << setprecision(3) << (parent_DNA_tail_distance  + evd_tDNA + 0.9) << endl;
				 DNA_H2B_distances         << setw(10) << setprecision(3) << (       DNA_tail_distance  + evd_tDNA + 0.9) << endl;
#endif
			       }
			     
							 
			     
                             if((tailE[k]>=44)&&(tailE[k]<=46))
			       {			
				 
				 totdata_H2A2++;
				 
				 switch(interaction_type)
				   {
				   case 100:cores_inter_H2A2_1++;
				     break;		   
				   case 101:parent_core_inter_H2A2_1++;
				     break;
				   case 102:DNA_inter_H2A2_1++;
				     break;
				   case 103:parent_DNA_inter_H2A2_1++;
				     break;										   										   
				   case 105:tails_inter_H2A2_1++;
				     break;
				   case 200:free_inter_H2A2_1++;
				     break;										   
				   }
#ifdef WRITE_TAIL_DISTANCES
				 parent_core_H2A2_distances << setw(10) << setprecision(3) << (parent_core_tail_distance + evd_tc) << endl;
				 core_H2A2_distances        << setw(10) << setprecision(3) << (       core_tail_distance + evd_tc) << endl;
				 parent_DNA_H2A2_distances  << setw(10) << setprecision(3) << (parent_DNA_tail_distance  + evd_tDNA + 0.9) << endl;
				 DNA_H2A2_distances         << setw(10) << setprecision(3) << (       DNA_tail_distance  + evd_tDNA + 0.9)  << endl;							
#endif
			       }
			     
							 
                             if((tailE[k]>=47)&&(tailE[k]<=49))
			       {				
							 
				 totdata_H2A2++;
				 
				 
				 switch(interaction_type)
				   {
				   case 100:cores_inter_H2A2_2++;
				     break;		   
				   case 101:parent_core_inter_H2A2_2++;
				     break;
				   case 102:DNA_inter_H2A2_2++;
				     break;
				   case 103:parent_DNA_inter_H2A2_2++;
				     break;										   										   
				   case 105:tails_inter_H2A2_2++;
				     break;
				   case 200:free_inter_H2A2_2++;
				     break;										   
				   }
#ifdef WRITE_TAIL_DISTANCES
				 parent_core_H2A2_distances << setw(10) << setprecision(3) << (parent_core_tail_distance + evd_tc) << endl;
				 core_H2A2_distances        << setw(10) << setprecision(3) << (       core_tail_distance + evd_tc) << endl;
				 parent_DNA_H2A2_distances  << setw(10) << setprecision(3) << (parent_DNA_tail_distance  + evd_tDNA + 0.9) << endl;
				 DNA_H2A2_distances         << setw(10) << setprecision(3) << (       DNA_tail_distance  + evd_tDNA + 0.9)  << endl;							
#endif
			       }
			     

	
							
                        }     // end - for (k = 0; k<10; k++)   /  - tails from the core
                    }      //  end -  for(j = 0; j<n_c; j++)		
		    
					/*
					  cout << "free_inter_H3   = "  << 100*(free_inter_H3_1 + free_inter_H3_2)/(totdata/5)  << endl;
					  cout << "free_inter_H4   = "  << 100*(free_inter_H4_1 + free_inter_H4_2)/(totdata/5) << endl;
					  cout << "free_inter_H2A1 = "  << 100*(free_inter_H2A1_1 + free_inter_H2A1_2)/(totdata/5) << endl;
					  cout << "free_inter_H2B  = "  << 100*(free_inter_H2B_1 + free_inter_H2B_2)/(totdata/5) << endl;							
					  cout << "free_inter_H2A2 = "  << 100*(free_inter_H2A2_1 + free_inter_H2A2_2)/(totdata/5) << endl;							
					  std::cin.get();
					*/

					
                }  // end if((simulation_frame>omit_frames)				
                else
                {
                    if((simulation_frame % 10)==0)
                        printf("Frame %d omitted!\n",simulation_frame);
                }
            }      // end - if((current_line % (lines_per_frame))==0)
	    //////////
	    /// Modification: TONI (LHref)
	  }  // end - if(!simulation_file.eof())	  
	    //////////
	}	   // end - while ((simulation_file >> frame[current_line][0] >> frame[current_line][1] >> frame[current_line][2])&&(totdata<1400)) 


		simulation_file.close();
       
	   #ifdef WRITE_TAIL_DISTANCES
  	    parent_core_H3_distances.close();
		core_H3_distances.close();
		parent_DNA_H3_distances.close();
		DNA_H3_distances.close();

  	    parent_core_H4_distances.close();
		core_H4_distances.close();
		parent_DNA_H4_distances.close();
		DNA_H4_distances.close();

  	    parent_core_H2A1_distances.close();
		core_H2A1_distances.close();
		parent_DNA_H2A1_distances.close();
		DNA_H2A1_distances.close();

  	    parent_core_H2B_distances.close();
		core_H2B_distances.close();
		parent_DNA_H2B_distances.close();
		DNA_H2B_distances.close();

  	    parent_core_H2A2_distances.close();
		core_H2A2_distances.close();
		parent_DNA_H2A2_distances.close();
		DNA_H2A2_distances.close();
	   #endif 
					
    }
	
	printf("Number of frames : %d\n",simulation_frame);
	
	#ifdef CHECK_ROTATIONS
	//core_orientations_file.close();
	#endif
	
	
	//for(i=0;i<2;i++)
	//{
	//
	//}

	
	if(simulation_frame>1)
	{
	
	
	    printf("totdata = %d\n",totdata);
		
	    printf("Interactions with tails!\n");
	
		sprintf(file_name,"tail_tail_interactions_H3_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H3_file(file_name);
		//tail_interactions_H3_file << setw(10) << setprecision(3) << float(float(tails_inter_H3_1 + tails_inter_H3_2)/(float(totdata)/50)) << endl;   			
		tail_interactions_H3_file << setw(10) << setprecision(3) << float(float(tails_inter_H3_1 + tails_inter_H3_2)/(float(totdata_H3))) << endl;   			
		tail_interactions_H3_file.close();
		
		
		sprintf(file_name,"tail_tail_interactions_H4_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H4_file(file_name);
		tail_interactions_H4_file << setw(10) << setprecision(3) << float(float(tails_inter_H4_1 + tails_inter_H4_2)/(float(totdata_H4))) << endl;   
		tail_interactions_H4_file.close();


		sprintf(file_name,"tail_tail_interactions_H2A1_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2A1_file(file_name);
		///////////
		/// Modification: TONI
		/// The end line was missing
		tail_interactions_H2A1_file << setw(10) << setprecision(3) << float(float(tails_inter_H2A1_1 + tails_inter_H2A1_2)/(float(totdata_H2A1))) << endl;   
		//tail_interactions_H2A1_file << setw(10) << setprecision(3) << float(float(tails_inter_H2A1_1 + tails_inter_H2A1_2)/(float(totdata_H2A1))) << " ";   
		///////////
		tail_interactions_H2A1_file.close();


		sprintf(file_name,"tail_tail_interactions_H2B_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2B_file(file_name);
		tail_interactions_H2B_file << setw(10) << setprecision(3) << float(float(tails_inter_H2B_1 + tails_inter_H2B_2)/(float(totdata_H2B))) << endl;   
		tail_interactions_H2B_file.close();


		sprintf(file_name,"tail_tail_interactions_H2A2_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2A2_file(file_name);
		tail_interactions_H2A2_file << setw(10) << setprecision(3) << float(float(tails_inter_H2A2_1 + tails_inter_H2A2_2)/(float(totdata_H2A2))) << endl;   
		tail_interactions_H2A2_file.close();





        /* --------------------------------- Free tails! ----------------------------- */ 

	    printf("Free tails!\n");
	
		sprintf(file_name,"free_tails_H3_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream free_tail_interactions_H3_file(file_name);
		//free_tail_interactions_H3_file << setw(10) << setprecision(3) << float(float(free_inter_H3_1 + free_inter_H3_2)/(float(totdata)/50)) << endl;  								
		free_tail_interactions_H3_file << setw(10) << setprecision(3) << float(float(free_inter_H3_1 + free_inter_H3_2)/(float(totdata_H3))) << endl;
		free_tail_interactions_H3_file.close();
		
		
		sprintf(file_name,"free_tails_H4_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream free_tail_interactions_H4_file(file_name);
		free_tail_interactions_H4_file << setw(10) << setprecision(3) << float(float(free_inter_H4_1 + free_inter_H4_2)/(float(totdata_H4))) << endl;   
		free_tail_interactions_H4_file.close();


		sprintf(file_name,"free_tails_H2A1_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream free_tail_interactions_H2A1_file(file_name);
		///////////
		/// Modification: TONI
		/// The end line was missing
		free_tail_interactions_H2A1_file << setw(10) << setprecision(3) << float(float(free_inter_H2A1_1 + free_inter_H2A1_2)/(float(totdata_H2A1))) << endl;   
		//free_tail_interactions_H2A1_file << setw(10) << setprecision(3) << float(float(free_inter_H2A1_1 + free_inter_H2A1_2)/(float(totdata_H2A1))) << " ";   
		free_tail_interactions_H2A1_file.close();


		sprintf(file_name,"free_tails_H2B_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream free_tail_interactions_H2B_file(file_name);
		free_tail_interactions_H2B_file << setw(10) << setprecision(3) << float(float(free_inter_H2B_1 + free_inter_H2B_2)/(float(totdata_H2B))) << endl;   
		free_tail_interactions_H2B_file.close();


		sprintf(file_name,"free_tails_H2A2_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream free_tail_interactions_H2A2_file(file_name);
		free_tail_interactions_H2A2_file << setw(10) << setprecision(3) << float(float(free_inter_H2A2_1 + free_inter_H2A2_2)/(float(totdata_H2A2))) << endl;   
		free_tail_interactions_H2A2_file.close();


		
				

        /* ------------- Interactions with cores ------------------- */
				
		printf("Interactions with cores!\n");
		
		sprintf(file_name,"tail_core_interactions_H3_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H3_C_file(file_name);
		//tail_interactions_H3_C_file << setw(10) << setprecision(3) << float(float(cores_inter_H3_1 + cores_inter_H3_2)/(float(totdata)/50)) << endl;                   
		tail_interactions_H3_C_file << setw(10) << setprecision(3) << float(float(cores_inter_H3_1 + cores_inter_H3_2)/(float(totdata_H3))) << endl;                   
		tail_interactions_H3_C_file.close();
		
		
		sprintf(file_name,"tail_core_interactions_H4_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H4_C_file(file_name);
		tail_interactions_H4_C_file << setw(10) << setprecision(3) << float(float(cores_inter_H4_1 + cores_inter_H4_2)/(float(totdata_H4))) << endl;   
		tail_interactions_H4_C_file.close();


		sprintf(file_name,"tail_core_interactions_H2A1_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2A1_C_file(file_name);
		tail_interactions_H2A1_C_file << setw(10) << setprecision(3) << float(float(cores_inter_H2A1_1 + cores_inter_H2A1_2)/(float(totdata_H2A1))) << endl;   
		tail_interactions_H2A1_C_file.close();


		sprintf(file_name,"tail_core_interactions_H2B_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2B_C_file(file_name);
		tail_interactions_H2B_C_file << setw(10) << setprecision(3) << float(float(cores_inter_H2B_1 + cores_inter_H2B_2)/(float(totdata_H2B))) << endl;   
		tail_interactions_H2B_C_file.close();


		sprintf(file_name,"tail_core_interactions_H2A2_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2A2_C_file(file_name);
		tail_interactions_H2A2_C_file << setw(10) << setprecision(3) << float(float(cores_inter_H2A2_1 + cores_inter_H2A2_2)/(float(totdata_H2A2))) << endl;   
		tail_interactions_H2A2_C_file.close();



        /* ----------------- Interactions with parent core ------------------- */
		
		printf("Interactions with parent cores!\n");
		
		sprintf(file_name,"tail_parent_core_interactions_H3_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H3_PC_file(file_name);
		//tail_interactions_H3_PC_file << setw(10) << setprecision(3) << float(float(parent_core_inter_H3_1 + parent_core_inter_H3_2)/(float(totdata)/50)) << endl;                   
		tail_interactions_H3_PC_file << setw(10) << setprecision(3) << float(float(parent_core_inter_H3_1 + parent_core_inter_H3_2)/(float(totdata_H3))) << endl;                   
		tail_interactions_H3_PC_file.close();
		
		
		sprintf(file_name,"tail_parent_core_interactions_H4_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H4_PC_file(file_name);
		tail_interactions_H4_PC_file << setw(10) << setprecision(3) << float(float(parent_core_inter_H4_1 + parent_core_inter_H4_2)/(float(totdata_H4))) << endl;   
		tail_interactions_H4_PC_file.close();


		sprintf(file_name,"tail_parent_core_interactions_H2A1_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2A1_PC_file(file_name);
		tail_interactions_H2A1_PC_file << setw(10) << setprecision(3) << float(float(parent_core_inter_H2A1_1 + parent_core_inter_H2A1_2)/(float(totdata_H2A1))) << endl;   
		tail_interactions_H2A1_PC_file.close();


		sprintf(file_name,"tail_parent_core_interactions_H2B_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2B_PC_file(file_name);
		tail_interactions_H2B_PC_file << setw(10) << setprecision(3) << float(float(parent_core_inter_H2B_1 + parent_core_inter_H2B_2)/(float(totdata_H2B))) << endl;   
		tail_interactions_H2B_PC_file.close();


		sprintf(file_name,"tail_parent_core_interactions_H2A2_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2A2_PC_file(file_name);
		tail_interactions_H2A2_PC_file << setw(10) << setprecision(3) << float(float(parent_core_inter_H2A2_1 + parent_core_inter_H2A2_2)/(float(totdata_H2A2))) << endl;   
		tail_interactions_H2A2_PC_file.close();
			
		
		
		/* ------------- Interactions with DNA ------------------- */
				
		printf("Interactions with DNA!\n");
		
		sprintf(file_name,"tail_DNA_interactions_H3_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H3_DNA_file(file_name);
		//tail_interactions_H3_DNA_file << setw(10) << setprecision(3) << float(float(DNA_inter_H3_1 + DNA_inter_H3_2)/(float(totdata)/50)) << endl;                   
		tail_interactions_H3_DNA_file << setw(10) << setprecision(3) << float(float(DNA_inter_H3_1 + DNA_inter_H3_2)/(float(totdata_H3))) << endl;
		tail_interactions_H3_DNA_file.close();
		
		
		sprintf(file_name,"tail_DNA_interactions_H4_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H4_DNA_file(file_name);
		tail_interactions_H4_DNA_file << setw(10) << setprecision(3) << float(float(DNA_inter_H4_1 + DNA_inter_H4_2)/(float(totdata_H4))) << endl;   
		tail_interactions_H4_DNA_file.close();


		sprintf(file_name,"tail_DNA_interactions_H2A1_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2A1_DNA_file(file_name);
		tail_interactions_H2A1_DNA_file << setw(10) << setprecision(3) << float(float(DNA_inter_H2A1_1 + DNA_inter_H2A1_2)/(float(totdata_H2A1))) << endl;   
		tail_interactions_H2A1_DNA_file.close();


		sprintf(file_name,"tail_DNA_interactions_H2B_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2B_DNA_file(file_name);
		tail_interactions_H2B_DNA_file << setw(10) << setprecision(3) << float(float(DNA_inter_H2B_1 + DNA_inter_H2B_2)/(float(totdata_H2B))) << endl;   
		tail_interactions_H2B_DNA_file.close();


		sprintf(file_name,"tail_DNA_interactions_H2A2_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2A2_DNA_file(file_name);
		tail_interactions_H2A2_DNA_file << setw(10) << setprecision(3) << float(float(DNA_inter_H2A2_1 + DNA_inter_H2A2_2)/(float(totdata_H2A2))) << endl;   
		tail_interactions_H2A2_DNA_file.close();

		
		

		/* ------------- Interactions with parent DNA ------------------- */
				
		printf("Interactions with parent DNA!\n");
		
		sprintf(file_name,"tail_parent_DNA_interactions_H3_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H3_parent_DNA_file(file_name);
		//tail_interactions_H3_parent_DNA_file << setw(10) << setprecision(3) << float(float(parent_DNA_inter_H3_1 + parent_DNA_inter_H3_2)/(float(totdata)/50)) << endl;
		tail_interactions_H3_parent_DNA_file << setw(10) << setprecision(3) << float(float(parent_DNA_inter_H3_1 + parent_DNA_inter_H3_2)/(float(totdata_H3))) << endl;
		tail_interactions_H3_DNA_file.close();
		
		
		sprintf(file_name,"tail_parent_DNA_interactions_H4_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H4_parent_DNA_file(file_name);
		tail_interactions_H4_parent_DNA_file << setw(10) << setprecision(3) << float(float(parent_DNA_inter_H4_1 + parent_DNA_inter_H4_2)/(float(totdata_H4))) << endl;   
		tail_interactions_H4_parent_DNA_file.close();


		sprintf(file_name,"tail_parent_DNA_interactions_H2A1_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2A1_parent_DNA_file(file_name);
		tail_interactions_H2A1_parent_DNA_file << setw(10) << setprecision(3) << float(float(parent_DNA_inter_H2A1_1 + parent_DNA_inter_H2A1_2)/(float(totdata_H2A1))) << endl;   
		tail_interactions_H2A1_parent_DNA_file.close();


		sprintf(file_name,"tail_parent_DNA_interactions_H2B_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2B_parent_DNA_file(file_name);
		tail_interactions_H2B_parent_DNA_file << setw(10) << setprecision(3) << float(float(parent_DNA_inter_H2B_1 + parent_DNA_inter_H2B_2)/(float(totdata_H2B))) << endl;   
		tail_interactions_H2B_parent_DNA_file.close();


		sprintf(file_name,"tail_parent_DNA_interactions_H2A2_t%d_",trajectory_number);	
		sprintf(file_name,"%s%1.1f.txt",file_name,cut_off_value);
		ofstream tail_interactions_H2A2_parent_DNA_file(file_name);
		tail_interactions_H2A2_parent_DNA_file << setw(10) << setprecision(3) << float(float(parent_DNA_inter_H2A2_1 + parent_DNA_inter_H2A2_2)/(float(totdata_H2A2))) << endl;   
		tail_interactions_H2A2_parent_DNA_file.close();
		
		

	}
	
		//printf("%s\n",file_directory);
	
	return 0;
}



