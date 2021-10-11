/*
 *  tail_interactions_debug.h
 *  
 *
 *  Created by Ognjen Perisic on 4/22/08. Corrected by Rosana Collepardo.
 *  Copyright 2009 Tamar Schlick Lab, New York University. All rights reserved. 
 *
 */

#define ALL_TAIL_PARTICLES


#ifdef ALL_TAIL_PARTICLES
  #define number_of_tails 50
#else
  #define number_of_tails 20
#endif  

#define ALL_FRAMES

#ifndef ALL_FRAMES
  #define temp_number_of_frames 5
#endif


#define TAIL_TAIL

#define multiplication_constant 1.2 

#define cut_off 1.8
#define DR_MAX 1.8
#define evd_tc 1.8
#define evd_tt 1.8
#define evd_tDNA 2.7

//#define WRITE_TAIL_DISTANCES






//#define CHECK_ROTATIONS
