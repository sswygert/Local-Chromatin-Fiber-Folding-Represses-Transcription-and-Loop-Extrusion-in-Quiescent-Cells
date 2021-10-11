%Created at Tamar Schlick Lab
function [positions] = load_MC(sim)


   input_file_name = sprintf('50traj_Log.dat')
   input_file=fopen(input_file_name,'r');



positions=fscanf(input_file,'%E %E %E\n',[3 inf]);
fclose(input_file);
   

   

