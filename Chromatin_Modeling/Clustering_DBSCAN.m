%Created by Stephanie Portillo at Tamar Schlick Lab

clear all;

 nlb=load('dim.in');
 
 cores=nlb(1);
   initial_core=1
      final_core=cores  

 
  links=nlb(2:cores+1);
  linkers=sum(links);  
  [positions] = load_MC;
   
      
 frame=0  
   
   
lines_per_frame=cores*4+linkers*4+cores*78;
          
number_of_frames = floor(length(positions)/lines_per_frame);
starting_frame = 1;
frame_number = starting_frame;
sample_number=0;
count = 0
si =0;

while(frame_number<=number_of_frames)
           
           sample_number = sample_number + 1;
           
           if(mod(sample_number,20)==0)
               disp('.')
           end
    
%  for k=1:number_of_frames;
    close all
    count = count + 1;
frame=positions(:,(frame_number-1)*lines_per_frame+1:frame_number*lines_per_frame)';


   core_positions=zeros(cores,3);
   core_orientations=zeros(cores*3,3);
   core_index=0;
   
%calc core positions   
   for i=1:cores        
      core_index=core_index+1;
      core_positions(core_index,1:3)=frame((i-1)*4+1,1:3);
      core_orientations((core_index-1)*3+1,1:3)=frame((i-1)*4+2,1:3);
      core_orientations((core_index-1)*3+2,1:3)=frame((i-1)*4+3,1:3);
      core_orientations((core_index-1)*3+3,1:3)=frame((i-1)*4+4,1:3);        
   end
   
    for i=1:cores             
       x = core_positions(:,1);
       y = core_positions(:,2);
       z = core_positions(:,3);
    end
   
 

temp_matrix = zeros(cores,cores);
 for i = 1:cores-1
     for j = i+1:cores
         dst = sqrt((x(i)-x(j))*(x(i)-x(j)) + (y(i)-y(j))*(y(i)-y(j)) + (z(i)-z(j))*(z(i)-z(j)));
            temp_matrix(i,j) = dst;
            temp_matrix(j,i) = dst;
     end
 end
 

    Z = dbscan(temp_matrix, 20,3, 'distance', 'precomputed');
    Z1=Z(find(Z~=-1));
    new_core=length(Z1);
    
     
       
   si = si + 1;
          Max_Z1 = max(Z1);
          C_array = zeros(1, Max_Z1);
              for i = 1:new_core
              C_array(Z1(i))=C_array(Z1(i))+1;
              end
             
             Av1(si) =mean(Max_Z1); 
             Av(si) = mean(C_array);
             Max(si) = max(C_array);
          
           % ----------------------
           
           
                        
             
          frame_number = frame_number + 1
          matrix = temp_matrix;
         
     
 end

       Av_cluter = mean(Av1);
       Std_cluster =std(Av1);
       Av_nuc1 = mean(Av);
       Std_nuc1 = std(Av);
       Av_nuc2 = mean(Max);
       Std_nuc2 = std(Max);
