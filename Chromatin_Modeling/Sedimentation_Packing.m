%Created at Tamar Schlick Lab

clear all;

 nlb=load('dim.in');
 
 cores=nlb(1);
   initial_core=1
      final_core=cores  

 
  links=nlb(2:cores+1);
  linkers=sum(links);
  S0 = 11.1; %sed coef ref for a single nucleosome without LH
  S1 = 12.0; %sed coef for a single nucleosome with LH
  R1 = 5.5; % spherical radius of the nucleosome in the sed coef formula
  LHnum = 12;
 
  LHconc = (1.*LHnum/cores);
  Sref = (S1 - S0)*LHconc + S0;
 
  
     
   [positions] = load_MC;
   
   
   
 frame=0  
   
   
lines_per_frame=cores*4+linkers*4+cores*78;
          
number_of_frames = floor(length(positions)/lines_per_frame);


FL=zeros(size(frame));
count = 0
    
 for k=1:number_of_frames;
    close all
    count = count + 1
    frame_number = k
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
   
 inend=[2 length(x)-1; 1 length(x); 2 length(x); 1 length(x)-1];
 init=inend(1,1);
last=inend(1,2);

values_x_i=[];
values_y_i=[];
values_y_i=[];


    i = (init:last)';
aa1 = csaps(i,x(init:last),.35);
values_x_i = ppval(aa1,i);                

aa2 = csaps(i,y(init:last),.35);
values_y_i = ppval(aa2,i);      

aa3 = csaps(i,z(init:last),.35);
values_z_i = ppval(aa3,i);   

seg = [];
korak = 2;          

ind=0;
for ii = 1:korak:last-init+1
    ind=ind+1;
    seg(ind,1) = values_x_i(ii);
    seg(ind,2) = values_y_i(ii);
    seg(ind,3) = values_z_i(ii);
end

FL_counter=0;   
for i=1:size(seg,1)-2
    FL_counter=FL_counter+norm(seg(i+1,:)-seg(i,:));
end



PackingRatio(k)=11*cores/FL_counter;

sedsum = 0;
 for j = 1:cores
     for k = 1:cores;
         if(j~=k);
            rj = core_positions(j,:);
            rk = core_positions(k,:);
            rjk = rk-rj;
            djk = sqrt(dot(rjk,rjk));
            sedsum = sedsum + 1/djk; 
        end
     end
end

frsc(1) = Sref*(1 + (R1/cores)*sedsum);
sc_list(count) = frsc(:);

             
 end

ave_packing_ratio = mean(PackingRatio)
sd_packing_ratio = std(PackingRatio)
ave_sed_coeff = mean(sc_list)
sd_sed_coeff = std(sc_list)
