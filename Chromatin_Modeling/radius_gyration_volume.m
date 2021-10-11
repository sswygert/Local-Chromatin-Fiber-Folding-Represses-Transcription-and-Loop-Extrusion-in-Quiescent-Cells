% % cell: Notes
% This m-file computes the radius of gyration (center of mass and central
% axis) of the cores
% for a given frame-sequence of a chromatin system.
% Created by Gavin Bascom at Tamar Schlick Lab

% The radius of gyration is computed from the standard formula:
% Rg^2 = 1/N sum_k{rk - rmean}
% This is the definition used for the distribution (and average within)

% The gyration tensorclear all is computed from the standard formula:
% Smn = 1/N sum_i{ri_m*ri_n}
% The volume is estimated using the AlphaShape function native to matlab
% NOTE: You should be graphing and checking the volumes drawn in that they
% can be sensitive to the value of alpha chosen. A high alpha will give the
% general volume of the entire shape (no filled in holes) whereas a low
% alpha will five more or less the cylindrical volume and is insensitive to
% global folding.
%clear all
load_data=1;
% % cell: input
%get number of cores
 nlb=load('dim.in');
 cores=nlb(1);
     
 %get linker length list
    links=nlb(2:cores+1);
  linkers=sum(links);
 

     sim=1;
     if (load_data)
   [positions] = load_MC(sim);
     end

%calc number of frames 
lines_per_frame = cores*4+linkers*4+cores*78
number_of_frames = floor(length(positions)/lines_per_frame);


 Rgave = zeros(cores,1);
 Rgstd = zeros(cores,1);
% 
 volcellave = zeros(cores,1);
 volcellstd = zeros(cores,1);
 volume_alphashape = zeros(cores,1);
 volume_alphashape_std = zeros(cores,1);
 position_array = zeros(cores*50,3);

yy_counter=0; 
linker_counter = 0;
bp_counter=0;
for yy=4:cores
end_bead = yy

%Load all frames, 
%in each successive frame, fill core_positions and dna_linker_position arrays,
%calc r_gyr and volume estimates
%load into list arrays
%end loop
%calc averages and stdev's

yy_counter=yy_counter+1;
%calc the genome position (in kb)
bp_counter=bp_counter+1;
linker_counter=linker_counter+nlb(yy_counter);
bp_index(bp_counter)=(((yy_counter*16)+linker_counter)*9)/1000;

for zzz=1:number_of_frames
    frame_number=zzz;
     frame=positions(:,(frame_number-1)*lines_per_frame+1:frame_number*lines_per_frame)';
  core_positions=zeros(cores,3);
   core_index=0;
   
   % ------------------- Nucleosome cores ------------------------------------------------------
   
   for i=1:cores        
      core_index=core_index+1;
      core_positions(core_index,1:3)=frame((i-1)*4+1,1:3);
      
   end
   
   % ---------------------------------------------------------------------------------------
   
   % ------------------ DNA linkers --------------------------------------------------------
   dna_linker_positions=zeros(linkers,3);
   dna_linker_index=0;
   for i=1:linkers
      dna_linker_index=dna_linker_index+1;
      poz=cores*4;
      dna_linker_positions(dna_linker_index,1:3)=frame(poz+(i-1)*4+1,1:3);
      
   end
   % ------------------TAILS---------------------------------------------------------------------
   
    poz = cores*4;
   i   = linkers;   
   current_position = poz + (i-1)*4 + 4 + 1;   
   end_tail = (current_position - 1) + cores*50;      
   tails = zeros(cores*50,3);
   tail_counter = 0;
   for i=current_position:end_tail
       tail_counter = tail_counter + 1;
       tails(tail_counter,1) = frame(i,1);
       tails(tail_counter,2) = frame(i,2);
       tails(tail_counter,3) = frame(i,3);       
   
   
   
   if(mod(tail_counter-1,50)==0)
            tailint=0;
        end
        if(tailint>=47)
            tails(tail_counter,4)=10;
            tailint=tailint+1;
        elseif(tailint>=44)
            tails(tail_counter,4)=9;
            tailint=tailint+1;
        elseif(tailint>=39)
            tails(tail_counter,4)=8;
            tailint=tailint+1;
        elseif(tailint>=34)
            tails(tail_counter,4)=7;
            tailint=tailint+1;
        elseif(tailint>=30)
            tails(tail_counter,4)=6;
            tailint=tailint+1;
        elseif(tailint>=26)
            tails(tail_counter,4)=5;
            tailint=tailint+1;
        elseif(tailint>=21)
            tails(tail_counter,4)=4;
            tailint=tailint+1;
        elseif(tailint>=16)
            tails(tail_counter,4)=3;
            tailint=tailint+1;
        elseif(tailint>=8)
            tails(tail_counter,4)=2;
            tailint=tailint+1;
        else
            tails(tail_counter,4)=1;
            tailint=tailint+1;
        end
   
   
   
   
   end   
   
%put core and linker and tail positions into a single array       
position_array=zeros(cores,3);
dna_counter=0;
i_counter=0;

for i=1:cores
    i_counter=i_counter+1;
    position_array(i_counter,1:3)=core_positions(i,1:3);
end

  
%calc geometric center
    rcm = mean(position_array(1:end_bead,1:3));
    
%calc Radius of gyration
    sumCM = 0;
    for j=1:end_bead
        v1 = position_array(j,:);
        v2 = rcm;
        v12 = v2-v1;
        r12sq = dot(v12,v12);
        sumCM = sumCM + r12sq;
    
    end
    Rg2cm = sumCM/(end_bead);
    Rgcm = sqrt(Rg2cm);
    %%% store
    Rgcm_list(zzz) = Rgcm;

    % Gyration tensor
    %%% Positions respect CM
    for j=1:end_bead
        frcores_xyzcm(j,:) = position_array(j,:) - rcm;
    end
    %%% Tensor
    vecrx = frcores_xyzcm(1:end_bead,1);
    vecry = frcores_xyzcm(1:end_bead,2);
    vecrz = frcores_xyzcm(1:end_bead,3);
    Sxx = dot(vecrx,vecrx)/end_bead;
    Syy = dot(vecry,vecry)/end_bead;
    Szz = dot(vecrz,vecrz)/end_bead;
    Sxy = dot(vecrx,vecry)/end_bead;
    Sxz = dot(vecrx,vecrz)/end_bead;
    Syz = dot(vecry,vecrz)/end_bead;
    Smn = [Sxx Sxy Sxz; Sxy Syy Syz; Sxz Syz Szz];   
    %%% Eigenvalues
    veceig = eig(Smn);
    %%% store
    eig_list(zzz,:) = veceig;
    
   clear k 
    %Use alphashape to get a more accurate volume
     k = alphaShape(position_array(1:end_bead, 1), position_array(1:end_bead, 2), position_array(1:end_bead,3),100);  
    volume_alphashape_list(zzz)=volume(k);
    
end

%yy_counter=yy_counter+1;
% gyration tensor parameters (output)
%%% Compute lists
vecx = eig_list(:,1);
vecy = eig_list(:,2);
vecz = eig_list(:,3);
Rgcmtens_list = sqrt(vecx + vecy + vecz);
asph_list = vecz - (vecx + vecy)/2;
acyl_list = vecy - vecx;
b2 = asph_list.*asph_list;
c2 = acyl_list.*acyl_list;
Rg4 = Rgcmtens_list.^4;
anis_list = (b2 + c2*(3/4))./Rg4;
volcell_list = sqrt(vecx.*vecy.*vecz);

Rgave(yy_counter) = mean(Rgcm_list(1:number_of_frames));
Rgstd(yy_counter) = std(Rgcm_list(1:number_of_frames));

volcellave(yy_counter) = mean(volcell_list(1:number_of_frames));
volcellstd(yy_counter) = std(volcell_list(1:number_of_frames));
volume_alphashape(yy_counter) = mean(volume_alphashape_list(1:number_of_frames));
volume_alphashape_std(yy_counter) = std(volume_alphashape_list(1:number_of_frames));


end

[volume_slope,intercept,MSE_volume]=logfit(bp_index(1:yy_counter),volume_alphashape(1:yy_counter),'loglog')
[rgyr_slope,intercept,MSE_Rgyr]=logfit(bp_index(1:yy_counter),Rgave(1:yy_counter),'loglog')

