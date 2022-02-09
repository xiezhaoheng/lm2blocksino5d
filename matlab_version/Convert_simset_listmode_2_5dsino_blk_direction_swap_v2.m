% switch to change orientation in BLK level instead of crystal level
% using getDefaultSinogramCrystalPairs_w_orientation instaed of getDefaultSinogramCrystalPairs to extract Cartesian coordinate of BLK or Crystal, remember to change the orientation of tc, simset output has 
%{ 
% Zhoujian's original defination
 * ^ y+
 *  |
 *  D   |   C
 *  +---|---+
 *  |   |   |
 *  ----|------> x+
 *  |   |   |
 *  +---|---+
 *  A   |   B

% Simset output defination,need to swap x-y and reverse the x axis
 * ^ x-
 *  |
 *  D   |   C
 *  +---|---+
 *  |   |   |
 *  ----|------> y+
 *  |   |   |
 *  +---|---+
 *  A   |   B

%}


clear
clc
p=genpath('/home/raid7/zhxie/Explorer_phantom/2021-07-01/5min_sino_orientation_phantom/PETsystem');  %in raid 7 /home/raid7/zhxie/Explorer_phantom/2021-07-01/5min_sino_orientation_phantom/PETsystem
addpath(p);
scanner_name = 'explorer2000mm_unitedimaging_bxp_120x77_24x5';
scanner = buildPET(scanner_name);
%xp = scanner.getDefaultSinogramCrystalPairs;
xp = scanner.getDefaultSinogramCrystalPairs_w_orientation; % Using revised version to order the listmode
tc = scanner.getCrystalTransaxialLocations; % tc first row is x axial, second row is y axial
np = scanner.system_parms.number_of_projections_per_angle;
xp_plot = reshape(xp, 2, np, size(xp,2)/np);
figure,scanner.plot("2d,geom")

hold on;
%% In BLK level

numrad = scanner.system_parms.number_of_projections_per_angle;
numang = scanner.getDefaultNumberOfAngles;
numxp2d = size(xp, 2);
numrings = scanner.getNumberOfCrystalRings;
num_tof_bin=27;
N_TX_BLK=2*numang;
N_TX_CRYS_PER_BLK= 7;
N_AX_CRYS_PER_BLK=12;
N_TIME_BIN_PER_TOF=7;
LM_add_Gap=true;
N_AX_Unit= 84;
n_ax_bnk=56;

%lm_lor_sino_non_tof= zeros( numxp2d,numrings,numrings, 'single');
lm_lor_sino= zeros( numxp2d,numrings,numrings,num_tof_bin, 'single');

blk_idx=zeros(N_TX_BLK,N_TX_BLK);

for i=1:numxp2d
        blk_idx(xp(1,i),xp(2,i))=i;
end


%figure,imshow_zj(blk_idx,'nc'),title('X and Y axis are transaxial crystal ID, entry is sinogram ID');
fname_in='/home/raid7/zhxie/Explorer_phantom/2021-07-01/5min_sino_orientation_phantom/Simset_generate_sino/new_water_material_1e8count_erode7e2_positron_false/f00000.iter1_trues.lm';
lm_cry_id = fread(fopen(fname_in, 'rb'), inf, 'int16=>int16');
fclose('all');
length_num = length(lm_cry_id)/5;
lm_cry_id = reshape(lm_cry_id, 5, length_num);

%lm_cry_id=lm_cry_id(:,1:50:end); %accelerate  debug process

length_num = size(lm_cry_id,2);
lm_cry_id=single(lm_cry_id);
lm_blk_id =single(zeros(5,length_num)); 
lm_blk_id_temp =lm_blk_id;


for i=1:length_num
    lm_blk_id(1,i) =fix( lm_cry_id(1,i)/N_TX_CRYS_PER_BLK)+1; %index start from 1
    lm_blk_id(3,i) = fix( lm_cry_id(3,i)/N_TX_CRYS_PER_BLK)+1; %index start from 1
    if LM_add_Gap==false
        lm_blk_id(2,i) = fix(lm_cry_id(2,i)/N_AX_CRYS_PER_BLK)+1;
        lm_blk_id(4,i) =  fix(lm_cry_id(4,i)/N_AX_CRYS_PER_BLK)+1;
    else %already add gap
        temp_axial1= lm_cry_id(2,i)-fix( lm_cry_id(2,i)/(N_AX_Unit+1));
        temp_axial2= lm_cry_id(4,i)-fix( lm_cry_id(4,i)/(N_AX_Unit+1));
        lm_blk_id(2,i) = fix(temp_axial1/N_AX_CRYS_PER_BLK)+1;
        lm_blk_id(4,i) =  fix(temp_axial2/N_AX_CRYS_PER_BLK)+1;
    end
     temp_tof=lm_cry_id(5,i)+3;
    %temp_tof=lm_cry_id(5,i);
    
    if  (temp_tof>0) || (temp_tof==0)
        lm_blk_id(5,i) = fix(temp_tof/N_TIME_BIN_PER_TOF);
    else
        lm_blk_id(5,i) = fix((temp_tof-N_TIME_BIN_PER_TOF)/N_TIME_BIN_PER_TOF);
    end
    
    %% switch the BLK list-mode
    event1_trans=tc(:,lm_blk_id(1,i));
    event2_trans=tc(:,lm_blk_id(3,i));
    
    event1_trans=flip(event1_trans);event1_trans(2)=-event1_trans(2);
    event2_trans=flip(event2_trans);event2_trans(2)=-event2_trans(2);
     
    lm_blk_id_temp(:,i) =lm_blk_id(:,i) ;

    if event1_trans(1)<event2_trans(1) % x_a<x_b, swap A and B
        lm_blk_id_temp(1,i) = lm_blk_id(3,i);
        lm_blk_id_temp(2,i) = lm_blk_id(4,i);
        lm_blk_id_temp(3,i) = lm_blk_id(1,i);
        lm_blk_id_temp(4,i) = lm_blk_id(2,i);
        lm_blk_id_temp(5,i) =- lm_blk_id(5,i);
    elseif round(event1_trans(1),4)==round(event2_trans(1),4)              
        if   event1_trans(2)<event2_trans(2) % y_a<y_b, swap A and B
            lm_blk_id_temp(1,i) = lm_blk_id(3,i);
            lm_blk_id_temp(2,i) = lm_blk_id(4,i);
            lm_blk_id_temp(3,i) = lm_blk_id(1,i);
            lm_blk_id_temp(4,i) = lm_blk_id(2,i);
            lm_blk_id_temp(5,i) = -lm_blk_id(5,i);
        
        end
    else
        continue
    end
    
    if  (lm_blk_id_temp(5,i)<14 )&& (lm_blk_id_temp(5,i)>-14 )&&(blk_idx(lm_blk_id_temp(1,i),lm_blk_id_temp(3,i))~=0)
        lm_lor_sino(blk_idx(lm_blk_id_temp(1,i),lm_blk_id_temp(3,i)),lm_blk_id_temp(2,i),lm_blk_id_temp(4,i),lm_blk_id_temp(5,i)+14)=lm_lor_sino(blk_idx(lm_blk_id_temp(1,i),lm_blk_id_temp(3,i)),lm_blk_id_temp(2,i),lm_blk_id_temp(4,i),lm_blk_id_temp(5,i)+14)+1;
    else

        continue
        
    end
    

end


lm_lor_sino=reshape(lm_lor_sino, [numrad numang numrings numrings num_tof_bin]);
lm_lor_sino=permute(lm_lor_sino,[1 2 4 3 5]); % From the sinogram, it look the axial dimension need to be swap
 data_out = sum(lm_lor_sino,4);
 data_out = sum(data_out,3);
 data_out=squeeze(data_out);

figure(3),subplot(131),imshow_zj(data_out(:,:,12),'nc'),title('SIMSET:TOF bin# 12'),
figure(3),subplot(132),imshow_zj(data_out(:,:,14),'nc'),title('SIMSET:TOF bin# 14'),
figure(3),subplot(133),imshow_zj(data_out(:,:,16),'nc'),title('SIMSET:TOF bin# 16'),
data_4d=sum(lm_lor_sino,5);
figure,imshow_zj(squeeze(sum(sum(data_4d,1),2)),'nc'),title('SIMSETsum radial and angular');


fname = './lm2sino_change_direction.5d';
fwrite(fopen(fname, 'w'), lm_lor_sino, 'uint64');
fclose('all');



