%This script give you the same result as C code(lm2blocksino5d_v2)
%Notice: need to change sign of TOF bin when we swap the eventA and eventB
% Any question, please contact Zhaoheng(zhxie@ucdavis.edu)
clear
clc

p=genpath('./PETsystem');  %
addpath(p);
fname_in='./f00000.iter1_trues.lm'; %input listmode-file directory
%%  for TOF sino rebin
scanner_name = 'explorer2000mm_unitedimaging_bxp_120x77_24x5';
scanner = buildPET(scanner_name);
xp = scanner.getDefaultSinogramCrystalPairs;
tc = scanner.getCrystalTransaxialLocations; % tc first row is x axial, second row is y axial
np = scanner.system_parms.number_of_projections_per_angle;
xp_plot = reshape(xp, 2, np, size(xp,2)/np);
%figure,scanner.plot("2d,geom")

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

lm_lor_sino= zeros( numxp2d,numrings,numrings,num_tof_bin, 'single');

blk_idx=zeros(N_TX_BLK,N_TX_BLK);

for i=1:numxp2d
    blk_idx(xp(1,i),xp(2,i))=i;
end

%%  Read Listmode
%figure,imshow_zj(blk_idx,'nc'),title('X and Y axis are transaxial crystal ID, entry is sinogram ID');

lm_cry_id = fread(fopen(fname_in, 'rb'), inf, 'int16=>int16');
fclose('all');
length_num = length(lm_cry_id)/5;
lm_cry_id = reshape(lm_cry_id, 5, length_num);

length_num = size(lm_cry_id,2);
lm_cry_id=single(lm_cry_id);
lm_blk_id =single(zeros(5,length_num));
lm_blk_id_temp =lm_blk_id;
%% 

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
    
    if  (temp_tof>0) || (temp_tof==0)
        lm_blk_id(5,i) = fix(temp_tof/N_TIME_BIN_PER_TOF);
    else
        lm_blk_id(5,i) = fix((temp_tof-N_TIME_BIN_PER_TOF)/N_TIME_BIN_PER_TOF);
    end
    
    %% swap the BLK list-mode

    lm_blk_id_temp(:,i) =lm_blk_id(:,i) ;
    

    if  (lm_blk_id_temp(5,i)<14 )&& (lm_blk_id_temp(5,i)>-14 )&&(blk_idx(lm_blk_id_temp(1,i),lm_blk_id_temp(3,i))~=0)
        
        lm_lor_sino(blk_idx(lm_blk_id_temp(1,i),lm_blk_id_temp(3,i)),lm_blk_id_temp(2,i),lm_blk_id_temp(4,i),lm_blk_id_temp(5,i)+14)=...
            lm_lor_sino(blk_idx(lm_blk_id_temp(1,i),lm_blk_id_temp(3,i)),lm_blk_id_temp(2,i),lm_blk_id_temp(4,i),lm_blk_id_temp(5,i)+14)+1;
        
    elseif (lm_blk_id_temp(5,i)<14 )&& (lm_blk_id_temp(5,i)>-14 )&&(blk_idx(lm_blk_id_temp(3,i),lm_blk_id_temp(1,i))~=0)  %swap the orientation of the LOR
        lm_lor_sino(blk_idx(lm_blk_id_temp(3,i),lm_blk_id_temp(1,i)),lm_blk_id_temp(2,i),lm_blk_id_temp(4,i),-lm_blk_id_temp(5,i)+14)=...
            lm_lor_sino(blk_idx(lm_blk_id_temp(3,i),lm_blk_id_temp(1,i)),lm_blk_id_temp(2,i),lm_blk_id_temp(4,i),-lm_blk_id_temp(5,i)+14)+1;
    else

        continue
        
    end
    
end

lm_lor_sino=reshape(lm_lor_sino, [numrad numang numrings numrings num_tof_bin]);
data_out = sum(lm_lor_sino,4);
data_out = sum(data_out,3);
data_out=squeeze(data_out);

figure,subplot(131),imshow_zj(data_out(:,:,12),'nc'),title('SIMSET:TOF bin# 12'),
hold on,subplot(132),imshow_zj(data_out(:,:,14),'nc'),title('SIMSET:TOF bin# 14'),
hold on,subplot(133),imshow_zj(data_out(:,:,16),'nc'),title('SIMSET:TOF bin# 16'),
data_4d=sum(lm_lor_sino,5);
figure,imshow_zj(squeeze(sum(sum(data_4d,1),2)),'nc'),title('SIMSETsum radial and angular');

fname = './lm2sino.5d';
fwrite(fopen(fname, 'w'), lm_lor_sino, 'uint64');
fclose('all');
