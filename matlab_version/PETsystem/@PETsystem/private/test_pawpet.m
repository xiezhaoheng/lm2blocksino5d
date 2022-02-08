%
scanner_name = 'pawPET';

% inveon 3d
%
switch scanner_name
    case 'pawPET'
        xtal_attn = 0.087;
        angle_offset = 0.0;
        FOV_diameter = 25; %mm
        num_of_doi_bins=4;
        scanner = build_scanner(30.0, ...
                                [0.5, 0.5, 0.5, 8], ...
                                [0, 0, 0], ...
                                32, 32, 2, 1, ...
                                num_of_doi_bins, angle_offset, 1);
    % don't use scanner.xtal_pairs!
end;

%
image_size = [100 100 32];
voxel_size = [0.15, 0.15, 0.5];
maxstep = 40;
angle_list = (0:(maxstep-1))/maxstep * 180.0;
nt = scanner.nxtal_trans;
na = scanner.nxtal_axial;

return;

%calc sensitivity
xss = zeros([image_size, 11]);
for n=1:11

	fprintf('processing %d ...\n', n);
	if n==1
		clear ll;
		[i1,i2,r1,r2,d1,d2]=ndgrid(0:(nt-1), nt:(nt*2-1), 0:(na-1),0:(na-1), ...
								   0:(num_of_doi_bins-1),0:(num_of_doi_bins-1));
		ll(1,:)=i1(:);
		ll(2,:)=r1(:);
		ll(3,:)=d1(:);
		ll(4,:)=i2(:);
		ll(5,:)=r2(:);
		ll(6,:)=d2(:);
		ll(7,:)=ones(size(ll(1,:))) * (n-1);
	else
		ll(7,:)=n-1;
	end

	if 1
		tic;
		x = bproj_pawpet_rect([], image_size, voxel_size, ...
			                 scanner.xtal_trans_location, scanner.xtal_ring_offsets, ...
			                 uint8(ll(:,:)), angle_list);
		toc;

		xss(:,:,:,n) = reshape(x, image_size);
	end

end

% form a complete sensitivity image
if 1
	s0 = xss(:,:,:,1);
	s2 = sum(xss(:,:,:,[2:end-1]),4);
	s3 = xss(:,:,:,end);
	
	senimg = s0 + s2 + s3 + imrotate(s2, 90) + imrotate(s0,90) + ...
			 imrotate(s3, 90) + flipdim(s2,2) + flipdim(imrotate(s2,90),2);
			 
	
	
end
