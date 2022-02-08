function [x, sbmat] = histo_recon(scanner, histo, image_size, voxel_size, sensitivity, num_of_os, num_of_ite, xinit)

% histo must be: data x maxstep

maxstep = size(histo, 2);
angle_list = (0:(maxstep-1))/maxstep * 180;
nt = scanner.nxtal_trans;
na = scanner.nxtal_axial;
nb = scanner.num_of_doi_bins;

%
[t1,t2,d1,d2,r1,r2] = ndgrid([0:(nt-1)],[nt:(2*nt-1)],[0:(nb-1)],[0:(nb-1)],[0:(na-1)],[0:(na-1)]);
lm = [t1(:),t2(:),d1(:),d2(:),r1(:),r2(:)]';
lm = reshape(lm, 6, nt*nt*nb*nb, na*na);
lm = nat2span(lm, sino_order_map(na));
lm = reshape(lm, 6, nt*nt*nb*nb*na*na);
lm = lm([1 5 4 2 6 3], :);
lm(7,:) = 0;
lm = uint8(lm);

%
fov_mask=circmask3d(image_size, image_size(1)/2.0-0.5, ...
					image_size(2)/2.0-0.5, image_size(1)/2-2.5);

%
blur_enabled = 1;
width=5;
fw=1.0;
f=fspecial('gaussian', width, fw/2.355);
if blur_enabled
	Rt = psfmat(f, image_size(1:2));
	Ra = speye(image_size(3));
	R = kron(Ra, Rt);
	fprintf('blur enabled: width=%d, fwhm=%.3f\n', width, fw);
end

%


if isempty(sensitivity)
	s_enabled = 0;
	disp 'not support sensitivity, calc on the fly!';
else
	s_enabled = 1;
end

%
if isempty(xinit)
	x = ones(image_size) .* fov_mask;
else
	x = xinit;
end

fprintf('ready, set, go! \n');

for n = 1 : num_of_ite

	fprintf('*** iteration #%d ***\n', n);

	if ~s_enabled
		if n==1
			sbmat = zeros([image_size, num_of_os]);
		end
	end

	for i = 1 : num_of_os
	
		ang = i:num_of_os:maxstep;
		fprintf('--- subset #%d ---\n', i);
		tic;
		
		% fp
		xb = zeros(image_size);
		if blur_enabled
%			x=imfilter(x, f);
			x=reshape(R*x(:), image_size);
		end
		x = x .* fov_mask;
		for j = 1 : length(ang)
			
			fprintf('processing %d in subset #%d(%d) at ite.#%d(%d) ...\n', ...
					j, i, num_of_os, n, num_of_ite)
			
			lm(7,:) = ang(j)-1;
			y = fproj_pawpet(x, image_size, voxel_size, ...
							 scanner.xtal_trans_location, ...
							 scanner.xtal_ring_offsets, lm, angle_list);
				
			yp = histo(:,ang(j)) ./ (y + 0.001); %yp(y==0) = 0;

			lm(7,:) = ang(j)-1;
			bp = bproj_pawpet(yp, image_size, voxel_size, ...
							  scanner.xtal_trans_location, ...
							  scanner.xtal_ring_offsets, lm, angle_list);
			if ~s_enabled				  
				if n==1				  
					lm(7,:) = ang(j)-1;
					s0 = bproj_pawpet([], image_size, voxel_size, ...
								  	  scanner.xtal_trans_location, ...
								  	  scanner.xtal_ring_offsets, lm, angle_list);
					sbmat(:,:,:,i) = sbmat(:,:,:,i) + reshape(s0, image_size);
				end
			end
			%
			xb(:) = xb(:) + bp(:);				  
		end
		if blur_enabled
%			xb = imfilter(xb, f);
			xb = reshape(R'*xb(:), image_size);
			if ~s_enabled
				if n==1
%					sbmat(:,:,:,i) = imfilter(sbmat(:,:,:,i), f);
					sbmat(:,:,:,i) = reshape(R'*col(sbmat(:,:,:,i)), image_size);
				end
			end
		end
		xb = xb .* fov_mask;
		
		% image update
		if ~s_enabled
			x = x ./ sbmat(:,:,:,i) .* xb; 
			x(sbmat(:,:,:,i) == 0) = 0;	
		else
			x = x ./ sensitivity(:,:,:,i) .* xb; 
			x(sensitivity(:,:,:,i) == 0) = 0;	
		end
		
		toc;
				
	end
	
	% save
	if mod(n, 200)==0
		save(sprintf('x%d.mat', n), 'x');
	end	
	
end

