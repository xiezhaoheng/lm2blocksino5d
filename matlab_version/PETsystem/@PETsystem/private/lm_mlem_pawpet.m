function x = lm_mlem_pawpet(scanner, ...
                        	lmdata, ...
                        	image_size, ...
                        	voxel_size, ...
                        	maxstep, ...
                        	sensitivity, ...
                        	num_of_os, ...
                        	nite, ...
                        	xinit)

%function x = lm_mlem_pawpet(scanner, 
%                        	lmdata, 
%                        	image_size, 
%                        	voxel_size, 
%                        	maxstep,
%                        	sensitivity, 
%                        	num_of_os, 
%                        	nite, 
%                        	xinit)
%
% Parameters:
%   scanner: a structure containing all geometry information about the targetted scanner
%   lmdata: listmode data [7 x num_of_events]
%   image_size: image size [i, j, k]
%   voxel_size: voxel size [i, j, k] (unit in mm)
%   maxstep: maximum number of rotation angles
%   sensitivity: sensitivity (same size as an image) (default: 1)
%   num_of_os: number of subsets (maybe temporal subsets) (default: 1)
%   nite: number of iterations (default: 1)
%   xinit: initial image (default: 1)
%
%
% jnzhou@04-18-2013
%

help(mfilename('fullname'));

if nargin < 1
    error('create your scanner first!');
end

if nargin < 2
    error('must specify listmode data!');    
else
    if isempty(lmdata)
        error('must specify listmode data!');
    end
end

num_of_events = size(lmdata,2);
fprintf('total # of counts: %d\n', num_of_events);

% check data

if any(lmdata(2,:)<0) | any(lmdata(2,:) >= scanner.nxtal_axial)
	error('invalid lmdata value:a!');
end

if any(lmdata(5,:)<0) | any(lmdata(5,:) >= scanner.nxtal_axial)
	error('invalid lmdata value:a!');
end

if any(lmdata(1,:)<0) | any(lmdata(1,:) >= scanner.nxtal_axial*2)
	error('invalid lmdata value:t!');
end

if any(lmdata(4,:)<0) | any(lmdata(4,:) >= scanner.nxtal_axial*2)
	error('invalid lmdata value:t!');
end

if any(lmdata(3,:)<0) | any(lmdata(3,:) >= scanner.num_of_doi_bins)
	error('invalid lmdata value:doi!');
end

if any(lmdata(6,:)<0) | any(lmdata(6,:) >= scanner.num_of_doi_bins)
	error('invalid lmdata value:doi!');
end

if nargin < 3
    error('must specify image size!');
else
    if isempty(image_size)
        error('must specify image size!');
    end
end

if nargin < 4
    error('must specify voxel size!');
else
    if isempty(voxel_size)
        error('must specify voxel size!');
    end
end

if nargin < 5
	error('must specify maximum number of rotation steps!');
else
	if isempty(maxstep) 
		disp('not specify maximum number of rotation steps, use default value');
		maxstep = 40;
	end
end

angle_list = (0:(maxstep-1))/maxstep * 180.0;

if nargin < 6
    disp('not specify sensitivity, use default value!');
    sensitivity = ones(image_size);
else
    if isempty(sensitivity)
        disp('not specify sensitivity, use default value!');
        sensitivity = ones(image_size);
    end
end

%
fw=2.0; width=5; blur_enabled=1;
if blur_enabled
	sensitivity=smooth3(sensitivity, 'gaussian', width, fw/2.355);
%	f=fspecial('gaussian', 9, fw/2.355);
%	sensitivity = imfilter(sensitivity, f);
end
%

if nargin < 7
    disp('not specify number of subsets, use default value!');
    num_of_os = 1;
end

if nargin < 8
    disp('not specify number of iterations, use default value!');
    nite = 1;
end

if nargin < 9
    disp('not specify initial image, use default value!');
    x = ones(image_size);
else
    if isempty(xinit)
        disp('not specify initial image, use default value!');
        x = ones(image_size);
    else
        x = xinit;
    end
end

fprintf('\nready, set, go! ...\n\n');

%
for n = 1 : nite
    
    fprintf('***** iteration #%d *****\n', n);
    % 
    %
    
    for m = 1 : num_of_os

        fprintf('\nprocessing subset #%d (%d, %d[%d]) ...\n', m, num_of_os, n, nite);
        
        % projectors
        fp = @(c)fproj_pawpet(c, image_size, voxel_size, ...
                           scanner.xtal_trans_location, ...
                           scanner.xtal_ring_offsets, ...
                           uint8(lmdata(:,m:num_of_os:end)), angle_list);
    
        bp = @(y)bproj_pawpet(y, image_size, voxel_size, ...
                           	  scanner.xtal_trans_location, ...
                           	  scanner.xtal_ring_offsets, ...
                           	  uint8(lmdata(:,m:num_of_os:end)), angle_list);        
        
        tic;
        %
        if blur_enabled
        	x=smooth3(x,'gaussian', width, fw/2.355);
%        	x = imfilter(x, f);
        end
        %
        yp = fp(x);
    
        ye = 1 ./ (yp + eps); %ye(yp == 0) = 0;

        xe = bp(ye); xe = reshape(xe, image_size);

		if blur_enabled
			xe=smooth3(xe,'gaussian', width, fw/2.355);
%			xe = imfilter(xe, f);
		end


        x = (x .* xe) ./ (sensitivity / num_of_os);
        x(sensitivity == 0) = 0;
        toc;
        %
        lk = -sum(sensitivity(:).*x(:)) + sum(log(yp(:) + 1e-30));
        fprintf('->likelihood = %.20f ...\n', lk);
        %
        
    end
    
end
                        
