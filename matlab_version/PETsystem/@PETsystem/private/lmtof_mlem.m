function x = lmtof_mlem(scanner, ...
                        lmdata, ...
                        image_size, ...
                        voxel_size, ...
                        tofinfo, ...
                        sensitivity, ...
                        mul_fac, ...
                        add_fac, ...
                        num_of_os, ...
                        nite, ...
                        xinit)
%function x = lmtof_mlem(scanner, 
%                        lmdata, 
%                        image_size, 
%                        voxel_size, 
%                        tofinfo,
%                        sensitivity, 
%                        mul_fac, 
%                        add_fac, 
%                        num_of_os, 
%                        nite, 
%                        xinit)
%
% Parameters:
%   scanner: a structure containing all geometry information about the targetted scanner
%   lmdata: listmode data [5 x num_of_events]
%   image_size: image size [i, j, k]
%   voxel_size: voxel size [i, j, k] (unit in mm)
%   tofinfo: time-of-flight information: [tw_resolution, tw_spacing] (unit in ps)
%   sensitivity: sensitivity (same size as an image) (default: 1)
%   mul_fac: multiplicative factors (default: 1)
%   add_fac: additive factors (default: 0)
%   num_of_os: number of subsets (maybe temporal subsets) (default: 1)
%   nite: number of iterations (default: 1)
%   xinit: initial image (default: 1)
%
%
% jnzhou@03-11-2013
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
    error('must specify TOF information!');
else
    if isempty(tofinfo)
        error('must specify TOF information!');
    end
end

if nargin < 6
    disp('not specify sensitivity, use default value!');
    sensitivity = ones(image_size);
else
    if isempty(sensitivity)
        disp('not specify sensitivity, use default value!');
        sensitivity = ones(image_size);
    end
end

if nargin < 7
    disp('not specify multipicative term, use default value!');
    mul_fac = ones(num_of_events, 1);
else
    if isempty(mul_fac) 
        disp('not specify multipicative term, use default value!');
        mul_fac = ones(num_of_events, 1);
    end
end

if nargin < 8
    disp('not specify additive term, use default value!');
    add_fac = zeros(num_of_events, 1);
else
    if isempty(add_fac) 
        disp('not specify additive term, use default value!');
        add_fac = zeros(num_of_events, 1);
    end
end

if nargin < 9
    disp('not specify number of subsets, use default value!');
    num_of_os = 1;
end

if nargin < 10
    disp('not specify number of iterations, use default value!');
    nite = 1;
end

if nargin < 11
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

fprintf('\nready to go ...\n\n');
                                         
%
for n = 1 : nite
    
    fprintf('***** iteration #%d *****\n', n);
    % 
    %
    
    for m = 1 : num_of_os

        fprintf('\nprocessing subset #%d (%d, %d) ...\n', m, num_of_os, n);
        
        % projectors
        fp = @(c)fproj_tof(c, image_size, voxel_size, ...
                           scanner.xtal_trans_location, ...
                           scanner.xtal_ring_offsets, ...
                           lmdata(:,m:num_of_os:end), tofinfo);
    
        bp = @(y)bproj_tof(y, image_size, voxel_size, ...
                           scanner.xtal_trans_location, ...
                           scanner.xtal_ring_offsets, ...
                           lmdata(:,m:num_of_os:end), tofinfo);        
        
        tic;
        yp = (fp(x) + add_fac(m:num_of_os:end)) ./ mul_fac(m:num_of_os:end);
    
        ye = 1 ./ yp; ye(yp == 0) = 0;

        xe = bp(ye); xe = reshape(xe, image_size);

        x = (x .* xe) ./ (sensitivity / num_of_os);
        x(sensitivity == 0) = 0;
        toc;
        %
        lk = -sum(sensitivity(:).*x(:)) + sum(log(yp(:) + 1e-30));
        fprintf('->likelihood = %.20f ...\n', lk);
        %
        
    end
    
end



