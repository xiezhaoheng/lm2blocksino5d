%
scanner_name = 'taperedPET';

% inveon 3d
%
switch scanner_name
    case 'taperedPET'
        xtal_attn = 0.087;
        angle_offset = 0.0;
        FOV_diameter = 30; %mm
        scanner = build_scanner(61.0, ...
                                [0.5, 0.865, 0.5, 13], ...
                                [0, 0, 0], ...
                                14, 14, 16, 1, ...
                                15, angle_offset, 71);
end;

%
image_size = [101 101 14];
voxel_size = [0.2, 0.2, 0.5];

xp=scanner.xtal_pairs;

%
ll = zeros(6, size(xp,2));
ll(1,:) = xp(1,:)-1;
ll(4,:) = xp(2,:)-1;

xb=0;
for n=1:15
    for m=1:15
    
        ll(3,:)=m-1;
        ll(6,:)=n-1;
        
        tic;
        x = bproj_doi([], image_size, voxel_size, ...
                      scanner.xtal_trans_location, scanner.xtal_ring_offsets, ...
                      uint16(ll(:,:)));
        toc;

        xb = xb + x;
    end
end
%

if 0

%
senimg = zeros(image_size);
for n = 1 : nring
    
    fprintf('processing #%d ...\n', n);
    x = bpmat(:,:,:,n); 

    for m = 1 : (nring-n+1) 
        
        st = (m-1)*axial_ratio;
        x_shift = zeros(image_size);           
        x_shift(:, :, (st+1):end) = x(:, :, 1:(end-st));
            
        senimg = senimg + x_shift;

        if n > 1
            senimg = senimg + flipdim(x_shift,3);
        end
            
    end        
        

end

end