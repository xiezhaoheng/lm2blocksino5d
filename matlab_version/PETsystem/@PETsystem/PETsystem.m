%A class to perform fast listmode- or sinogram-based PET image reconstruction
%
%
%Author: Jian Zhou
%Date: Aug. 1, 2013
%
% Add orientation criteria, 2021 Zhaoheng Xie
classdef PETsystem
    
    properties(GetAccess = 'public', SetAccess = 'private')
        name_tag; % scanner name
        system_parms; % structure that contains various system parameters
    end
    
    methods(Access = 'public')
        
        function obj = PETsystem(...
                tag, ...
                ring_diameter, ...
                crystal_size, ...
                crystal_gap_size, ...
                crystal_array_size, ...
                number_of_detector_modules_transaxial, ...
                number_of_detector_modules_axial, ...
                number_of_DOI_bins, ...
                detector_module_initial_angle_offset, ...
                detector_modula_axial_extra_offsets, ...
                number_of_projections_per_angle, ...
                tof_info)
            
            obj.name_tag = tag;
            obj.system_parms.ring_diameter = ring_diameter;            
            obj.system_parms.crystal_size = crystal_size;
            obj.system_parms.crystal_gap_size = crystal_gap_size;
            obj.system_parms.crystal_array_size = crystal_array_size;
            obj.system_parms.number_of_detector_modules_transaxial = ...
                number_of_detector_modules_transaxial;
            obj.system_parms.number_of_detector_modules_axial = ...
                number_of_detector_modules_axial;
            obj.system_parms.number_of_DOI_bins = number_of_DOI_bins;
            obj.system_parms.detector_module_initial_angle_offset = ...
                detector_module_initial_angle_offset;
            obj.system_parms.detector_module_axial_extra_offsets = ...
                detector_modula_axial_extra_offsets;
            obj.system_parms.number_of_projections_per_angle = ...
                number_of_projections_per_angle; 
            obj.system_parms.tof_info = tof_info;
            obj.system_parms.projector = 'Siddon';
            obj.system_parms.depth_ratio = 0.5;
        end
       
        function obj = setNumberOfDetectorModules(obj, num_trans, num_axial)
        	obj.system_parms.number_of_detector_modules_transaxial = num_trans;
        	obj.system_parms.number_of_detector_modules_axial = num_axial;
        end
       
		function obj = setCrystalArraySize(obj, crystal_array_size)
			if length(crystal_array_size) ~= 2
				error('must be: [t a]!');
			end
			obj.system_parms.crystal_array_size = crystal_array_size;
		end
       
        function obj = setCrystalSize(obj, crystal_size)
        	if length(crystal_size) ~= 4
        		error('must be 4 elements:[tf tr a d]!');
        	end
        	obj.system_parms.crystal_size = crystal_size;
        end
        
        function obj = setGapSize(obj, gap_size)
        	if length(gap_size) ~= 3
        		error('must be 3 elements: [tf tr a]!');
        	end
        	obj.system_parms.crystal_gap_size = gap_size;
        end
       
        function obj = setRingDiameter(obj, ring_diameter)
            obj.system_parms.ring_diameter = ring_diameter;
        end

        function obj = setProjector(obj, projector)
        %set projector type: `Siddon' or `Bresenham' or 'Linterp'
        %function obj = setProjector(projector)
        	obj.system_parms.projector = projector;
        end
        
        function obj = setTOFInfo(obj, tof_info)
            obj.system_parms.tof_info = tof_info;
        end

        function obj = setDepthRatio(obj, ratio)
        %set a ratio to determine average depth interaction point inside crystal
        %0.5 is default, meaning at the crystal center
        %function obj = setDepthRatio(ratio)
        	if ratio < 0 || ratio > 1.0
        		error('ratio must be in [0 1]!');
        	end
        	obj.system_parms.depth_ratio = ratio;
        end
        
        function obj = setNumberOfProjectionsPerAngle(obj, num_of_radial_bins)
        %set number of radial bins
        %function obj = setNumberOfProjectionsPerAngle(num_of_radial_bins)
        	obj.system_parms.number_of_projections_per_angle = num_of_radial_bins;
        end
        
        function obj = setDetectorModuleInitialAngleOffset(obj, angle)
        %function obj = setDetectorModuleInitialAngleOffset(obj, angle)
        %set angle offset (deg)
        	obj.system_parms.detector_module_initial_angle_offset = angle;
        end
       
        function dmo = getDetectorModuleAxialOffsets(obj)
        %calculate detector module axial offsets
		%function dmo = getDetectorModuleAxialOffsets        
            crystal_axial_pitch = obj.system_parms.crystal_size(3) + ...
                                  obj.system_parms.crystal_gap_size(3);
            detector_module_axial_size = crystal_axial_pitch * ...
                                         obj.system_parms.crystal_array_size(2);
            
            nb = obj.system_parms.number_of_detector_modules_axial;
            dmo = (- nb * 0.5 + 0.5 + (0 : (nb-1))) * ...
                   detector_module_axial_size + ...
                   obj.system_parms.detector_module_axial_extra_offsets;            
        end
        
        function nr = getNumberOfCrystalRings(obj)
        %get number of crystal rings
        %function nr = getNumberOfCrystalRings
            nr = obj.system_parms.crystal_array_size(2) * ...
                obj.system_parms.number_of_detector_modules_axial;
        end
        
        function na = getDefaultNumberOfAngles(obj)
        %get number of projection angles in default mode
        %function na = getDefaultNumberOfAngles
            na = obj.system_parms.crystal_array_size(1) * ...
                obj.system_parms.number_of_detector_modules_transaxial / 2;
        end
        
        function ro = getCrystalRingOffsets(obj)   
        %get crystal ring axial offsets
        %function ro = getCrystalRingOffsets   
            crystal_axial_pitch = obj.system_parms.crystal_size(3) + ...
                                  obj.system_parms.crystal_gap_size(3);
            dm_offsets = getDetectorModuleAxialOffsets(obj);
            nb = obj.system_parms.crystal_array_size(2);
            xt_centers = (- nb * 0.5 + 0.5 + (0 : (nb-1))) * crystal_axial_pitch;
            ro = [];          
            for n = 1 : obj.system_parms.number_of_detector_modules_axial
                ro = [ro , xt_centers + dm_offsets(n)];
            end
        end
        
        function tc = getCrystalTransaxialLocations(obj)            
        %calculate crystal bin transaxial coordinates (bin means DOI bin)
        %function tc = getCrystalTransaxialLocations
            tfs = (obj.system_parms.crystal_size(1) + obj.system_parms.crystal_gap_size(1));
            trs = (obj.system_parms.crystal_size(2) + obj.system_parms.crystal_gap_size(2));
            nxtal_trans = obj.system_parms.crystal_array_size(1);
            df = (-nxtal_trans*0.5 + 0.5 + (0 : nxtal_trans-1)) * tfs;
            dr = (-nxtal_trans*0.5 + 0.5 + (0 : nxtal_trans-1)) * trs;
            
            num_of_doi_bins = obj.system_parms.number_of_DOI_bins;
            xtal_loc = zeros(2, num_of_doi_bins, nxtal_trans);
            xtal_size_depth = obj.system_parms.crystal_size(4);
            for n=1: nxtal_trans
                l = df(n)-dr(n);
                dl = l / num_of_doi_bins;
                cl = l * 0.5 - ((0:num_of_doi_bins-1) + 0.5) * dl;
                h = xtal_size_depth;
                dh = h / num_of_doi_bins;
                ch = h * 0.5 - ((0:num_of_doi_bins-1) + 0.5) * dh;
                xtal_loc(:, :, n) = [cl(:) + dr(n), ch(:)]';
            end
            
            R = obj.system_parms.ring_diameter * 0.5;
            nblock_trans = obj.system_parms.number_of_detector_modules_transaxial;
            tc = zeros(2, num_of_doi_bins, nxtal_trans, nblock_trans);
            for n=1:nblock_trans                
                % always start at 3 o'clock
                a0 = (n-1) * 2*pi / nblock_trans  + ...
                    obj.system_parms.detector_module_initial_angle_offset * pi / 180; 
                xs = xtal_loc(1,:,:);
                ys = xtal_loc(2,:,:);
    
                t = pi * 0.5 + a0;
                x = xs * cos(t) - ys * sin(t) + (R + xtal_size_depth*obj.system_parms.depth_ratio) * cos(a0);
                y = xs * sin(t) + ys * cos(t) + (R + xtal_size_depth*obj.system_parms.depth_ratio) * sin(a0);
    
                tc(1,:,:,n) = x;
                tc(2,:,:,n) = y;    
            end
            tc = reshape(tc, 2, num_of_doi_bins, nxtal_trans*nblock_trans);            
            tc = squeeze(tc);
        end
        
        function plot(obj, options)            
        % plot scanner geometry with options: 2d,trans,simple; 2d,trans,simple,+lor; 3d,simple
        %function plot(options)
            switch lower(options)
            	case '2d,geom'
            		line_color = [0, 0, 1];
            		c_width_front = obj.system_parms.crystal_size(1) + obj.system_parms.crystal_gap_size(1);
            		c_width_back = obj.system_parms.crystal_size(2) + obj.system_parms.crystal_gap_size(2);
            		b_width_front = obj.system_parms.crystal_array_size(1) * c_width_front;
            		b_width_back = obj.system_parms.crystal_array_size(1) * ...
            			 	(obj.system_parms.crystal_size(2) + obj.system_parms.crystal_gap_size(2));
            		b_depth = obj.system_parms.crystal_size(end);
            		for n = 1 : (obj.system_parms.crystal_array_size(1)+1)
            			pfront(n,:) = [-b_depth/2, b_width_front/2 - (n-1)*c_width_front];
            			pback(n,:) = [b_depth/2, b_width_back/2 - (n-1)*c_width_back];
            		end 
            		
            		nblock_trans = obj.system_parms.number_of_detector_modules_transaxial;        
					a0 = obj.system_parms.detector_module_initial_angle_offset;      
					cc = obj.system_parms.ring_diameter * 0.5 + b_depth / 2;
					plot(0, 0, '+'); hold on;
            		for m = 1 : nblock_trans
            		            		
            			a = 360 / nblock_trans * pi / 180 * (m - 1) + a0 * pi / 180;        		
 						r = [cos(a), -sin(a); sin(a), cos(a)];					
 						t = [cc * cos(a); cc * sin(a)];
 					           		            		
	            		pf = r * pfront' + repmat(t, 1, size(pfront,1));
	            		pb = r * pback' + repmat(t, 1, size(pback, 1));
	            		for n = [1:(obj.system_parms.crystal_array_size(1)+1)]
	            			line([pf(1,n), pb(1,n)], [pf(2,n), pb(2,n)], 'color', line_color); 
	            		end
	            		line([pf(1,1), pf(1,end)], [pf(2,1), pf(2,end)], 'color', line_color);
	            		line([pb(1,1), pb(1,end)], [pb(2,1), pb(2,end)], 'color', line_color);
	            		
            		end
                	maxc = ceil((cc+b_depth)/100)*100;
            		line([-maxc,maxc], [0, 0], 'linestyle', '--', 'color', 'k');
            		line([0, 0], [-maxc,maxc], 'linestyle', '--', 'color', 'k');
            		hold off;
            		axis equal;
					set(gca, 'drawMode', 'fast', ...
                		'xtickmode', 'manual', ...
                		'ytickmode', 'manual');
            		xlim([-maxc,maxc]);
            		ylim([-maxc,maxc]);
            		set(gca, 'xtick', [-maxc,0,maxc]);
            		set(gca, 'ytick', [-maxc,0,maxc]);
            		xlabel('x (mm)', 'fontsize', 12, 'fontweight', 'bold');
            		ylabel('y (mm)', 'fontsize', 12, 'fontweight', 'bold');
            		set(gca, 'fontsize', 12);
            		title(sprintf('%s - 2D Geometry', obj.name_tag), 'fontweight', 'bold');
%            		grid on;
%            		set(gca, 'gridlinestyle', '--');
            		
                case '2d,trans,simple'
                    tc = getCrystalTransaxialLocations(obj);
                    if size(tc, 3) > 1
                        for n = 1 : size(tc, 2)
                            plot(squeeze(tc(1,n,:)), squeeze(tc(2,n,:)), '.');
                            if n==1
                                hold on;
                            end
                        end
                        
                        R = obj.system_parms.ring_diameter * 0.5 - 25;
                        nblock_trans = obj.system_parms.number_of_detector_modules_transaxial;
                        for n = 1 : nblock_trans
                        	a = 360/nblock_trans * (n-1) + ...
                        		obj.system_parms.detector_module_initial_angle_offset;
                        	text(R*cos(a*pi/180), R*sin(a*pi/180), num2str(n), ...
                        	'horizontalalignment', 'center');
                        end
                        axis square;
                        
                        hold off;
                    else
                        plot(tc(1,:), tc(2,:), '.');
                        R = obj.system_parms.ring_diameter * 0.5 - 25;
                        nblock_trans = obj.system_parms.number_of_detector_modules_transaxial;
                        for n = 1 : nblock_trans
                        	a = 360/nblock_trans * (n-1) + ...
                        		obj.system_parms.detector_module_initial_angle_offset;
                        	text(R*cos(a*pi/180), R*sin(a*pi/180), num2str(n), ...
                        	'horizontalalignment', 'center');
                        end
                        title(sprintf('%s - 2D Geometry', obj.name_tag), 'fontweight', 'bold');
                        axis square;
                    end
                
                case '2d,trans,simple,+lor'                    
                    tc = getCrystalTransaxialLocations(obj);
                    if size(tc, 3) > 1
                        for n = 1 : size(tc, 2)
                            plot(squeeze(tc(1,n,:)), squeeze(tc(2,n,:)), '.');
                            if n==1
                                hold on;
                            end
                        end
                        
                        axis square;
                        hold off;
                    else
                        plot(tc(1,:), tc(2,:), '.');
                        axis square;
                    end
                    xp = getDefaultSinogramCrystalPairs(obj);
                    np = obj.system_parms.number_of_projections_per_angle;
                    xp = reshape(xp, 2, np, size(xp,2)/np);
                    hold on;
                    for i = 1 : floor(size(xp,3)) : size(xp, 3)
                    for n = 1 : np
                        xtal_pair = xp(:, n, i);
                        p0 = tc(:, xtal_pair(1));
                        p1 = tc(:, xtal_pair(2));
                        line([p0(1), p1(1)], [p0(2), p1(2)]);
                    end
                    end
                    hold off;
                    title('Only show LORs at the 1st angle');
                    
                case '3d,simple'
                    tc = getCrystalTransaxialLocations(obj);
                    ro = getCrystalRingOffsets(obj);
                    for z = 1 : length(ro)
                        if size(tc, 3) > 1
                            for n = 1 : size(tc, 2)
                                plot3(squeeze(tc(1,n,:)), ...
                                    squeeze(tc(2,n,:)), ...
                                    ro(z) * ones(size(tc,3),1), '.');
                                if n==1 && z==1
                                    hold on;
                                end
                            end
                            axis square;
                        else
                            plot3(tc(1,:), tc(2,:), ro(z) * ones(size(tc,2),1), '.');
                            axis square;
                            if z==1
                                hold on;
                            end
                        end                        
                    end
                    hold off;
                otherwise
                    error('unknown options! valid options: "2d,geom", "2d,trans,simple", "2d,trans,simple,+lor", "3d,simple"');
            end            
        end        
    end
    
    methods
        function xtal_pairs = getDefaultSinogramCrystalPairs(obj)
        % create crystal pairs according my own numbering scheme
        %function xtal_pairs = getDefaultCrystalPairs
        % 
            
            % always equal to total number of crystals per ring divided by
            % 2
            nxtal_trans = obj.system_parms.crystal_array_size(1);
            xtal_num_trans_total = obj.system_parms.number_of_detector_modules_transaxial * ...
                            obj.system_parms.crystal_array_size(1);
            
            num_of_angles = xtal_num_trans_total / 2;
            
            if obj.system_parms.number_of_projections_per_angle > (fix(xtal_num_trans_total/2)*2 - 2)
                error('too many projections!');
            end

            % case when detector module is not located exactly at 3-clock
            if obj.system_parms.detector_module_initial_angle_offset ~= 0
                disp('The 1st detector module is rotated by an angle!');
                disp('NOTE: For crystal pairing, this angle offset is never used directly!');
                disp('but always assume it is equal to 180 / (number_of_detector_modules_transaxial)');
                id0 = 0;
                id1 = num_of_angles;    
                nr = fix(num_of_angles/2) * 2 - 2;    
                h0 = zeros(nr,1);
                h1 = zeros(nr,1);
    
                for i=0:nr-1
       
                    id=id0 + fix(nr/2) - i - 1;
                    if id < 0
                        id = id + xtal_num_trans_total;
                    end
        
                h0(i+1)=id;        
                id = id1-fix(nr/2) + i;
                h1(i+1)=id;        
                end
            else
    
                id0 = fix(nxtal_trans / 2);
                odd = mod(nxtal_trans, 2) ~= 0;
    
                if odd
                    id1 = fix(xtal_num_trans_total / 2) + fix(nxtal_trans / 2);
                    nr = fix(num_of_angles / 2) * 2 - 1;
                else
                    id1 = fix(xtal_num_trans_total / 2) + fix(nxtal_trans / 2) - 1;
                    nr = fix(num_of_angles / 2) * 2;
                end

                h0 = zeros(nr,1);
                h1 = zeros(nr,1);

                for i=0:nr-1        
                    if odd
                        id = id0 + fix(nr/2) - i;
                    else
                        id = id0 + fix(nr/2) -1 - i;
                    end
        
                    if id < 0
                        id = id + xtal_num_trans_total; 
                    end
        
                    h0(i+1) = id;        
                    if odd
                        id = id1 - fix(nr/2) + i;
                    else
                        id = id1 - fix(nr/2) + 1 + i;
                    end        
                    h1(i+1) = id;        
                end
            end

            % pairing
            c = 1; k = 1;
            while c < length(h0)
                xtal_id1 = h0(c);
                xtal_id2 = h1(c);    
                xp_first_angle(:,k) = [xtal_id1; xtal_id2];
                k = k + 1;    
                if (c+1) <= length(h1)        
                    xtal_id1 = h0(c);
                    xtal_id2 = h1(c+1);        
                    xp_first_angle(:,k) = [xtal_id1; xtal_id2];
                    k = k + 1;
                end    
                c = c + 1;
            end

            %
            nn = size(xp_first_angle,2);
            num_of_projs_per_angle = obj.system_parms.number_of_projections_per_angle;
            xtal_pairs = zeros(2, num_of_angles * num_of_projs_per_angle);
            k = 1;
            for i=1:num_of_angles
                for j=1:num_of_projs_per_angle
                    pp = xp_first_angle(:, fix(nn / 2) - fix(num_of_projs_per_angle / 2) + j);        
                    p0 = mod(pp(1) + i-1, xtal_num_trans_total);
                    p1 = mod(pp(2) + i-1, xtal_num_trans_total);        
                    xtal_pairs(:,k) = [p0; p1];
                    k = k + 1;    
                end        
            end
            xtal_pairs = xtal_pairs + 1;
        end
    end
        
    methods(Access = 'public')
        function rp = getDefaultSinogramPlaneOrder(obj)
        %get default plane arrangement in terms of ring pair (different from Micxx-gram)
        %function rp = getDefaultSinogramPlaneOrder
        %
            nring = obj.getNumberOfCrystalRings();
            rp = zeros(2, nring*nring);
            offset = 0;
            for n = 1 : nring
                if n==1
                    rp(:,1:nring) = [1:nring; 1:nring];
                    offset = offset + nring;
                else
                    r_odd = [1:(nring-n+1); n:(nring)];
                    r_even = nring-r_odd+1;
                    nr = (nring-n+1)*2;
                    rp(:, offset + (1:2:nr)) = r_odd;
                    rp(:, offset + (2:2:nr)) = r_even;
                    offset = offset + nr;
                end
            end
        end
    end
    
       methods
        function xtal_pairs = getDefaultSinogramCrystalPairs_w_orientation(obj)
        % create crystal pairs according my own numbering scheme
        % % create crystal pairs according zhoujian's scheme and UIH
        % oritation criteria
        %function xtal_pairs = getDefaultCrystalPairs
        % 
            
            % always equal to total number of crystals per ring divided by
            % 2
            nxtal_trans = obj.system_parms.crystal_array_size(1);
            xtal_num_trans_total = obj.system_parms.number_of_detector_modules_transaxial * ...
                            obj.system_parms.crystal_array_size(1);
            
            num_of_angles = xtal_num_trans_total / 2;
            
            if obj.system_parms.number_of_projections_per_angle > (fix(xtal_num_trans_total/2)*2 - 2)
                error('too many projections!');
            end

            % case when detector module is not located exactly at 3-clock
            if obj.system_parms.detector_module_initial_angle_offset ~= 0
                disp('The 1st detector module is rotated by an angle!');
                disp('NOTE: For crystal pairing, this angle offset is never used directly!');
                disp('but always assume it is equal to 180 / (number_of_detector_modules_transaxial)');
                id0 = 0;
                id1 = num_of_angles;    
                nr = fix(num_of_angles/2) * 2 - 2;    
                h0 = zeros(nr,1);
                h1 = zeros(nr,1);
    
                for i=0:nr-1
       
                    id=id0 + fix(nr/2) - i - 1;
                    if id < 0
                        id = id + xtal_num_trans_total;
                    end
        
                h0(i+1)=id;        
                id = id1-fix(nr/2) + i;
                h1(i+1)=id;        
                end
            else
    
                id0 = fix(nxtal_trans / 2);
                odd = mod(nxtal_trans, 2) ~= 0;
    
                if odd
                    id1 = fix(xtal_num_trans_total / 2) + fix(nxtal_trans / 2);
                    nr = fix(num_of_angles / 2) * 2 - 1;
                else
                    id1 = fix(xtal_num_trans_total / 2) + fix(nxtal_trans / 2) - 1;
                    nr = fix(num_of_angles / 2) * 2;
                end

                h0 = zeros(nr,1);
                h1 = zeros(nr,1);

                for i=0:nr-1        
                    if odd
                        id = id0 + fix(nr/2) - i;
                    else
                        id = id0 + fix(nr/2) -1 - i;
                    end
        
                    if id < 0
                        id = id + xtal_num_trans_total; 
                    end
        
                    h0(i+1) = id;        
                    if odd
                        id = id1 - fix(nr/2) + i;
                    else
                        id = id1 - fix(nr/2) + 1 + i;
                    end        
                    h1(i+1) = id;        
                end
            end

            % pairing
            c = 1; k = 1;
            while c < length(h0)
                xtal_id1 = h0(c);
                xtal_id2 = h1(c);    
                xp_first_angle(:,k) = [xtal_id1; xtal_id2];
                k = k + 1;    
                if (c+1) <= length(h1)        
                    xtal_id1 = h0(c);
                    xtal_id2 = h1(c+1);        
                    xp_first_angle(:,k) = [xtal_id1; xtal_id2];
                    k = k + 1;
                end    
                c = c + 1;
            end

            %
            nn = size(xp_first_angle,2);
            num_of_projs_per_angle = obj.system_parms.number_of_projections_per_angle;
            xtal_pairs = zeros(2, num_of_angles * num_of_projs_per_angle);
            k = 1;
            tc = obj.getCrystalTransaxialLocations();

            for i=1:num_of_angles
                for j=1:num_of_projs_per_angle

                    pp = xp_first_angle(:, fix(nn / 2) - fix(num_of_projs_per_angle / 2) + j);
                    p0 = mod(pp(1) + i-1, xtal_num_trans_total);
                    p1 = mod(pp(2) + i-1, xtal_num_trans_total);    
%                    if p0==0 || p1==0
%                      disp('debug');
%                  end
                    event1_trans = tc(:, p0+1);
                    event2_trans = tc(:, p1+1);
                    %%%% simset listmode do not match with zhoujian format
                    event1_trans=flip(event1_trans);event1_trans(2)=-event1_trans(2);
                    event2_trans=flip(event2_trans);event2_trans(2)=-event2_trans(2);
                     %%%
                    if event1_trans(1)>event2_trans(1)
                        xtal_pairs(:,k) = [p0; p1];
                    elseif event1_trans(1)<event2_trans(1)
                        xtal_pairs(:,k) = [p1; p0];
                    elseif round(event1_trans(1),4)==round(event2_trans(1),4)
                        if  event1_trans(2)>event2_trans(2)
                            xtal_pairs(:,k) = [p0; p1];
                        else
                            xtal_pairs(:,k) = [p1; p0];
                        end
                    else
                        xtal_pairs(:,k) = [p0; p1];
                    end
            
                    k = k + 1;
                end
            end
            xtal_pairs = xtal_pairs + 1;
        end
       end
           
    methods(Access = 'public')
        function prjs = doListModeForwardProjectionNonTOF(obj, image, image_size, voxel_size, lmdata)
        %do listmode-based forward projection for non-TOF case
        %function prjs = doListModeForwardProjectionNonTOF(image, image_size, voxel_size, lmdata)
        %             
            image = reshape(image, image_size);
            %
            if (~isa(lmdata, 'int16')) && (~isa(lmdata, 'uint16'))
                error('invalid data type: lmdata must be int16 or uint16!');
            end
            %
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            switch obj.system_parms.projector
            	case 'Siddon'
        	    	prjs = fproj_mt(image, image_size, voxel_size, tc, to, lmdata);
    	        case 'Bresenham'
	            	prjs = fproj_mt_bresenham(image, image_size, voxel_size, tc, to, lmdata);
	            case 'Linterp'
	            	prjs = fproj_mt_linterp(image, image_size, voxel_size, tc, to, lmdata);	
	            otherwise
	            	error('invalid projector!');
            end
        end

        function prjs = doListModeForwardProjectionTOF(obj, image, image_size, voxel_size, lmdata)
        %do listmode-based forward projection for TOF case
        %function prjs = doListModeForwardProjectionNonTOF(image, image_size, voxel_size, lmdata)
        %             
            image = reshape(image, image_size);
            %
            if ~isa(lmdata, 'int16')
                error('invalid data type: lmdata must be int16');
            end
            %
            if isempty(obj.system_parms.tof_info)
                error('unable to run TOF projector on non-TOF scanner!');
            end
            
            %
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            switch obj.system_parms.projector
    	        case 'Siddon'
        		    prjs = fproj_tof_mt(image, image_size, voxel_size, tc, to, lmdata, ...
            		    obj.system_parms.tof_info);
	            case 'Bresenham'
	        	    prjs = fproj_tof_mt_bresenham(image, image_size, voxel_size, tc, to, lmdata, ...
    	        	    obj.system_parms.tof_info);
    	        case 'Linterp'
	        	    prjs = fproj_tof_mt_linterp(image, image_size, voxel_size, tc, to, lmdata, ...
    	        	    obj.system_parms.tof_info);	    
    	        otherwise
    	        	error('invalid projector!');
            end
        end
        
        function image = doListModeBackProjectionNonTOF(obj, prjs, image_size, voxel_size, lmdata)
        %do listmode-based backprojection for non-TOF case
        %function image = doListModeBackProjectionNonTOF(projs, image_size, voxel_size, lmdata)
        %
        %
            if (~isa(lmdata, 'int16')) && (~isa(lmdata, 'uint16'))
                error('invalid data type: lmdata must be int16 or uint16!');
            end
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            switch obj.system_parms.projector
	            case 'Siddon'
    		        image = bproj_mt(prjs, image_size, voxel_size, tc, to, lmdata);            
            	case 'Bresenham'
            		image = bproj_mt_bresenham(prjs, image_size, voxel_size, tc, to, lmdata);
            	case 'Linterp'
            		image = bproj_mt_linterp(prjs, image_size, voxel_size, tc, to, lmdata);	
            	otherwise
            		error('invalid projector!');
            end
            image = reshape(image, image_size);
        end
        
        function image = doListModeBackProjectionNonTOFSingleThread(obj, prjs, image_size, voxel_size, lmdata)
        %do listmode-based backprojection for non-TOF case using single CPU
        %function image = doListModeBackProjectionNonTOFSingleThread(projs, image_size, voxel_size, lmdata)
        %
        %
            if (~isa(lmdata, 'int16')) && (~isa(lmdata, 'uint16'))
                error('invalid data type: lmdata must be int16 or uint16!');
            end
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            switch obj.system_parms.projector
	            case 'Siddon'
    		        image = bproj_st(prjs, image_size, voxel_size, tc, to, lmdata);            
            	case 'Bresenham'
            		image = bproj_st_bresenham(prjs, image_size, voxel_size, tc, to, lmdata);
            	case 'Linterp'
            		image = bproj_st_linterp(prjs, image_size, voxel_size, tc, to, lmdata);	
            	otherwise
            		error('invalid projector!');
            end
            image = reshape(image, image_size);
        end
        
        function image = doSquaredListModeBackProjectionNonTOF(obj, prjs, image_size, voxel_size, lmdata)
        %do listmode-based backprojection where the weight along each LOR is squared (Currently Siddon only)
        %function image = doSquaredListModeBackProjectionNonTOF(projs, image_size, voxel_size, lmdata)
        %
        %
            if (~isa(lmdata, 'int16')) && (~isa(lmdata, 'uint16'))
                error('invalid data type: lmdata must be int16 or uint16!');
            end
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            image = bproj_sq_mt(prjs, image_size, voxel_size, tc, to, lmdata);            
            image = reshape(image, image_size);
        end

        function image = doSquaredListModeBackProjectionTOF(obj, prjs, image_size, voxel_size, lmdata)
        %do listmode-based backprojection where the weight along each LOR is squared (Currently Siddon only!)
        %function image = doSquaredListModeBackProjectionTOF(projs, image_size, voxel_size, lmdata)
        %
        %
            if (~isa(lmdata, 'int16')) && (~isa(lmdata, 'uint16'))
                error('invalid data type: lmdata must be int16 or uint16!');
            end
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            image = bproj_sq_tof_mt(prjs, image_size, voxel_size, tc, to, lmdata, obj.system_parms.tof_info);            
            image = reshape(image, image_size);
        end
        
        function image = doListModeBackProjectionTOF(obj, prjs, image_size, voxel_size, lmdata)
        %do listmode-based backprojection for TOF case
        %function image = doListModeBackProjectionTOF(projs, image_size, voxel_size, lmdata)
        %
        %
            if ~isa(lmdata, 'int16')
                error('invalid data type: lmdata must be int16 or uint16!');
            end
            
            if isempty(obj.system_parms.tof_info)
                error('unable to run TOF projector on non-TOF scanner!');
            end
            
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            switch obj.system_parms.projector
	            case 'Siddon'
    	    	    image = bproj_tof_mt(prjs, image_size, voxel_size, tc, to, lmdata, ...
        	    	    obj.system_parms.tof_info);
            	case 'Bresenham'
    	    	    image = bproj_tof_mt_bresenham(prjs, image_size, voxel_size, tc, to, lmdata, ...
        	    	    obj.system_parms.tof_info);
        	    case 'Linterp'
    	    	    image = bproj_tof_mt_linterp(prjs, image_size, voxel_size, tc, to, lmdata, ...
        	    	    obj.system_parms.tof_info);	    
            	otherwise
            		error('invalid projector');
            end
            image = reshape(image, image_size);
            
        end
        
        function senimg = calculatePureGeometricSensitivity(obj, image_size, voxel_size, max_ringdiff)
        %compute a pure geometrical sensitivity image (currently fast_version is not available)
        %function senimg = calculatePureGeometricSensitivity(image_size, voxel_size, max_ringdiff)
        %
            nring = obj.getNumberOfCrystalRings();
        	if nargin < 4
        		max_ringdiff = nring-1;
        	end
        	if isempty(nring)
        		max_ringdiff = nring-1;
        	end
        	
        	if max_ringdiff > nring-1 || max_ringdiff < 0
        		error('invalid maximum ring difference!');
        	end
        	        	
            fprintf('calculating pure geometrical sensitivity for maximum ring difference = %d...\n', ...
            	max_ringdiff);
            fprintf('image size: %d x %d x %d, voxel size: %f x %f x %f (mm^3)\n', ...
                image_size, voxel_size);
        
            xp = int16(obj.getDefaultSinogramCrystalPairs());
            
            fast_version = false; % not used any more
            if fast_version % fast version
            	disp('use fast calculation');
            	axial_ratio = floor(obj.system_parms.crystal_size(3) / voxel_size(3));
            	fprintf('assumed axial ratio = %.2f ...\n', axial_ratio);
            	if (axial_ratio * voxel_size(3)) ~= obj.system_parms.crystal_size(3)
	            	error('NOTE: voxel axial size must be proportional to crystal axial pitch size!');
            	end
            	%
	            i = int16(0);            
	            lmdata = repmat(i, 5, size(xp,2));
    	        lmdata(1,:) = xp(1,:)-1;
    	        lmdata(3,:) = xp(2,:)-1;
    	        lmdata(4,:) = i;
    	        % in debugging
    	        for n = 1 : nring
    	        %	fprintf('processing ring #%d ...\n', n);    	        	
    	        %	lmdata(4,:) = int16(n-1);
    	        %	s0 = obj.doListModeBackProjectionNonTOF([], image_size, voxel_size, lmdata);
    	        %	s1 = s0;
    	        %	for m = 1 : (nring-n+1)    	        		
				 %       st = (m-1)*axial_ratio;
				 %       ss = zeros(image_size);           
				 %       ss(:, :, (st+1):end) = s0(:, :, 1:(end-st));            
				 %       s1 = s1 + ss;    	        	
    	        %	end    	        	
    	        	%
		        %	if n == 1
			    %    	senimg = s1;
			    %    else
			    %        senimg = s1 + flipdim(s1,3);
			    %   end
    	        end
            else % slow version
	            i = int16(0);            
	            m = repmat(int16(0:nring-1),size(xp,2),1);
            
	            lmdata = repmat(i, 5, size(xp,2) * nring);
    	        lmdata(1,:) = repmat(xp(1,:)-1, 1, nring);
    	        lmdata(3,:) = repmat(xp(2,:)-1, 1, nring);
    	        lmdata(4,:) = m(:);
            
    	        senimg = 0;
    	        for n = 1 : nring
	    	        fprintf('processing ring #%d ...\n', n); 
    	            lmdata(2,:) = int16(n-1);
    	            if max_ringdiff < nring-1
       		            lm0 = lmdata(:, abs(lmdata(2,:)-lmdata(4,:)) < max_ringdiff);
    		            s0 = obj.doListModeBackProjectionNonTOF([], image_size, voxel_size, lm0);
    	            else
	    	            s0 = obj.doListModeBackProjectionNonTOF([], image_size, voxel_size, lmdata);
    	            end
    	            senimg = senimg + s0;
    	        end
            
    	        senimg = reshape(senimg, image_size);            
            end
        end






        function senimg = calculatePureGeometricSensitivity_Att_NormlmBP_zhou(obj, raw_data, lmdata, att_image, att_image_size, att_voxel_size, image_size, voxel_size)
            
            att_image = reshape(att_image, att_image_size);
            
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            
            
            %           fprintf('calculating pure geometrical line integral for...\n');
            %
            %           raw_data = obj.doListModeForwardProjectionNonTOF(cylinderimage, cylinderimage_size, cylindervoxel_size, lmdata);
            
            
            fprintf('calculating forward projection through att_image \n');
            
            prjs_att = fproj_mt(att_image, att_image_size, att_voxel_size, tc, to, lmdata);
            
            whos 

            t = size(prjs_att)

            whos prjs_att            
            
            ratio_MC_Ana = reshape(raw_data, t(1), t(2));
            
            % clear raw_data t

            whos 
            
            fprintf('calculating pure geometrical sensitivity ...\n');
            fprintf('image size: %d x %d x %d, voxel size: %f x %f x %f (mm^3)\n', image_size, voxel_size);
            
            senimg  = obj.doListModeBackProjectionNonTOF( exp(-prjs_att) ./ double(ratio_MC_Ana) , image_size, voxel_size, lmdata);
            
            senimg(senimg==Inf) = 0;
            
            fprintf('reshape(senimg, image_size) ...\n');
            senimg = reshape(senimg, image_size);
            
        end
        
        




        
    end









    
    methods
        function sa = calculateApproximatedSolidAngle(obj)

            nrad = obj.system_parms.number_of_projections_per_angle;
            nang = obj.getDefaultNumberOfAngles;
            nring = obj.getNumberOfCrystalRings;
            xp = obj.getDefaultSinogramCrystalPairs;
            xp = int16(xp);
            %
            rp = obj.getDefaultSinogramPlaneOrder;
            rp = int16(rp);

            i = int16(0);
            lm = repmat(i, 5, size(xp,2) * nring * nring); 
            lm(1,:) = repmat(xp(1,:)-1, 1, nring * nring);
            ii = repmat(rp(1,:)-1, size(xp,2), 1);
            lm(2,:) = ii(:);
            lm(3,:) = repmat(xp(2,:)-1, 1, nring * nring);
            ii = repmat(rp(2,:)-1, size(xp,2), 1);
            lm(4,:) = ii(:);
            tc = obj.getCrystalTransaxialLocations;
            rc = obj.getCrystalRingOffsets;
            p1 = tc(:,lm(1,:)+1);
            p1(3,:) = rc(lm(2,:)+1);
            p2 = tc(:,lm(3,:)+1);
            p2(3,:) = rc(lm(4,:)+1);
            sa = 1./sum((p1-p2).^2,1);

        end
    end

    % sinogram-based projector
    methods
        function sino = doSinogramForwardProjection(obj, ...
                image, image_size, voxel_size, angle_subsets, fast)
        %do sinogram-based forward projection        
        %function sino = doSinogramForwardProjection(image, image_size, voxel_size, angle_subsets, fast)
        %        
            if nargin < 5
                fast = false;
            else
                if isempty(fast)
                    fast = false;
                end
            end
            
            nrad = obj.system_parms.number_of_projections_per_angle;
            nang = obj.getDefaultNumberOfAngles;
            nring = obj.getNumberOfCrystalRings;
            xp = obj.getDefaultSinogramCrystalPairs;
            xp = int16(xp);
            % pick angle subsets
            if ~isempty(angle_subsets)
                xp = reshape(xp, 2, nrad, nang);
                xp = xp(:, :, angle_subsets);
                xp = reshape(xp, 2, nrad * length(angle_subsets));
            end
            %
            rp = obj.getDefaultSinogramPlaneOrder;
            rp = int16(rp);
            
            if fast
                i = int16(0);
                lm = repmat(i, 5, size(xp,2) * nring * nring); 
                lm(1,:) = repmat(xp(1,:)-1, 1, nring * nring);
                ii = repmat(rp(1,:)-1, size(xp,2), 1);
                lm(2,:) = ii(:);
                lm(3,:) = repmat(xp(2,:)-1, 1, nring * nring);
                ii = repmat(rp(2,:)-1, size(xp,2), 1);
                lm(4,:) = ii(:);
                sino = obj.doListModeForwardProjectionNonTOF(image, image_size, voxel_size, lm);
            else
                i = int16(0);
                lm = repmat(i, 5, size(xp,2) * nring); 
                lm(1,:) = repmat(xp(1,:)-1, 1, nring);
                lm(3,:) = repmat(xp(2,:)-1, 1, nring);
                sino = zeros(size(xp,2)*nring, nring);
                
                rp = reshape(rp, 2, nring, nring);                
                for n = 1 : nring
                    r = rp(:, :, n);
                    %
                    ii = repmat(r(1,:)-1, size(xp,2), 1);
                    lm(2,:) = ii(:);
                    ii = repmat(r(2,:)-1, size(xp,2), 1);
                    lm(4,:) = ii(:);
                    %
                    s = obj.doListModeForwardProjectionNonTOF(image, image_size, voxel_size, lm);
                    sino(:, n) = s;
                end
            end
            
            nb = obj.system_parms.number_of_projections_per_angle;
            if isempty(angle_subsets)
                na = obj.getDefaultNumberOfAngles;
            else
                na = length(angle_subsets);
            end
            sino = reshape(sino, nb, na, nring*nring);
            
        end
        
        function image = doSinogramBackProjection(obj, ...
                sino, image_size, voxel_size, angle_subsets, fast)
        %do sinogram-based backprojection                    
        %function image = doSinogramBackProjection(sino, image_size, voxel_size, angle_subsets, fast)
            
            if nargin < 5
                fast = false;
            else
                if isempty(fast)
                    fast = false;
                end
            end
                        
            nring = obj.getNumberOfCrystalRings;
            nrad = obj.system_parms.number_of_projections_per_angle;
            nang = obj.getDefaultNumberOfAngles;
            xp = obj.getDefaultSinogramCrystalPairs;
            xp = int16(xp);
            % pick angle subsets
            if ~isempty(angle_subsets)
                xp = reshape(xp, 2, nrad, nang);
                xp = xp(:, :, angle_subsets);
                xp = reshape(xp, 2, nrad * length(angle_subsets));
            end
            %
            rp = obj.getDefaultSinogramPlaneOrder;
            rp = int16(rp);
            image = 0;
            
            if fast
                i = int16(0);
                lm = repmat(i, 5, size(xp,2) * nring * nring); 
                lm(1,:) = repmat(xp(1,:)-1, 1, nring * nring);
                ii = repmat(rp(1,:)-1, size(xp,2), 1);
                lm(2,:) = ii(:);
                lm(3,:) = repmat(xp(2,:)-1, 1, nring * nring);
                ii = repmat(rp(2,:)-1, size(xp,2), 1);
                lm(4,:) = ii(:);                
                image = obj.doListModeBackProjectionNonTOF(sino, image_size, voxel_size, lm);
            else
                i = int16(0);
                lm = repmat(i, 5, size(xp,2) * nring); 
                lm(1,:) = repmat(xp(1,:)-1, 1, nring);
                lm(3,:) = repmat(xp(2,:)-1, 1, nring);
                if ~isempty(sino)
                    sino = reshape(sino, size(xp,2)*nring, nring);
                end
                
                rp = reshape(rp, 2, nring, nring);                
                for n = 1 : nring
                    r = rp(:, :, n);
                    %
                    ii = repmat(r(1,:)-1, size(xp,2), 1);
                    lm(2,:) = ii(:);
                    ii = repmat(r(2,:)-1, size(xp,2), 1);
                    lm(4,:) = ii(:);
                    %
                    if ~isempty(sino)
                        bp = obj.doListModeBackProjectionNonTOF(sino(:,n), ...
                            image_size, voxel_size, lm);
                    else
                        bp = obj.doListModeBackProjectionNonTOF([], ...
                            image_size, voxel_size, lm);
                    end
                    image = image + bp;
                end
            end
            image = reshape(image, image_size);
        end
              	
       	function image = calculateFisherInformationMatrixDiagonal(obj, ...
                			diag_weight, image_size, voxel_size, fast)
		%compute the diagonal of the Fisher information matrix for a non-TOF system
       	%function image = calculateFisherInformationMatrixDiagonal(diag_weight, image_size, voxel_size, fast)
       	%
		%
		%                			
            nring = obj.getNumberOfCrystalRings;
            nrad = obj.system_parms.number_of_projections_per_angle;
            nang = obj.getDefaultNumberOfAngles;
            xp = obj.getDefaultSinogramCrystalPairs;
            xp = int16(xp);
            
            rp = obj.getDefaultSinogramPlaneOrder;
            rp = int16(rp);
            image = 0;
            
            if fast
            	i = int16(0);
            	lm = repmat(i, 5, size(xp,2) * nring * nring); 
            	lm(1,:) = repmat(xp(1,:)-1, 1, nring * nring);
            	ii = repmat(rp(1,:)-1, size(xp,2), 1);
           		lm(2,:) = ii(:);
          		lm(3,:) = repmat(xp(2,:)-1, 1, nring * nring);
           		ii = repmat(rp(2,:)-1, size(xp,2), 1);
           		lm(4,:) = ii(:);                
           		image = obj.doSquaredListModeBackProjectionNonTOF(diag_weight, ...
           			image_size, voxel_size, lm);
           	else
           		error('not implemented!');
           	end
           	image = reshape(image, image_size);
        end
        
        function lm = dumpSinogramBasedLORindex(obj)        
        	nring = obj.getNumberOfCrystalRings;
            xp = obj.getDefaultSinogramCrystalPairs;
            xp = int16(xp);
            %
            rp = obj.getDefaultSinogramPlaneOrder;
            rp = int16(rp);        
        	i = int16(0);
            lm = repmat(i, 5, size(xp,2) * nring * nring);
            lm(1,:) = repmat(xp(1,:)-1, 1, nring * nring);
            ii = repmat(rp(1,:)-1, size(xp,2), 1);
            lm(2,:) = ii(:);
            lm(3,:) = repmat(xp(2,:)-1, 1, nring * nring);
            ii = repmat(rp(2,:)-1, size(xp,2), 1);
            lm(4,:) = ii(:);
        end
        
    end
    
    methods 
        
        function mask = getImageMask(obj, image_size, radius_unit_in_pixel)
        %create a circular mask
        %function mask = getImageMask(image_size, radius_unit_in_pixel)
        %
            t1 = (0:image_size(1)-1) - image_size(1)/2 + 0.5;
            t2 = (0:image_size(2)-1) - image_size(2)/2 + 0.5;
            [I, J] = ndgrid(t1, t2);
            m0 = sqrt(I.^2 + J.^2) < radius_unit_in_pixel;
            if length(image_size) == 3
                mask = repmat(m0(:), 1, image_size(3));
                mask = reshape(mask, image_size);
            else
                mask = m0;
            end
        end
        
        function G = getSimple2DMatrix(obj, image_size, voxel_size)
        %create a 2d system matrix
        %function G = getSimple2DMatrix(image_size, voxel_size)
        %
        	if obj.getNumberOfCrystalRings > 1
        		error('sorry, only for 2D scanner!');
        	end
        	if length(image_size)<3
        		if length(image_size)<2
        			error('invalid image dimension!');
        		else
	        		image_size(3) = 1;
        		end
        	end
        	if length(voxel_size)<3
        		if length(voxel_size)<2
        			error('invalid voxel size!');
        		else
        			voxel_size(3) = 1;
        		end
        	end
        	xp = obj.getDefaultSinogramCrystalPairs;
            xp = int16(xp);
            lm = int16(zeros(5, size(xp,2)));
            lm(1,:) = xp(1,:)-1;
            lm(3,:) = xp(2,:)-1; 
            nz0 = 30000000;           
            ii=zeros(nz0, 1);
            jj=zeros(nz0, 1);
            ss=zeros(nz0, 1);
            offset = 0;
            for n = 1 : size(lm,2)
            	if mod(n, obj.system_parms.number_of_projections_per_angle) == 0
            		fprintf('processing #%d ... (nnz=%d)\n', n, offset);
            	end
            	bp = obj.doListModeBackProjectionNonTOFSingleThread([], ...
                            image_size, voxel_size, lm(:,n));
                lor_nz = nnz(bp); 
				if lor_nz > 0
	                [i] = find(bp(:));
					ii([1:lor_nz] + offset) = i;
					ss([1:lor_nz] + offset) = bp(i);
					jj([1:lor_nz] + offset) = n*ones(size(i));
					offset = offset + lor_nz;
				end                
            end
            ii=ii(1:offset);
            jj=jj(1:offset);
            ss=ss(1:offset);            
            G = sparse(jj, ii, ss, size(lm,2), prod(image_size(1:2)));
    	end    
        
        function xmat = reconSinogramMLEM(obj, sino, image_size, voxel_size, ...
                nrm_fac, rnd_fac, resolution_model, ...
                number_of_subsets, number_of_iterations, stepsize, ...
                initial_guess)
        %do sinogram-based image reconstruction
        %function xmat = reconSinogramMLEM(sino, image_size, voxel_size, ...
        %        nrm_fac, rnd_fac, resolution_model, ...
        %        number_of_subsets, number_of_iterations, stepsize, ...
        %        initial_guess)
            
            fprintf('\nSinogram-based Reconstruction started at %s\n\n', datestr(now));
            
            if isempty(sino)
                error('need a sinogram!');
            end
            
            if isempty(initial_guess)
                x = ones(image_size);
            else
                x = initial_guess;
            end
            
            if number_of_subsets == 1
                disp('calculating sensitivity image ...');
                s0 = obj.doSinogramBackProjection(nrm_fac, ...
                    image_size, voxel_size, [], true);
                if  isempty(resolution_model)
                	senimg = s0;
                else
	                senimg = resolution_model.doBackwardBlurring(s0);
                end
            end                        
            
            nang = obj.getDefaultNumberOfAngles;
            
            % save output
            xmat = zeros([image_size, floor(number_of_iterations / stepsize)]);
            
            for n = 1 : number_of_iterations
                
                tic;
                fprintf('\n*** processing iteration #%d (%d in total) ***\n\n', n, number_of_iterations);
                
                for k = 1 : number_of_subsets
                    
                    fprintf('processing subset #%d (%d in total) ...\n', k, number_of_subsets);
                    
                    angle_subsets = k:number_of_subsets:nang;
                    
                    %
                    if ~isempty(nrm_fac)
                        sub_nrm_fac = nrm_fac(:,angle_subsets,:);
                    else
                        sub_nrm_fac = [];
                    end
                    if ~isempty(rnd_fac)
                        sub_rnd_fac = rnd_fac(:,angle_subsets,:);
                    else
                        sub_rnd_fac = 0;
                    end
                    %
                    if number_of_subsets > 1
	                    s0 = obj.doSinogramBackProjection(sub_nrm_fac, ...
    	                        image_size, voxel_size, angle_subsets, true);
    	                if isempty(resolution_model)
    	                	senimg = s0;
    	                else
	        	            senimg = resolution_model.doBackwardBlurring(s0);
        	            end
                    end
                    %
                    if 	isempty(resolution_model)
                    	xf = x;
                    else
	                    xf = resolution_model.doForwardBlurring(x);
                    end
                    
                    y = obj.doSinogramForwardProjection(xf, ...
                        image_size, voxel_size, angle_subsets, true);
                    
                    sub_sino = sino(:, angle_subsets, :);
                    
                    if isempty(sub_nrm_fac)
                        y = y + sub_rnd_fac;
                    else
                        y = y .* sub_nrm_fac + sub_rnd_fac;
                    end
                    yp = (sub_sino ./ y) .* sub_nrm_fac; yp(y == 0) = 0;
                    %
                    xb = obj.doSinogramBackProjection(yp, ...
                        image_size, voxel_size, angle_subsets, true);
                    
                    if isempty(resolution_model)
                    	xbf = xb;
                    else
                    	xbf = resolution_model.doBackwardBlurring(xb);
                    end
                    
                    % update
                    x = x ./ senimg .* xbf; 
                    x(senimg==0) = 0;
                    
                end
                toc;
                
                if mod(n, stepsize) == 0
                	disp('saved');
                    xmat(:,:,:,n/stepsize) = x;
                end
            end
            
        end
        
        function xmat = reconListModeMLEM(obj, lmdata, image_size, voxel_size, ...
                                       mul_fac, add_fac, sensitivity, resolution_model, ...
                                       number_of_subsets, number_of_iterations, stepsize, ...
                                       initial_guess, is_tof)
        %do listmode-based image reconstruction	
        %function xmat = reconListModeMLEM(lmdata, image_size, voxel_size, ...
        %                               mul_fac, add_fac, sensitivity, resolution_model, ...
        %                               number_of_subsets, number_of_iterations, stepsize, ...
        %                               initial_guess, is_tof)
        %
        	fprintf('\nListmode-based Reconstruction started at %s\n\n', datestr(now));
        	
            if ~isa(lmdata, 'int16') && ~isa(lmdata, 'uint16')
                error('invalid data type: lmdata must be int16 or uint16');
            end
            
            if isempty(mul_fac)
            	disp('[-]no multiplicative factors specified! use default value!');
            end

            if isempty(add_fac)
            	disp('[-]no additive factors specified! use default value!');
            	%add_fac = zeros(size(lmdata,2), 1);
            end
            
            if isempty(sensitivity)
                disp('[-]no sensitivity image specified! use default value!');
                sensitivity = ones(image_size);
            else
                disp('[+]sensitivity is detected!');
            end
            
            if ~isempty(resolution_model)
                senimg = resolution_model.doBackwardBlurring(sensitivity);
                disp('[+]resolution model is detected!');
            else            	
            	disp('[-]no resolution model specified! ignored!');
            	senimg = sensitivity;
            end
            
            if isempty(is_tof)
                is_tof = false;
                disp('[-]no TOF flag specified! treated as non-TOF recon!');
            end
                        
            if ~isempty(initial_guess)
                x = initial_guess;
            	disp('[+]use initial guess!');
            else
            	disp('[-]no initial guess specified! use default value!');
                x = ones(image_size);
            end
            
            fprintf('[+]number of prompts in total: %d\n', size(lmdata,2));
            fprintf('[+]image size: %d x %d x %d\n', image_size);
            fprintf('[+]voxel size: %.5f x %.5f x %.5f (mm^3)\n', voxel_size);

            % save output
            xmat = zeros([image_size, floor(number_of_iterations / stepsize)]);
            
            fprintf('\nready, set, go, OS-LM-MLEM ...\n');
            
            for n = 1 : number_of_iterations
                tic; 
                fprintf('\n*** processing iteration #%d (%d in total) ***\n\n', n, number_of_iterations);
                for k = 1 : number_of_subsets
                    
                    fprintf('processing subset #%d (%d in total) ...\n', k, number_of_subsets);
                    
                    sub_lmdata = lmdata(:, k:number_of_subsets:end);
                    
                    if isempty(mul_fac)
                        sub_mul_fac = 1.0;
                    else
                        sub_mul_fac = mul_fac(k:number_of_subsets:end);
                    end
                    
                    if isempty(add_fac)
                        sub_add_fac = 0.0;
                    else
                        sub_add_fac = add_fac(k:number_of_subsets:end);
                    end
                    
                    if ~isempty(resolution_model)
                        xf = resolution_model.doForwardBlurring(x);                    
                    else
                    	xf = x;
                    end
                    
                    if is_tof
                        yp = obj.doListModeForwardProjectionTOF(xf, ...
                                image_size, voxel_size, sub_lmdata);
                    else
                        yp = obj.doListModeForwardProjectionNonTOF(xf, ...
                                image_size, voxel_size, sub_lmdata);
                    end
                    
                    e = sub_mul_fac ./ (yp + sub_add_fac); 
                    e((yp==0) & (sub_add_fac==0)) = 0.0;
                    
                    if is_tof
                        xb = obj.doListModeBackProjectionTOF(e, image_size, ...
                                voxel_size, sub_lmdata);
                    else
                       xb = obj.doListModeBackProjectionNonTOF(e, image_size, ...
                                voxel_size, sub_lmdata);
                    end
                    
                    if ~isempty(resolution_model)
                        xf = resolution_model.doBackwardBlurring(xb);
                    else
                    	xf = xb;
                    end
                    
                    x(:) = x(:) ./ (senimg(:) / number_of_subsets) .* xf(:);
                    x(senimg == 0.0) = 0.0;    
                    
                end
                if mod(n, stepsize) == 0
                	disp('saved');
                    xmat(:,:,:,n/stepsize) = x;
                end                                
                toc;
            end
            
            %if isempty(xmat)
            %    xmat = x;
            %end            
        end

		function xmat = reconSinogramMAPEM(obj, sino, image_size, voxel_size, ...
                nrm_fac, rnd_fac, resolution_model, regularizer, ...
                number_of_subsets, number_of_iterations, stepsize, ...
                initial_guess, b0)
        %do sinogram-based MAP-EM image reconstruction
        %function xmat = reconSinogramMAPEM(obj, sino, image_size, voxel_size, ...
        %        nrm_fac, rnd_fac, resolution_model, regularizer, ...
        %        number_of_subsets, number_of_iterations, stepsize, ...
        %        initial_guess, b0)
            
            fprintf('\nSinogram-based MAP-EM Reconstruction started at %s\n\n', datestr(now));
            
            if isempty(sino)
                error('need a sinogram!');
            end
            
            if isempty(initial_guess)
                x = ones(image_size);
            else
                x = initial_guess;
            end
                        
            if number_of_subsets == 1
                disp('calculating sensitivity image ...');
                s0 = obj.doSinogramBackProjection(nrm_fac, ...
                    image_size, voxel_size, [], true);
                if  isempty(resolution_model)
                	senimg = s0;
                else
	                senimg = resolution_model.doBackwardBlurring(s0);
                end
            end                        
            
            nang = obj.getDefaultNumberOfAngles;
            
            % save output
            xmat = zeros([image_size, floor(number_of_iterations / stepsize)]);
            
            
            for n = 1 : number_of_iterations
                
                tic;
                fprintf('\n*** processing iteration #%d (%d in total) ***\n\n', n, number_of_iterations);
                
                for k = 1 : number_of_subsets
                    
                    fprintf('processing subset #%d (%d in total) ...\n', k, number_of_subsets);
                    
                    angle_subsets = k:number_of_subsets:nang;
                    
                    %
                    if ~isempty(nrm_fac)
                        sub_nrm_fac = nrm_fac(:,angle_subsets,:);
                    else
                        sub_nrm_fac = [];
                    end
                    if ~isempty(rnd_fac)
                        sub_rnd_fac = rnd_fac(:,angle_subsets,:);
                    else
                        sub_rnd_fac = 0;
                    end
                    %
                    if number_of_subsets > 1
	                    s0 = obj.doSinogramBackProjection(sub_nrm_fac, ...
    	                        image_size, voxel_size, angle_subsets, true);
    	                if isempty(resolution_model)
    	                	senimg = s0;
    	                else
	        	            senimg = resolution_model.doBackwardBlurring(s0);
        	            end
                    end
                    %
                    if 	isempty(resolution_model)
                    	xf = x;
                    else
	                    xf = resolution_model.doForwardBlurring(x);
                    end
                    
                    y = obj.doSinogramForwardProjection(xf, ...
                        image_size, voxel_size, angle_subsets, true);
                    
                    sub_sino = sino(:, angle_subsets, :);
                    
                    if isempty(sub_nrm_fac)
                        y = y + sub_rnd_fac;
                    else
                        y = y .* sub_nrm_fac + sub_rnd_fac;
                    end
                    yp = (sub_sino ./ y) .* sub_nrm_fac; yp(y == 0) = 0;

                    % whos yp sub_sino y sub_nrm_fac

                    %
                    xb = obj.doSinogramBackProjection(yp, ...
                        image_size, voxel_size, angle_subsets, true);
                    
                    if isempty(resolution_model)
                    	xbf = xb;
                    else
                    	xbf = resolution_model.doBackwardBlurring(xb);
                    end
                    
                    % update
                    %
                    xem = x .* xbf;
                    if (b0 > 0) & (~isempty(regularizer)) % beta
                    
                    	if 1                    		
                    		[t1, t2] = regularizer.calculateGradientBasedOnOptimzationTransfer(reshape(x, image_size)); 
                    	%                   	
                	    	a = b0 * t1;
                    		b = -b0 * (t1 .* x - t2) + senimg;
                   	 		c = xem; % * number_of_subsets;
                   	 		%
                   	 		clear t1 t2;
                   	 		%
                    	%
	                    	det = sqrt(b.*b + 4 * a .* c);
    	                	d = 2*c ./ (b + det);
    	                	d(b + det <= 0) = 0;						
							x = d;
							x((a == 0) & (b == 0)) = 0;
							ii = (a == 0) & (b ~= 0);
							x(ii) = c(ii) ./ b(ii);	
						end
                    
                    else
                    
	                    x = xem ./ (senimg); 
    	                x(senimg==0) = 0;
                    
                    end
                    
                end
                toc;
                
                if mod(n, stepsize) == 0
                	disp('saved');
                    xmat(:,:,:,n/stepsize) = x;
                end
            end
            
        end    

		function xmat = reconSinogramBSREM(obj, sino, image_size, voxel_size, ...
                nrm_fac, rnd_fac, resolution_model, regularizer, ...
                number_of_subsets, number_of_iterations, stepsize, ...
                initial_guess, b0)
        %do sinogram-based BSREM image reconstruction
        %function xmat = reconSinogramBSREM(obj, sino, image_size, voxel_size, ...
        %        nrm_fac, rnd_fac, resolution_model, regularizer, ...
        %        number_of_subsets, number_of_iterations, stepsize, ...
        %        initial_guess, b0)
            
            fprintf('\nSinogram-based BSREM Reconstruction started at %s\n\n', datestr(now));
            
            if isempty(sino)
                error('need a sinogram!');
            end
            
            if isempty(initial_guess)
                x = ones(image_size);
            else
                x = initial_guess;
            end
            
            if number_of_subsets == 1
                disp('calculating sensitivity image ...');
                s0 = obj.doSinogramBackProjection(nrm_fac, ...
                    image_size, voxel_size, [], true);
                if  isempty(resolution_model)
                	senimg = s0;
                else
	                senimg = resolution_model.doBackwardBlurring(s0);
                end
            end                        
            
            nang = obj.getDefaultNumberOfAngles;
            
            % save output
            xmat = zeros([image_size, floor(number_of_iterations / stepsize)]);
            
            %
            c0 = 1/15;
            
            for n = 1 : number_of_iterations
                
                alpha = 1/ (n * c0 + 1);
                
                tic;
                fprintf('\n*** processing iteration #%d (%d in total) ***\n\n', n, number_of_iterations);
                
                for k = 1 : number_of_subsets
                    
                    fprintf('processing subset #%d (%d in total) ...\n', k, number_of_subsets);
                    
                    angle_subsets = k:number_of_subsets:nang;
                    
                    %
                    if ~isempty(nrm_fac)
                        sub_nrm_fac = nrm_fac(:,angle_subsets,:);
                    else
                        sub_nrm_fac = [];
                    end
                    if ~isempty(rnd_fac)
                        sub_rnd_fac = rnd_fac(:,angle_subsets,:);
                    else
                        sub_rnd_fac = 0;
                    end
                    %
                    if number_of_subsets > 1
	                    s0 = obj.doSinogramBackProjection(sub_nrm_fac, ...
    	                        image_size, voxel_size, angle_subsets, true);
    	                if isempty(resolution_model)
    	                	senimg = s0;
    	                else
	        	            senimg = resolution_model.doBackwardBlurring(s0);
        	            end
                    end
                    %
                    if 	isempty(resolution_model)
                    	xf = x;
                    else
	                    xf = resolution_model.doForwardBlurring(x);
                    end
                    
                    y = obj.doSinogramForwardProjection(xf, ...
                        image_size, voxel_size, angle_subsets, true);
                    
                    sub_sino = sino(:, angle_subsets, :);
                    
                    if isempty(sub_nrm_fac)
                        y = y + sub_rnd_fac;
                    else
                        y = y .* sub_nrm_fac + sub_rnd_fac;
                    end
                    yp = (sub_sino ./ y) .* sub_nrm_fac; yp(y == 0) = 0;
                    %
                    xb = obj.doSinogramBackProjection(yp, ...
                        image_size, voxel_size, angle_subsets, true);
                    
                    if isempty(resolution_model)
                    	xbf = xb;
                    else
                    	xbf = resolution_model.doBackwardBlurring(xb);
                    end
                    
                    % calc gradient for penalty
                    grad_p = regularizer.calculateGradient(x);
                    
                    
                    % update
                    %                    
                    x = x + alpha * (x ./ senimg / number_of_subsets) .* (xbf - b0 / number_of_subsets * grad_p);
                    x = max(x, 0.0);
                    
                end
                toc;
                
                if mod(n, stepsize) == 0
                	disp('saved');
                    xmat(:,:,:,n/stepsize) = x;
                end
            end
            
        end

        function xmat = reconListModeMAPEM(obj, lmdata, image_size, voxel_size, ...
                                       	   mul_fac, add_fac, sensitivity, resolution_model, regularizer, ...
                                       	   number_of_subsets, number_of_iterations, stepsize, ...
                                       	   initial_guess, is_tof, b0)
        %do listmode-based image reconstruction	
        %function xmat = reconListModeMAPEM(obj, lmdata, image_size, voxel_size, ...
        %                               	mul_fac, add_fac, sensitivity, resolution_model, regualrizer, ...
        %                               	number_of_subsets, number_of_iterations, stepsize, ...
        %                               	initial_guess, is_tof, b0)
        %
        	fprintf('\nListmode-based MAP-EM Reconstruction started at %s\n\n', datestr(now));
        	
            if ~isa(lmdata, 'int16') && ~isa(lmdata, 'uint16')
                error('invalid data type: lmdata must be int16 or uint16');
            end
            
            if isempty(mul_fac)
            	disp('[-]no multiplicative factors specified! use default value!');
            end

            if isempty(add_fac)
            	disp('[-]no additive factors specified! use default value!');
            	%add_fac = zeros(size(lmdata,2), 1);
            end
            
            if isempty(sensitivity)
                disp('[-]no sensitivity image specified! use default value!');
                sensitivity = ones(image_size);
            else
                disp('[+]sensitivity is detected!');
            end
            
            if ~isempty(resolution_model)
                senimg = resolution_model.doBackwardBlurring(sensitivity);
                disp('[+]resolution model is detected!');
            else            	
            	disp('[-]no resolution model specified! ignored!');
            	senimg = sensitivity;
            end
            
            if ~isempty(regularizer) 
            	disp('[+]regularizer is detected!');
            else
            	disp('[-]no regularizer is found! ignored! ML recon only!');
            end
            
            if isempty(is_tof)
                is_tof = false;
                disp('[-]no TOF flag specified! treated as non-TOF recon!');
            end
            
            if isempty(b0)
            	b0 = 0;
            	disp('[-]no hyperparameter specified! use 0!');
            else
            	fprintf('[+]regularizer strength: %.3e\n', b0);
            end
                        
            if ~isempty(initial_guess)
                x = initial_guess;
            	disp('[+]use initial guess!');
            else
            	disp('[-]no initial guess specified! use default value!');
                x = ones(image_size);
            end
            
            fprintf('[+]number of prompts in total: %d\n', size(lmdata,2));
            fprintf('[+]image size: %d x %d x %d\n', image_size);
            fprintf('[+]voxel size: %.5f x %.5f x %.5f (mm^3)\n', voxel_size);

            % save output
            xmat = zeros([image_size, floor(number_of_iterations / stepsize)]);
            
            fprintf('\nready, set, go, OS-LM-MLEM ...\n');
            
            for n = 1 : number_of_iterations
                tic; 
                fprintf('\n*** processing iteration #%d (%d in total) ***\n\n', n, number_of_iterations);
                for k = 1 : number_of_subsets
                    
                    fprintf('processing subset #%d (%d in total) ...\n', k, number_of_subsets);
                    
                    sub_lmdata = lmdata(:, k:number_of_subsets:end);
                    
                    if isempty(mul_fac)
                        sub_mul_fac = 1.0;
                    else
                        sub_mul_fac = mul_fac(k:number_of_subsets:end);
                    end
                    
                    if isempty(add_fac)
                        sub_add_fac = 0.0;
                    else
                        sub_add_fac = add_fac(k:number_of_subsets:end);
                    end
                    
                    if ~isempty(resolution_model)
                        xf = resolution_model.doForwardBlurring(x);                    
                    else
                    	xf = x;
                    end
                    
                    if is_tof
                        yp = obj.doListModeForwardProjectionTOF(xf, ...
                                image_size, voxel_size, sub_lmdata);
                    else
                        yp = obj.doListModeForwardProjectionNonTOF(xf, ...
                                image_size, voxel_size, sub_lmdata);
                    end
                    
                    e = sub_mul_fac ./ (yp + sub_add_fac); 
                    e((yp==0) & (sub_add_fac==0)) = 0.0;
                    
                    if is_tof
                        xb = obj.doListModeBackProjectionTOF(e, image_size, ...
                                voxel_size, sub_lmdata);
                    else
                       	xb = obj.doListModeBackProjectionNonTOF(e, image_size, ...
                                voxel_size, sub_lmdata);
                    end
                    
                    if ~isempty(resolution_model)
                        xf = resolution_model.doBackwardBlurring(xb);
                    else
                    	xf = xb;
                    end
                    
                    %
                    xem = x .* xf;
                    if (b0 > 0) & (~isempty(regularizer)) % beta
                    
                    if 1
                    	[t1, t2] = regularizer.calculateGradientBasedOnOptimzationTransfer(reshape(x, image_size)); 
                    	%                   	
                    	a = b0 * t1;
                    	b = -b0 * (t1 .* x - t2) + senimg;
                    	c = xem * number_of_subsets;
                    	%
                    	det = sqrt(b.*b + 4 * a .* c);
                    	d = 2*c ./ (b + det);
                    	d(b + det <= 0) = 0;						
						x = d;
						x((a == 0) & (b == 0)) = 0;
						ii = (a == 0) & (b ~= 0);
						x(ii) = c(ii) ./ b(ii);			
					else % test!
%						b=(xm - senimg/number_of_subsets/b0);
%						c=-xem/b0;
%						det = sqrt(b.^2 - 4*c);
%						x = (b + det)/2;
%						x(det==0) = 0;
						xm = medfilt3(x, [3 3 3]);
						d = (x - xm);% ./ xm; d(xm == 0) = 0;
						x(:) = xem(:) ./ (senimg(:) / number_of_subsets + b0 * d(:));
						x = max(x,0);
					end
                    
                    else                 
                    %
	                    x(:) = xem(:) ./ (senimg(:) / number_of_subsets);
    	                x(senimg == 0.0) = 0.0;                        
                    end
                end
                if mod(n, stepsize) == 0
                	disp('saved');
                    xmat(:,:,:,n/stepsize) = x;
                end                                
                toc;
            end
            
            %if isempty(xmat)
            %    xmat = x;
            %end            
        end
    end
    
    methods
        function xmat = reconListModeMAPEM2(obj, lmdata, image_size, voxel_size, ...
                                            mul_fac, add_fac, sensitivity, resolution_model, regularizer, ...
                                            number_of_subsets, number_of_iterations, stepsize, ...
                                            initial_guess, is_tof, b0, g0)
        %do listmode-based image reconstruction (hybrid penalty)
        %function xmat = reconListModeMAPEM2(obj, lmdata, image_size, voxel_size, ...
        %                                   mul_fac, add_fac, sensitivity, resolution_model, regualrizer, ...
        %                                   number_of_subsets, number_of_iterations, stepsize, ...
        %                                   initial_guess, is_tof, b0)
        %
            fprintf('\nListmode-based PL-ADMM Reconstruction (hybrid penalty), %s\n\n', datestr(now));
            
            if ~isa(lmdata, 'int16') && ~isa(lmdata, 'uint16')
                error('invalid data type: lmdata must be int16 or uint16');
            end
            
            if isempty(mul_fac)
                disp('[-]no multiplicative factors specified! use default value!');
            end

            if isempty(add_fac)
                disp('[-]no additive factors specified! use default value!');
                %add_fac = zeros(size(lmdata,2), 1);
            end
            
            if isempty(sensitivity)
                disp('[-]no sensitivity image specified! use default value!');
                sensitivity = ones(image_size);
            else
                disp('[+]sensitivity is detected!');
            end
            
            if ~isempty(resolution_model)
                senimg = resolution_model.doBackwardBlurring(sensitivity);
                disp('[+]resolution model is detected!');
            else                
                disp('[-]no resolution model specified! ignored!');
                senimg = sensitivity;
            end
            
            if ~isempty(regularizer) 
                disp('[+]regularizer is detected!');
            else
                disp('[-]no regularizer is found! ignored! ML recon only!');
            end
            
            if isempty(is_tof)
                is_tof = false;
                disp('[-]no TOF flag specified! treated as non-TOF recon!');
            end
            
            if isempty(b0)
                b0 = 0;
                disp('[-]no hyperparameter specified! use 0!');
            else
                fprintf('[+]regularizer strength: %.3e, %.3e\n', b0, g0);
            end
                        
            if ~isempty(initial_guess)
                x = initial_guess;
                disp('[+]use initial guess!');
            else
                disp('[-]no initial guess specified! use default value!');
                x = ones(image_size);
            end
            
            fprintf('[+]number of prompts in total: %d\n', size(lmdata,2));
            fprintf('[+]image size: %d x %d x %d\n', image_size);
            fprintf('[+]voxel size: %.5f x %.5f x %.5f (mm^3)\n', voxel_size);

            % save output
            xmat = zeros([image_size, floor(number_of_iterations / stepsize)]);
            
            %
            rho = 1e8,
            z = zeros(image_size); 
            eta = zeros(image_size);

            fprintf('\nready, set, go, OS-LM-MLEM ...\n');
            
            for n = 1 : number_of_iterations

                %
                v = mirt_idctn(z) + eta;

                tic; 
                fprintf('\n*** processing iteration #%d (%d in total) ***\n\n', n, number_of_iterations);
                for k = 1 : number_of_subsets
                    
                    fprintf('processing subset #%d (%d in total) ...\n', k, number_of_subsets);
                    
                    sub_lmdata = lmdata(:, k:number_of_subsets:end);
                    
                    if isempty(mul_fac)
                        sub_mul_fac = 1.0;
                    else
                        sub_mul_fac = mul_fac(k:number_of_subsets:end);
                    end
                    
                    if isempty(add_fac)
                        sub_add_fac = 0.0;
                    else
                        sub_add_fac = add_fac(k:number_of_subsets:end);
                    end
                    
                    if ~isempty(resolution_model)
                        xf = resolution_model.doForwardBlurring(x);                    
                    else
                        xf = x;
                    end
                    
                    if is_tof
                        yp = obj.doListModeForwardProjectionTOF(xf, ...
                                image_size, voxel_size, sub_lmdata);
                    else
                        yp = obj.doListModeForwardProjectionNonTOF(xf, ...
                                image_size, voxel_size, sub_lmdata);
                    end
                    
                    e = sub_mul_fac ./ (yp + sub_add_fac); 
                    e((yp==0) & (sub_add_fac==0)) = 0.0;
                    
                    if is_tof
                        xb = obj.doListModeBackProjectionTOF(e, image_size, ...
                                voxel_size, sub_lmdata);
                    else
                        xb = obj.doListModeBackProjectionNonTOF(e, image_size, ...
                                voxel_size, sub_lmdata);
                    end
                    
                    if ~isempty(resolution_model)
                        xf = resolution_model.doBackwardBlurring(xb);
                    else
                        xf = xb;
                    end
                    
                    %
                    xem = x .* xf;

                    %
                    [t1, t2] = regularizer.calculateGradientBasedOnOptimzationTransfer(reshape(x, image_size)); 

                    a = b0 * t1 + rho;
                    b = -b0 * (t1 .* x - t2) + senimg - rho * v;
                    c = xem * number_of_subsets;
                    %
                    det = sqrt(b.*b + 4 * a .* c);
                    d = 2*c ./ (b + det);
                    d(b + det <= 0) = 0;                        
                    x = d;
                    x((a == 0) & (b == 0)) = 0;
                    ii = (a == 0) & (b ~= 0);
                    x(ii) = c(ii) ./ b(ii); 
                    x = max(x,0);

%                    x(:) = xem(:) ./ (senimg(:) / number_of_subsets);
%                    x = reshape(x, image_size);
                    x(senimg == 0.0) = 0.0;                                            

                end 

                % update z and eta
                u = mirt_dctn(x - eta);
                z = sign(u) .* max(abs(u) - g0 / rho, 0); 
                eta = eta - (x - mirt_idctn(z));

                if mod(n, stepsize) == 0
                    disp('saved');
                    xmat(:,:,:,n/stepsize) = x;
                end        
                                        
                toc;
            end
            
            %if isempty(xmat)
            %    xmat = x;
            %end            
        end
    end

    methods
    
        function h = computeHistogram(obj, lm, data)
            if (size(lm,1)~= 4 && size(lm,1)~= 5)
                error('unknown listmode data format, should be 4xN or 5xN (TOF) (currently not support DOI bins)');
            end
            
            M = obj.system_parms.crystal_array_size(1) * obj.system_parms.number_of_detector_modules_transaxial;
            nring = obj.getNumberOfCrystalRings;
        
            % remove unexpected prompts
            ii=((lm(1,:)<M) & (lm(1,:)>=0)) & ...
               ((lm(3,:)<M) & (lm(3,:)>=0)) & ...
               ((lm(2,:)<nring) & (lm(2,:)>=0)) & ...
               ((lm(4,:)<nring) & (lm(4,:)>=0));
               
            if length(ii) ~= size(lm, 2)
                warning('invalid listmode data found, # = %d', size(lm,2)-length(ii));
            end
            
            lm = lm(:,ii);
            idx = sub2ind([M, M, nring, nring], ...
                  double(lm(1,:))+1, double(lm(3,:))+1, double(lm(2,:))+1, double(lm(4,:))+1);
            h = zeros(M, M, nring, nring);
            
            for i=1:length(idx)
                h(idx(i)) = h(idx(i)) + data(i);
            end

        end






    	function h = l2h(obj, lm, flag)
    	%convert listmode data to histograms
    	%function h = l2h(lmdata, flag)
    		if nargin < 3
    			flag = 0;
    		end
    	
    		if (size(lm,1)~= 4 && size(lm,1)~= 5)
				error('unknown listmode data format, should be 4xN or 5xN (TOF) (currently not support DOI bins)');
			end
			
			M = obj.system_parms.crystal_array_size(1) * obj.system_parms.number_of_detector_modules_transaxial;
			nring = obj.getNumberOfCrystalRings;
		
			% remove unexpected prompts
			ii=((lm(1,:)<M) & (lm(1,:)>=0)) & ...
			   ((lm(3,:)<M) & (lm(3,:)>=0)) & ...
			   ((lm(2,:)<nring) & (lm(2,:)>=0)) & ...
			   ((lm(4,:)<nring) & (lm(4,:)>=0));
			   
			if length(ii) ~= size(lm, 2)
				warning('invalid listmode data found, # = %d', size(lm,2)-length(ii));
			end
			
			lm = lm(:,ii);
			idx = sub2ind([M, M, nring, nring], ...
				  double(lm(1,:))+1, double(lm(3,:))+1, double(lm(2,:))+1, double(lm(4,:))+1);
			h = zeros(M, M, nring, nring);
						
			fprintf('# of counts being histogramed = %d\n', size(lm, 2));
			 			
			if flag
			
				ss = sort(idx, 'ascend');
				nz=0;
				c=1;
				while (~isempty(ss))

					dd = diff(ss);
					iii = [1, find(dd>0)+1]; nz=nz + length(iii);
					h(ss(iii)) = h(ss(iii)) + 1;

					ss(iii)=-1;
					ss=ss(ss>0);
					if 0
					   fprintf('#%d: %d (%d)\n', c, length(iii), size(lm,2));
					end
					c=c+1;
	
				end
				fprintf('nz=%d\n', nz);
			
			else
				for n = 1 : length(idx)
					h(idx(n)) = h(idx(n)) + 1;
				end
			end
    	end
    	
    	function s = h2s(obj, h, flag)
    	%convert histogram to singoram (h: MxMxnringxnring)
    	%function s = h2s(h)
    		%fprintf('%d\n', nargin);
    		if nargin < 3
    			flag = false;
    		end
    		xp = obj.getDefaultSinogramCrystalPairs;
    		idx = sub2ind([size(h,1), size(h,2)], xp(1,:), xp(2,:));
    		s = zeros(size(xp,2), size(h,3), size(h,4));
    		for n = 1 : size(h,4)
				for m = 1 : size(h,3)
					if flag % alread joint!
						hh = h(:,:,m,n);
					else
						hh = h(:,:,m,n) + h(:,:,n,m)';
					end
					s(:, m, n) = hh(idx);
				end
			end
			s = reshape(s, obj.system_parms.number_of_projections_per_angle, ...
						obj.getDefaultNumberOfAngles, size(h,3), size(h,4));
			disp('NOTE: natural order!');
    	end







        %% convert to TOF histgram and TOF sinogram

        function h = l2h_tof(obj, lm, tof_nbin, tof_nbin_abs_max, flag)
        %convert TOF listmode data to TOF histograms
        %function h = l2h(lmdata, flag)
            if nargin < 3
                flag = 0;
            end
        
            if (size(lm,1)~= 4 && size(lm,1)~= 5)
                error('unknown listmode data format, should be 4xN or 5xN (TOF) (currently not support DOI bins)');
            end
            
            M = obj.system_parms.crystal_array_size(1) * obj.system_parms.number_of_detector_modules_transaxial;
            nring = obj.getNumberOfCrystalRings;
        
            % remove unexpected prompts
            ii=((lm(1,:)<M) & (lm(1,:)>=0)) & ...
               ((lm(3,:)<M) & (lm(3,:)>=0)) & ...
               ((lm(2,:)<nring) & (lm(2,:)>=0)) & ...
               ((lm(4,:)<nring) & (lm(4,:)>=0));
               
            if length(ii) ~= size(lm, 2)
                warning('invalid listmode data found, # = %d', size(lm,2)-length(ii));
            end
            
            % tof_nbin_abs_max = abs(max(lm(5,:)));
            % if abs(min(lm(5,:))) > tof_nbin_abs_max
            %     tof_nbin_abs_max = abs(min(lm(5,:)));
            % end

            lm = lm(:,ii);
            idx = sub2ind([M, M, nring, nring, tof_nbin], ...
                  double(lm(1,:))+1, double(lm(3,:))+1, double(lm(2,:))+1, double(lm(4,:))+1, double(lm(5,:)+tof_nbin_abs_max+1));
            h = zeros(M, M, nring, nring, tof_nbin);
                        
            fprintf('# of counts being histogramed = %d\n', size(lm, 2));
                        
            if flag
            
                ss = sort(idx, 'ascend');
                nz=0;
                c=1;
                while (~isempty(ss))

                    dd = diff(ss);
                    iii = [1, find(dd>0)+1]; nz=nz + length(iii);
                    h(ss(iii)) = h(ss(iii)) + 1;

                    ss(iii)=-1;
                    ss=ss(ss>0);
                    if 0
                       fprintf('#%d: %d (%d)\n', c, length(iii), size(lm,2));
                    end
                    c=c+1;
    
                end
                fprintf('nz=%d\n', nz);
            
            else
                for n = 1 : length(idx)
                    h(idx(n)) = h(idx(n)) + 1;
                end
            end
        end
        


        % function s = h2s_tof(obj, h, tof_nbin, tof_nbin_abs_max, flag)
        % %convert TOF histogram to TOF singoram (h: MxMxnringxnring)
        % %function s = h2s(h)
        %     %fprintf('%d\n', nargin);
        %     if nargin < 3
        %         flag = false;
        %     end
        %     xp = obj.getDefaultSinogramCrystalPairs;
        %     idx = sub2ind([size(h,1), size(h,2)], xp(1,:), xp(2,:));
        %     s = zeros(size(xp,2), size(h,3), size(h,4), tof_nbin);
        %     for n = 1 : size(h,4)
        %         for m = 1 : size(h,3)
        %             if flag % alread joint!
        %                 hh = h(:,:,m,n, :);
        %             else
        %                 hh = h(:,:,m,n, :) + h(:,:,n,m, :)';
        %             end
        %             s(:, m, n, :) = hh(idx);
        %         end
        %     end
        %     s = reshape(s, obj.system_parms.number_of_projections_per_angle, ...
        %                 obj.getDefaultNumberOfAngles, size(h,3), size(h,4), );
        %     disp('NOTE: natural order!');
        % end




    	
    	function mask = getHistogramMask(obj)
    		xp = obj.getDefaultSinogramCrystalPairs;
    		nx = obj.system_parms.number_of_detector_modules_transaxial * obj.system_parms.crystal_array_size(1);
    		mask = zeros(nx);
    		for n=1:size(xp,2)
    			mask(xp(1,n),xp(2,n))=1;
    			mask(xp(2,n),xp(1,n))=1;
    		end
    	end
    	
    end
    
end
