function [x,fval,xmat]=ez_lm_pcg(scanner, lmdata, image_size, voxel_size, maxstep, R, senimg, ci, ri, nite, x0, precond, chat, flag, step)

if nargin < 8
    chat = 1;    
end

%if isempty(R)
%    R=speye(prod(image_size(1:2)));
%end

% default normalization
if isempty(ci)
	ci=ones(size(lmdata,2),1);
end

% default randoms
if isempty(ri)
	ri=0;
end

angle_list = (0:(maxstep-1))/maxstep * 180.0;

fp = @(c)fproj_pawpet(c, image_size, voxel_size, ...
                      scanner.xtal_trans_location, ...
                      scanner.xtal_ring_offsets, ...
                      uint8(lmdata(:,:)), angle_list);

bp = @(xxx)bproj_pawpet(xxx, image_size, voxel_size, ...
                      	scanner.xtal_trans_location, ...
                        scanner.xtal_ring_offsets, ...
                        uint8(lmdata(:,:)), angle_list);

%
% f(x) = -senimg'*x + log(G*x + r)
% g = -senimg + G'*(1 ./ (G*x + r))
%

% initial image
if isempty(x0)
    x=double(senimg>0);
else
    x=x0(:);
end

%
cc = 1;
xmat = zeros(prod(image_size), nnz(mod([1:nite], step)==0));

xxx = reshape(x, image_size);
yyy = fp(xxx);

yp=ci(:).*(yyy(:)) + ri(:);
for n=1:nite
    
    if chat
        fprintf('----start iteration #%d----\n',n);
    end
    
    x_old = x;
    yp=max(yp,eps);   
    
    if flag 
        fval(n)=-senimg(:)'*x(:) + sum(log(yp));
    end
  
  	bbb=bp(ci(:)./yp);
  	g = bbb(:) - senimg(:);
%    g=(((ci(:)./yp)'*G)*R)' - senimg;    
    % preconditioner (apply your own preconditioner here!!!)
    
    if isempty(precond)
        m=(x + 1e-8); %./(max(senimg,eps));
        p=m(:).*g(:);
    else
        if size(precond,2)==0
            error('precond should be a matrix!');
        end
        p=precond * g;
    end
%
%    m = (G.^2)'*(1./(max(yp,eps)));
%
    
    if n==1
        d=p;
    else
%        gamma=((g-g_old)'*p)/(g_old'*p_old); 
        gamma=max(0,((g-g_old)'*p)/(g_old'*p_old));
        d=p+gamma*d;
        if chat
            fprintf('gamma=%f\n', gamma);
        end
    end
    g_old=g;
    p_old=p;

    % modify dir to maintain a continous update 
    % (both numeric precision and direction feasibility should be checked up!)
	d((x(:)<mean(x(:))*1e-6) & (d<0)) = 0;	
    
    % check conj grad 
    c = abs(d'*g) / (norm(d,2) * norm(g,2));
    if chat
        fprintf('c=%f\n',c);
    end
    if c < 0.0001
        d = g;
        % modify dir again
        d((x(:)<mean(x(:))*1e-6) & (d<0)) = 0;	
    end
        
    %
    ddd = d(:);
    yyy = fp(reshape(ddd, image_size));
    yd = ci(:).*yyy(:);
%    yd=ci(:).*(G*(R*d));
    if chat
        fprintf('===the 1st linesearch (exp amax=inf)===\n');
    end
    a=lsearch(yp,yd,senimg,d,inf,chat);
    if chat
        fprintf('=== 1st line search: a=%f ===\n\n', a);
    end
    x_new = x(:) + a * d(:);

    %
    %
    
    % trig the 2nd linesearch if negative values occur
    if nnz(x_new<0)>0
        x_new(x_new<0)=0;
        d=x_new(:) - x(:);
        ddd = d;
        yyy = fp(reshape(ddd, image_size));
        yd = ci(:).*yyy(:);
%        yd=ci(:).*(G*(R*d));
        if chat
            fprintf('===the 2nd linesearch (exp amax=1.0)===\n');
        end
        a=lsearch(yp,yd,senimg,d,1.0,chat);
        if chat
            fprintf('=== 2nd line search: a=%f ===\n\n', a);
        end
        x_new=x(:) + a*d(:);
    end
    
    % update
    x = x_new; %x + a*d;
    yp = yp+a*yd;

    if mod(n, step)==0
	    xmat(:, cc) = x(:);
	    cc = cc + 1;
	    fprintf('likelihood = %.12f\n', fval(n));
	end
        
    %
    if chat
%        fprintf('---- iteration #%d is done! likelihood=%f----\n\n',n, lk(n));
    end
end

%plot(1:nite,lk,'o-',1:nite,-lk_em); %imshow(sq(m));

function a=lsearch(yp,yd,senimg,d,amax,chat)
%senimg: sensitivity image
%d: dir
%yp: G*x
%yd: G*d
%x <= x+alpha*d
%
idx=yd<0;
if chat
    fprintf('determining maximum linesearch step ... ');
end
if nnz(idx)>0
%    yp(idx)./yd(idx)
    amax=min(-max(yp(idx)./yd(idx)), amax);
end
if chat
    fprintf('max a: %f\n', amax);
end

if chat
    fprintf('update step ... \n');
end
a=0;
for n=1:6
    ym=yp + a*yd;
    ym=max(ym, eps);
    lk1=-senimg(:)'*d(:) + yd'*(1./ym); %yd'*(y./ym-1);
    lk2=-yd'*((1.*yd) ./ ym.^2);
    a=a - lk1/lk2;
    if chat
    fprintf('lk1=%f,lk2=%f,a=%f\n', lk1, lk2,a);
    end
    if a>amax
        a=amax;
        if chat
            fprintf('break! ->%f\n', a);
        end
        break;
    end
end
