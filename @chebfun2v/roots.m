function varargout = roots( F, varargin )
%ROOTS   Find the common zeros of a CHEBFUN2V object.
%   r = ROOTS(F) finds the common zeros of the two bivariate functions F(1) and
%   F(2) in their domain of definition under the assumption that the solution
%   set is zero-dimensional. R is a matrix with two columns storing the x- and
%   y-values of the solutions. This script is also called by the syntax
%   ROOTS(f,g), where f and g are CHEBFUN2 objects.
%
%   [x, y] = ROOTS(F) returns the x- and y-values as two separate columns.
%
%   Currently, if the maximum degree of F(1) and F(2) is greater than 200 then
%   an algorithm based on Marching squares is employed, and an algorithm based
%   on a resultant method is used otherwise (see [1]).
%
%   ROOTS(F, 'ms') or ROOTS(F, 'marchingsquares') always employs the marching
%   squares algorithm.
%
%   ROOTS(F, 'resultant') always employs the algorithm based on the hidden
%   variable resultant method.
%
%   [1] Y. Nakatsukasa, V. Noferini, and A. Townsend, Computing the common zeros
%   of two bivariate functions via Bezout resultants, (2013).
%
% See also CHEBFUN2/ROOTS, CHEBFUN/ROOTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Maximum degree for resultant method:
max_degree = 200;  

% Empty check:
if ( isempty(F) )
    varargout = {[]};
    return
end

f = F.components{1}; 
g = F.components{2};
[nf, mf] = length(f);
[ng, mg] = length(g);
% Maximum degree:
dd = max([mf, nf, mg, ng]); 

validArgs = {'ms', 'marchingsquares', 'resultant'};
if ( (nargin > 1) && ~any(strcmpi(varargin{1}, validArgs)) )
    error('CHEBFUN:CHEBFUN2V:roots:badInput', ...
        'Unrecognised optional argument.');
end


if ( isempty(varargin) )
    % If rootfinding method has not been defined, then use resultant if
    % degrees are small: 
    if ( dd <= max_degree ) 
        [xroots, yroots] = roots_resultant(F);
    else
        [xroots, yroots] = roots_marchingSquares(F);
        xroots = xroots.'; 
        yroots = yroots.';
    end
elseif ( strcmpi(varargin{1}, 'resultant') )
    % If the user wants the resultant method, then use it: 
    [xroots, yroots] = roots_resultant(F);
elseif ( any( strcmpi(varargin{1}, {'ms', 'marchingsquares'} ) ) )
    % If the user wants the marching squares method, then use it: 
    [xroots, yroots] = roots_marchingSquares(F);
    xroots = xroots.'; 
    yroots = yroots.';
else
    % Print error. Unknown method. 
    error('CHEBFUN2V:ROOTS:METHOD', 'Unknown rootfinding method.') 
end

if ( nargout <= 1 )
    varargout{1} = [xroots ; yroots].';
else
    varargout = {xroots, yroots};
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%      RESULTANT METHOD        %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xroots,yroots] = roots_resultant(F)

% extract out the two CHEBFUN2 objects.
f = F.components{1}; g = F.components{2};

% Useful parameters
rect = f.domain;  % rectangular domain of CHEBFUN2.
[nf,mf]=length(f);[ng,mg]=length(g);
dd = max([mf nf mg ng]); % max degree
max_degree = min( 16 , dd ); % subdivision threshold degree
reg_tol = 1e-15;   % Regularization (for Bezoutian)
domain_overlook = 1e-10; % start by looking for zeros on larger domain.
donewton = 0;   % Newton polishing?
multitol = sqrt(eps); % for multiple roots we do not go for more than sqrt(eps);

% Consider a slightly large domain to make sure we get all the roots along
% edges.
% domain width (attempt to make scale invariant)
xwid = diff(rect(1:2))/2; ywid = diff(rect(3:4))/2;
xmax = rect(2) + xwid*domain_overlook;
xmin = rect(1) - xwid*domain_overlook;
ymax = rect(4) + ywid*domain_overlook;
ymin = rect(3) - ywid*domain_overlook;

subdividestop = 2*(1/2)^((log(16)-log(dd))/log(.79)); % subdivision threshold
subdividestop = min(subdividestop,1/4);

% initial scaling to O(1)
f = f/abs(f.pivotValues(1));
g = g/abs(g.pivotValues(1));

[xroots,yroots] = subrootsreptwo(f,g,xmin,xmax,ymin,ymax,xwid,ywid,reg_tol,max_degree,subdividestop);
% ballpark obtained, next do local bezoutian refinement


[fx, fy] = gradient(chebfun2(f,rect)); 
[gx, gy] = gradient(chebfun2(g,rect));

%%%%% subdivision for accuracy when dynamical range is an issue %%%%%%%%
% find region in which roots might have been missed
F = chebcoeffs2(f); 
G = chebcoeffs2(g);
F = rot90(F, -2);
G = rot90(G, -2);

xpts = linspace(xmin,xmax,2*max(size(F,2),size(G,2)));
ypts = linspace(ymin,ymax,2*max(size(F,1),size(G,1)));

[xx,yy] = meshgrid(xpts,ypts);
FX = fx(xx,yy);FY = fy(xx,yy); GY = gy(xx,yy);GX = gx(xx,yy);
FG = abs(FX.*GY-FY.*GX);  % Bezout conditioning
FpG = abs(f(xx,yy))+abs(g(xx,yy)); % function values

FGG = (log10(FpG+eps)<-6 & log10(FG/max(max(FG))+eps)<-12); % dangerous region (Bezout might have missed)
if sum(sum(FGG))>0, % Bezout could have missed these, do local search
    xx = zeros(1,size(FGG,1)*size(FGG,2)); yy = xx;
    ip = 1;
    for i=1:size(FGG,1)
        for j=1:size(FGG,2)
            if FGG(i,j)==1,
                s = svd([FX(i,j) FY(i,j);GX(i,j) GY(i,j)]);
                Jacobs = 1./s(end); % conditioning of original problem
                if ( Jacobs < 1e15 ), % give up if solution is too ill conditioned
                    xx(ip) = xpts(j);        yy(ip) = ypts(i);
                    ip = ip+1;
                end
            end
        end
    end
    
    xdis = (xpts(2)-xpts(1)); ydis = (ypts(2)-ypts(1));
    xx = xx(1:ip-1); yy = yy(1:ip-1);
    m = blockingxy(xx,yy,max(xdis,ydis)*1.1);
    for i = 1:max(m)
        xnow = xx(m==i); ynow = yy(m==i);
        [xt,yt] = subrootsreptwo(f,g,min(xnow)-xdis,max(xnow)+xdis,min(ynow)-ydis,max(ynow)+ydis,xwid,ywid,reg_tol,max_degree); % local (second) bez
        xroots = [xroots xt];  yroots = [yroots yt];
    end
end


% xroots,yroots should contain solutions (condnum<1e14) and maybe more.
tolblock = xwid*min(3e-3,1/(10*dd)); % clustering threshold
m = blockingxy(xroots,yroots,tolblock);
z = xroots+1i*yroots;
xroots = [];yroots = [];

tolbdefault = min(5*xwid*sqrt(multitol),tolblock); % local bezout width
for i = 1:max(m)
    znow = z(m==i); xnow = real(znow); ynow = imag(znow);
    J = [fx(mean(xnow),mean(ynow)) fy(mean(xnow),mean(ynow));gx(mean(xnow),mean(ynow)) gy(mean(xnow),mean(ynow))];
    s = svd(J);
    Jacobs = 1./s(end);     % original conditioning
    Jacobsbez = 1/abs(det(J)); % Bezout conditioning
    tolb = max(tolbdefault,eps*100*Jacobsbez);
    tolb = min(tolb,xwid*1e-2); % give up if Jacobs too large
    [xrootsnow,yrootsnow] = subrootsreptwo(f,g,min(xnow)-tolb,max(xnow)+tolb,min(ynow)-tolb,max(ynow)+tolb,xwid,ywid,0,max_degree); % local (second) bez no subdivision
       
    xxx = []; yyy = [];
    if length(xrootsnow)>=1,
        
        if Jacobs<1e10, %only for well-conditioned roots
            tolbb = max(xwid*multitol,Jacobs*1000*eps);
            for j = 1:length(xrootsnow)
                [xx,yy] = subrootsreptwo(f,g,min(xrootsnow(j))-tolbb,max(xrootsnow(j))+tolbb,min(yrootsnow(j))-tolbb,max(yrootsnow(j))+tolbb,xwid,ywid,0,max_degree/2); % local (second) bez
                xxx = [xxx xx]; yyy = [yyy yy];
            end
            
        else
            xxx = xrootsnow;
            yyy = yrootsnow;
        end
        
    end
    xroots = [xroots xxx]; yroots = [yroots yyy];
end

% finally collapse close roots
xrootsnow = [];    yrootsnow = [];
mm = blockingxy(xroots,yroots,min(xwid*sqrt(eps)*10));
for ii = 1:max(mm)
    xtmp = xroots(mm==ii);  ytmp = yroots(mm==ii);
    res = zeros(1,length(xtmp));
    for ij = 1:length(xtmp)
        res(ij) = norm([feval(f,xtmp(ij),ytmp(ij)) feval(g,xtmp(ij),ytmp(ij))]);
    end
    [mres,IX] = min(res);
    if min(mres)<1e-10,
        xrootsnow = [xrootsnow xtmp(IX)];  yrootsnow = [yrootsnow ytmp(IX)];
    end
end
xroots = xrootsnow; yroots = yrootsnow;

if donewton % optional Newton update (don't do by default)
    [xroots,yroots] = newtonupdate(f,g,xroots,yroots);
end


% discard roots safely outside initial domain
xmax = xmax - xwid*domain_overlook; xmin = xmin + xwid*domain_overlook; % original domain
ymax = ymax - ywid*domain_overlook; ymin = ymin + ywid*domain_overlook;

ii = find(xroots<xmax+xwid*1e-15 & xroots>xmin-xwid*1e-15 & yroots<ymax+ywid*1e-15 & yroots>ymin-ywid*1e-15);
xroots = xroots(ii); yroots = yroots(ii);
for i=1:length(xroots)
    if xroots(i)>xmax, xroots(i) = xmax; end % push outliers to boundary
    if xroots(i)<xmin, xroots(i) = xmin; end
    if yroots(i)>ymax, yroots(i) = ymax; end
    if yroots(i)<ymin, yroots(i) = ymin; end
end

end



function [xroots,yroots] = subrootsreptwo(f,g,xmin,xmax,ymin,ymax,xwid,ywid,tolreg,maxd,subdividestop)
% execute subdivision and if degrees small enough, run bezout
xroots = [];yroots = [];
if xmin==xmax || ymin==ymax, return , end % empty domain

approx_tol = 1e-13;  % Approximation accuracy.
magicnum = 0.004849834917525; % 'magic number'
honournum = -0.0005194318842611; % 'honourary number'
if exist('subdividestop','var')==0, subdividestop = 0.008;end

ff = cheb2(f,[xmin xmax ymin ymax],approx_tol,maxd+2); % sample at a bit more than maxd points
gg = cheb2(g,[xmin xmax ymin ymax],approx_tol,maxd+2);
fcoef = ff.coeffs; gcoef = gg.coeffs;

% subdivision test
if (xmax-xmin)>xwid*subdividestop || (ymax-ymin)>ywid*subdividestop % subdivide only if domain is not too small
    if (size(fcoef,1)>maxd && size(fcoef,2)>maxd) ||...
            (size(gcoef,1)>maxd && size(gcoef,2)>maxd)
        % subdivide in both x,y
        xmed = (xmax+xmin)/2 - magicnum * (xmax-xmin)/2;
        ymed = (ymin+ymax)/2 - honournum * (ymax-ymin)/2;
        [xroots1,yroots1] = subrootsreptwo(f,g,xmin,xmed,ymin,ymed,xwid,ywid,tolreg,maxd,subdividestop);
        [xroots12,yroots12] = subrootsreptwo(f,g,xmed,xmax,ymin,ymed,xwid,ywid,tolreg,maxd,subdividestop);
        [xroots2,yroots2] = subrootsreptwo(f,g,xmin,xmed,ymed,ymax,xwid,ywid,tolreg,maxd,subdividestop);
        [xroots22,yroots22] = subrootsreptwo(f,g,xmed,xmax,ymed,ymax,xwid,ywid,tolreg,maxd,subdividestop);
        xroots = [xroots xroots1 xroots2 xroots12 xroots22];
        yroots = [yroots yroots1 yroots2 yroots12 yroots22];
        return
    elseif (size(fcoef,1)>maxd)||(size(gcoef,1)>maxd) % subdivide in y
        ymed = (ymin+ymax)/2 - magicnum * (ymax-ymin)/2;
        [xroots1,yroots1] = subrootsreptwo(f,g,xmin,xmax,ymin,ymed,xwid,ywid,tolreg,maxd,subdividestop);
        [xroots12,yroots12] = subrootsreptwo(f,g,xmin,xmax,ymed,ymax,xwid,ywid,tolreg,maxd,subdividestop);
        xroots = [xroots xroots1 xroots12];
        yroots = [yroots yroots1 yroots12];
        return
    elseif (size(fcoef,2)>maxd)||(size(gcoef,2)>maxd) % subdivide in x
        xmed = (xmax+xmin)/2 - honournum * (xmax-xmin)/2;
        [xroots1,yroots1] = subrootsreptwo(f,g,xmin,xmed,ymin,ymax,xwid,ywid,tolreg,maxd,subdividestop);
        [xroots12,yroots12] = subrootsreptwo(f,g,xmed,xmax,ymin,ymax,xwid,ywid,tolreg,maxd,subdividestop);
        xroots = [xroots xroots1 xroots12];
        yroots = [yroots yroots1 yroots12];
        return
    end
else
    ff = cheb2(f,[xmin xmax ymin ymax],approx_tol,round(1.5*maxd)); % didn't resolve, sample at more than maxd points but not too much
    gg = cheb2(g,[xmin xmax ymin ymax],approx_tol,round(1.5*maxd));
end

F = ff.coeffs; G = gg.coeffs;
if isempty(F); F = 0; end
if isempty(G); G = 0; end
if (2-eps*10)*abs(F(end,end))>sum(sum(abs(F))) || (2-eps*10)*abs(G(end,end))>sum(sum(abs(G)))
    % 'no roots here!'
else
    if min(min(size(F)),min(size(G)))<=1,
        if length(F)<=1 && length(G)<=1, % constant
            if (abs(F)<=1e-15) && (abs(G)<=1e-15), % flat 0; just say middle point is 0
                xroots = (xmax+xmin)/2; yroots = (ymax+ymin)/2;
            else
                xroots = [];yroots=[];
            end
        else
            %    'no roots here!'
            [xroots,yroots] = onevar(F,G,xmin,xmax,ymin,ymax);
        end
    else
        
        [xroots,yroots] = runbezval(F,G,xmin,xmax,ymin,ymax,tolreg);
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xroots,yroots] = runbezval(F,G,xmin,xmax,ymin,ymax,tol)
% run bezout (core code):
% form Bezout matrix polynomial and solve eigenproblem, then find other variable
% form Bezoutian at magic number and find null space and diagonal balancing
% regularization tolerance

mc = size(F,1)-1+size(G,1)-1+1; % # of sampling points (degree of B(y))
mc = max(mc,2);

magicnum = 0.004849834917525; % 'magic number'

doswap = 0;

% swap x and y if appropriate
if max(size(F,1),size(G,1))*(size(F,2)+size(G,2)) > max(size(F,2),size(G,2))*(size(F,1)+size(G,1))
    F = F'; G = G';     doswap = 1;
    mc = size(F,1)-1+size(G,1)-1+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% form B and regularize %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B0 = formBez(F,G,magicnum); % sample Bezoutian at 'magic number'

% now form Bezoutians at Chebyshev points
x = chebpts(mc);
B = zeros(size(B0,1),size(B0,1),mc);
Bsum = zeros(size(B0));
for i = 1:mc
    B(:,:,i) = formBez(F,G,x(i));
    Bsum = Bsum+abs(B(:,:,i));
    %        B(:,:,i) = Btmp(k+1:end,k+1:end); % extract regular part
end

if tol == 0,  % no regularization
    k = 0;
else
    %%%%%%%%%% REGULARIZATION parameter %%%%%%%%%%%%%%% extract lower-right part of Bezoutian %%%
    d = diag(B0);
    for i = 1:size(B0,1),   d(i) = max(max(abs(Bsum(1:i,1:i))));    end
    m = find(abs(d)/max(abs(d))> tol);
    if isempty(m) 
        k = 0;
    else
        offd = zeros(1,m(1)-1);
        for i = 1:m(1)-1
            offd(i) = max(max(abs(Bsum(i+1:end,1:i))));
        end
        mm = find(offd/max(abs(d))<sqrt(tol));
        if isempty(mm), m(1) = 1; else m(1) = max(mm); end
        k = max(m(1)-1,0);
    end
end
%regsize = [length(B0) k]

if length(size(B))<3,
    ei = [];
else
    B = B(k+1:end,k+1:end,:); % regularize
    
    ns = size(B);   BB = reshape(B,ns(1),ns(1)*ns(3));
    CC = matrixChebfft(BB);
    for i = 1:mc
        B(:,:,i) = CC(:,(i-1)*size(B,1)+1:i*size(B,1));
    end
    
    % cutoff negligible B
    nrmB = norm(B(:,:,end),'fro');
    for ii=1:size(B,3)
        if norm(B(:,:,ii),'fro')/nrmB > 10*eps,    break;    end
    end
    B = B(:,:,ii:end);
    ns = size(B);
    
    %Bori = B;
    % diagonal balancing (optional)
    [Dori] = balancecong(B0(k+1:end,k+1:end));
    for i = 1:size(B,3)
        B(:,:,i) = Dori * B(:,:,i) * Dori;
        B(:,:,i) = rot90(B(:,:,i),2); % fliplr,ud
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% done forming B, now solve polyeig %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % colleague QZ (most stable but slow)
    nrm = norm(BB,'fro')/ns(1); % scale as suggested as Van Dooren
    ei = chebT1rtsmatgep(B/nrm);
    %[ei,AG,BG] = chebT1rtsmatgep(B/nrm);    ei = stairqz(AG,BG);     %    semi-staircase
    
end
yreal = sort(real(ei(abs(real(ei))<=1+10*eps & abs(imag(ei))<sqrt(eps)*10)));

% finally obtain x-values
[xroots,yroots] = xunivariate(F,G,xmin,xmax,ymin,ymax,yreal,doswap);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   UNIVARIATE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xval,yval] = xunivariate(F,G,xmin,xmax,ymin,ymax,y,doswap)
% compute x-vals from y
% collapse if near-multiple y-values are present

magicnum = 0.004849834917525; % 'magic number'
ep = 10*eps; % tolerance for collapsing
dist = abs(y(1:end)-[y(2:end);2]);
ii = dist>ep;
y = y(ii);


xreal = magicnum*ones(size(y));
xtmp = [];ytmp = [];
for j = 1:length(y)
    if xreal(j)==magicnum
        
        % compute coefficients via vandermonde vector
        vy = vandercheb(y(j),size(G,1));     coef1 = vy'*G;
        vy = vandercheb(y(j),size(F,1));     coef2 = vy'*F;
        
        % find roots of two univariate polys via colleague
        if length(coef1)<=1 || norm(coef1)==0, rts1=[];
        else
            rts1 = chebT1rts(coef1(find(abs(coef1)>0,1,'first'):end)');
            rts1 = real(sort(rts1(abs(imag(rts1))<=10*1e-8 & abs(rts1)<=1+10*eps)));
        end
        if length(coef2)<=1 || norm(coef2)==0, rts2=[];
        else
            rts2 = chebT1rts(coef2(find(abs(coef2)>0,1,'first'):end)');
            rts2 = real(sort(rts2(abs(imag(rts2))<=10*1e-8 & abs(rts2)<=1+10*eps)));
        end
        % collapse nearby roots if any (deal with multiple roots)
        dist = abs(rts1(1:end)-[rts1(2:end);2]);
        ii = dist>10*ep;  rts1 = rts1(ii);
        dist = abs(rts2(1:end)-[rts2(2:end);2]);
        ii = dist>10*ep;  rts2 = rts2(ii);
        rts = sort([rts1;rts2]);    res = ones(size(rts));
        
        % compute function values
        for i = 1:length(rts)
            res(i) = max([abs(coef1*vandercheb(rts(i),size(G,2))),abs(coef2*vandercheb(rts(i),size(F,2)))]);
        end
        
        % adopt depending on the interval size
        if xmax-xmin > 2e-3,
            %rtscan = rts(find(res<1000*sqrt(eps))); % candidates for x-values
            rtscan = rts((res<100*sqrt(eps))); % candidates for x-values
        elseif xmax-xmin > 1e-7, % local bez
            rtscan = rts((res<1e-10));
        else                     % very local (final) bez
            rtscan = rts((res<1e-12));
        end
        
        % collect/collapse candidate roots
        if ~isempty(rtscan)
            skipnext = 0;
            for i = 1:length(rtscan)-1
                if skipnext
                    skipnext = 0;
                else
                    if ( abs(rtscan(i)-rtscan(i+1))<10*ep ),
                        skipnext = 1;
                        vali = max([abs(coef1*vandercheb(rtscan(i),size(G,2))),  abs(coef2*vandercheb(rtscan(i),size(F,2)))]);
                        valip= max([abs(coef1*vandercheb(rtscan(i+1),size(G,2))),abs(coef2*vandercheb(rtscan(i+1),size(F,2)))]);
                        if ( vali < valip )
                            xtmp = [xtmp;rtscan(i)];
                        else
                            xtmp = [xtmp;rtscan(i+1)];
                        end
                    else
                        xtmp = [xtmp;rtscan(i)];
                    end
                    ytmp = [ytmp;y(j)];
                end
            end
            % last i
            if skipnext % do nothing
            else
                xtmp = [xtmp;rtscan(end)];     ytmp = [ytmp;y(j)];
            end
        end
    end
end

if doswap ==0 % recover if swapped initially
    xval  = (xmin+xmax)/2+(xmax-xmin)/2*xtmp';
    yval  = (ymin+ymax)/2+(ymax-ymin)/2*ytmp';
else
    xval  = (xmin+xmax)/2+(xmax-xmin)/2*ytmp';
    yval  = (ymin+ymax)/2+(ymax-ymin)/2*xtmp';
end

return
end

%%%%%%%%%%%%%%%%%%%%%% KEEP IN CASE I NEED AT SOME STAGE %%%%%%%%%%%%%%%%% 
% function [xroots,yroots,sproots] = localrefine(f,g,xmin,xmax,ymin,ymax,xwid,ywid,tolreg,maxd)
% % execute subdivision and if degrees small enough, run bezout
%
% approx_tol = 1e-13; % cheb2 cutoff tolerance
% ff = cheb2(f,[xmin xmax ymin ymax],approx_tol,maxd+2); % sample at a bit more than maxd points
% gg = cheb2(g,[xmin xmax ymin ymax],approx_tol,maxd+2);
% F = ff.coeffs; G = gg.coeffs;
%
%     if isempty(F),
%             v = vandercheb(0,size(G,2));
%             G = G*v;
%             rx = chebT1rts(G(find(abs(G)>0,1,'first'):end)')';
%             rx = real(sort(rx(abs(imag(rx))<=1e-8 & abs(rx)<=1+10*eps)));
%             yroots  = (ymin+ymax)/2+(ymax-ymin)/2*rx;
%             xroots = (xmin+xmax)/2*ones(size(yroots));
%         return,
%     end
%     if isempty(G),
%             v = vandercheb(0,size(F,2));
%             F = F*v;
%             rx = chebT1rts(F(find(abs(F)>0,1,'first'):end)')';
%             rx = real(sort(rx(abs(imag(rx))<=1e-8 & abs(rx)<=1+10*eps)));
%             yroots  = (ymin+ymax)/2+(ymax-ymin)/2*rx;
%             xroots = (xmin+xmax)/2*ones(size(yroots));
%         return,
%     end
%
%
%     if min(min(size(F)),min(size(G)))<=1,
%         if length(F)<=1 && length(G)<=1, % constant
%             if (abs(F)<=1e-15) && (abs(G)<=1e-15), % flat 0; just say middle point is 0
%                 xroots = (xmax+xmin)/2; yroots = (ymax+ymin)/2;
%             else
%                 xroots = [];yroots=[];
%             end
%         else
%             %    'no roots here!'
%             %[xmin xmax ymin ymax],keyboard
%             [xroots,yroots] = onevar(F,G,xmin,xmax,ymin,ymax);
%         end
%     else
%    [xroots,yroots] = runbezval(F,G,xmin,xmax,ymin,ymax,tolreg);        % bez not syl
%     end
% end


function [xroots,yroots] = onevar(F,G,xmin,xmax,ymin,ymax)
% when either F or G is constant
doswap = 0;

if min(min(size(F)))==1,
    if size(F,2) ==1, FF = F'; GG = G'; doswap = 1;
    else FF = F; GG = G;
    end
else
    if size(G,2) ==1, FF = G'; GG = F'; doswap = 1;
    else FF = G; GG = F;
    end
end

if length(FF)==1
    if abs(FF)<1e-15,
        if min(size(GG))==1,
            rx = chebT1rts(GG(find(abs(GG)>0,1,'first'):end)')';
            rx = real(sort(rx(abs(imag(rx))<=1e-8 & abs(rx)<=1+10*eps)));
            rr = rx;
        else % GG is matrix. just take x=0 and find y.
            %gcheb = chebfun2(g,[xmin xmax ymin ymax]);                roots(gcheb);
            v = vandercheb(0,size(GG,2));
            GG = GG*v;
            rx = chebT1rts(GG(find(abs(GG)>0,1,'first'):end)')';
            rx = real(sort(rx(abs(imag(rx))<=1e-8 & abs(rx)<=1+10*eps)));
            rr = rx;
        end
        if size(GG,1)>size(GG,2)
            yroots =rr; xroots = zeros(size(rr));
        else
            xroots =rr; yroots = zeros(size(rr));
        end
    else
        xroots = []; yroots = [];
    end
else
    rx = chebT1rts(FF(find(abs(FF)>0,1,'first'):end)');
    rx = real(sort(rx(abs(imag(rx))<=1e-8 & abs(rx)<=1+10*eps)));
    xroots = []; yroots = [];
    for j=1:length(rx)
        vy = vandercheb(rx(j),size(GG,2)); coef = GG*vy;
        % ry = roots(chebfun(coef,'coeffs'))';
        if length(coef)<=1,
            if norm(coef)<1e-15,
                ry = 0;
            else ry=[];
            end
        else
            rts2 = chebT1rts(coef(find(abs(coef)>0,1,'first'):end));
            rts2 = real(sort(rts2(abs(imag(rts2))<=1e-8 & abs(rts2)<=1+10*eps)));
            ry = rts2';
        end
        yroots = [yroots ry];
        xroots = [xroots rx(j)*ones(size(ry))];
    end
end
if doswap
    xt = xroots;  xroots = yroots; yroots = xt;
end
xroots  = (xmin+xmax)/2+(xmax-xmin)/2*xroots;
yroots  = (ymin+ymax)/2+(ymax-ymin)/2*yroots;
end


function [r,A,B] = chebT1rtsmatgep(c)
% CHEBT1RTSMATGEP(c), finds via QZ the roots of a MATRIX polynomial
% expressed in a ChebT basis
% by using the colleague matrix *pencil* of the first kind.  F can be a vector of
% coefficients or a chebfun. Coefficients are ordered  highest degree down.
%
% c is a nxnxk array.  k-1 is the degree, n is the matrix size.

k=length(c(1,1,:)); n=length(c(:,:,1));

if size(c,3) ==2,  % linear case
    r = eig(c(:,:,2),-c(:,:,1));
    return
end

for ii=2:k
    c(:,:,ii)=c(:,:,ii)*(-.5); % coefficients
end
c(:,:,3) = c(:,:,3)+.5*c(:,:,1);

oh = .5*ones(n*(k-2),1);
% form colleague matrix A,B:
A = diag(oh,n)+diag(oh,-n);
A(end-n+1:end,end-2*n+1:end-n) = eye(n);

for ii=1:k-1
    A(1:n,(ii-1)*n+1:ii*n) = c(:,:,ii+1);
end
B=eye(size(A)); B(1:n,1:n)=c(:,:,1);

r = eig(A,B);% Compute roots

end

function r = chebT1rts(c)
% CHEBT1RTS(F), finds the roots of a polynomial expressed in a Cheb T basis
% by using the colleague matrix of the first kind.  F can be a vector of
% coefficients or a chebfun.

if length(c)<=1 
    r=[]; 
    return
end
if size(c,1) < size(c,2) % c is column vector
    c = c';
end 

if length(c)==2, % linear case
    r = -c(2)/c(1);  return
end

if c(1)==0
    c(1)=eps*max(c);  % perturb a zero leading coefficient by eps. 
end
c = -.5*c(end:-1:2)/c(1); 
c(end-1) = c(end-1)+.5;
oh = .5*ones(length(c)-1,1);
% Modified colleague matrix:
A = diag(oh,1)+diag(oh,-1);
A(end,end-1) = 1; A(1,:) = flipud(c);
r = eig(A);% Compute roots as eig(A)

end

function B = formBez(F,G,y)
% forms Bezoutian using DLP for F,G (matrices) at a specific value of y.

yv = vandercheb(y,max(size(F,1),size(G,1))); % Chebyshev-vandermonde form vector

ff = yv(end-size(F,1)+1:end)'*F;
gg = yv(end-size(G,1)+1:end)'*G;

if length(ff)<length(gg), gt = gg; gg=ff; ff = gt;end
if length(ff)==length(gg), ff = [0 ff]; shrink = 1; else shrink = 0; end

B = DLPforbez(ff,[zeros(1,length(ff)-length(gg)-1)';gg']); %B=(B+B')/2; input k

if shrink 
    B = B(2:end,2:end);
end

end

function [Y] = DLPforbez(AA,v)
% DLPforbez constructs the Bezoutian. Highly specified for Chebyshev
% biroots.
%
[n m] = size(AA); k=m/n-1; s=n*k;              % matrix size and degree
S = [zeros(1,k+1);(2*v)*AA];
R = S'-S;
% Bartel-Stewart algorithm on M'Y+YM=R, M is upper triangular.
Y = zeros(s);
if s ==1, Y = R(1,2); return, end
Y(1,:) = R(1,2:end);
Y(2,:) = R(2,2:end)+[Y(1,2:end-1) 2*Y(1,end) 0]+[0 Y(1,1:end-1)];
for i = 3:k                                    % backwards substitution
    Y(i,:) = R(i,2:end)-Y(i-2,1:end)+[Y(i-1,2:end-1) 2*Y(i-1,end) 0]+[0 Y(i-1,1:end-1)];
end
Y(k,:) = Y(k,:)/2;
end


function D = matrixChebfft(A)
% First attempt and matrix chebfft. Given a set of matrix coefficients,
% this function is designed to return the set of matrix values.
% Assumption: The matrix coefficients are square.
n = size(A,1); k = size(A,2)/size(A,1);  % get matrix size and degree.

if ( abs( k - round(k) ) > 0 )
    error('CHEBFUN:CHEBFUN2V:roots:matrixChebfft:badDegree', ...
        'Degree must be integer');
end

D = A;
for jj = 1:n  % for each column of A
    B = A(:,jj:n:n*k);
    C = chebtech2.vals2coeffs(B.');   % convert first column of each coefficient to values.
    D(:,jj:n:n*k) = rot90(C, -1);     % assign to output.
end

end

function m = blockingxy(x,y,delta)
% BLOCKINGXY  Produce blocking pattern for data x,y within tolerance delta.
% Elements will be assigned numbers for each class.

if isempty(x), m=[]; return, end
[xx,IX] = sort(x); m = ones(size(x));

xxp = ones(size(x)); ppos = ones(size(x));
p = 1;
for i=1:length(x)-1
    if abs(xx(i)-xx(i+1))>delta,
        p = p+1;    ppos(p) = i+1;
    end
end
ppos = ppos(1:p); % done with x
q = 0;
for i=1:p
    q = q+1;
    if i<p
        m(IX(ppos(i):ppos(i+1)-1)) = q*ones(size(ppos(i):ppos(i+1)-1));
        ynow = y(IX(ppos(i):ppos(i+1)-1));
    else % i=p, end
        m(IX(ppos(i):end)) = q*ones(size(IX(ppos(i):end)));
        ynow = y(IX(ppos(i):end));
    end
    [yy,IY] = sort(ynow);
    for j = 1:length(yy)-1
        if abs(yy(j)-yy(j+1))>delta,
            q = q+1;
            m(IX(ppos(i)+IY(j+1:end)-1)) = q*ones(size(IX(ppos(i)+IY(j+1:end)-1)));
        end
    end
end
end

function v = vandercheb(x,n)
% v = vandercheb(x,n) forms unit vector in vandermonde form
v = zeros(n,1);
for i = 1:n 
    v(end-i+1) = real(cos((i-1)*acos(x)));
end
end

function D = balancecong(B)
% diagonal congruence balancing to make the diagonals of DBD equal.
D = eye(length(B));
for i = length(B)-1:-1:1;
    if norm(B(i,:))>0,
        D(i,i) = max(1,sqrt(norm(B(end,:))/norm(B(i,:))));
    else
        D(i,i) = D(i+1,i+1);
    end
end
end

function [xroots,yroots] = newtonupdate(f,g,xroots,yroots,itnum)
%%%%%%%%%%%% newton update %%%%%%%%%%%%%%%%
% Use one iteration of Newton to get 14-15 digits.
if nargin < 5
    itnum = 1;
end
f = chebfun2(f); g = chebfun2(g) ;

tol = 1e-15;
fx = diff(f,1,2); fy=diff(f); gx = diff(g,1,2); gy=diff(g);            % derivatives.
J = @(x,y) [feval(fx,x,y) feval(fy,x,y);feval(gx,x,y) feval(gy,x,y)];  % Jacobian
for jj=1:itnum
    r = [xroots' yroots'];
    for kk = 1:size(r,1)
        x0 = [r(kk,1),r(kk,2)].';dx=1; iter = 1;
        while ( norm(dx) > 10*tol && iter < 2 )
            dx = J(x0(1),x0(2)) \ -[feval(f,x0(1),x0(2));feval(g,x0(1),x0(2))];    % update
            x0 = dx + x0; iter = iter + 1;
        end
        r(kk,:) = x0;
    end
    xroots = r(:,1)'; yroots = r(:,2)';
end

end

function g = cheb2(f,varargin)
% Basic bivariate tensor product, nonadaptive constructor.

default_tol = 10*eps;
default_n = 300;

% Did we get any user defined domain?
if nargin == 2
    % Assume this is a user defined domain [a b c d]
    if numel(varargin{1}) > 1
        ends = varargin{1};
        tol = default_tol;
    elseif numel(varargin{1}) == 1
        tol = varargin{1};
    end
elseif nargin == 3
    ends = varargin{1};
    tol = varargin{2};
elseif nargin == 4
    ends = varargin{1};
    tol = varargin{2};
    n = varargin{3};
else
    % default to the unit interval [-1 1 -1 1]
    ends = [-1 1 -1 1];
    tol = default_tol;
    n = default_n;
end

if isstruct(f)
    g = f;
    return;
end

% evaluate the function on a grid.
if exist('n','var')==0,
    n = 300;
end

x = mypoints(n,ends(1:2)); y = mypoints(n,ends(3:4));
[xx, yy]=meshgrid(x,y); F = f(xx,yy);

% vertical scale for machine precision
vscl = max(1,max(abs(F(:))));  % don't go for more than absolute accuracy.

% Compute bivariate Chebyshev T coefficients.
C = chebfun2.vals2coeffs(F);
C = rot90(C, -2);

% Very simple truncation of the coefficients.
%m = find(max(abs(C))>100*eps*vscl,1,'first'); n = find(max(abs(C.'))>100*eps*vscl,1,'first');
m = find(max(abs(C))>tol*vscl,1,'first'); n = find(max(abs(C.'))>tol*vscl,1,'first');
C = C(n:end,m:end);

% Form cheb2 object.
g = struct('coeffs',C,'scl',vscl,'corners',ends);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%       Marching Squares       %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xroots, yroots] = roots_marchingSquares( f )

fx = f.components{1}; 
fy = f.components{2};
pref = chebfunpref;
tol = pref.cheb2Prefs.chebfun2eps;
num = 0; 
r = zeros(1,2); 
dom = fy.domain;

nf = f.nComponents; 
if ( nf > 2 )
    error('CHEBFUN:CHEBFUN2V:roots:zeroSurface', ...
        'CHEBFUN2 is unable to find zero surfaces.');
end

if ( length(fx) == 1 || length(fy) == 1 )   % one of them is of the form u(x)v(y)
    
    if ( length(fx) == 1 )
        % Find roots of fx. 
        rowz = roots(fx.rows); 
        colz = roots(fy.cols); 
        for jj = 1:length(rowz)
            rr = roots(fy(rowz(jj),:));
            for kk = 1:length(rr)
                r(num+1,:) = [rowz(jj) rr(kk)]; num=num+1;
            end
        end
        for jj = 1:length(colz)
            rr = roots(fy(:,colz(jj)));
            for kk = 1:length(rr)
                r(num+1,:) = [rr(kk) colz(jj)]; num=num+1;
            end
        end
    elseif ( length(fy) == 1 )
        % roots lie along these lines
        rowz = roots(fy.rows); 
        colz = roots(fy.cols); 
        for jj = 1:length(rowz)
            ff = fx(rowz(jj),:);
            rr = roots(ff);
            for kk = 1:length(rr)
                r(num+1,:) = [rowz(jj) rr(kk)]; num=num+1;
            end
        end
        for jj = 1:length(colz)
            rr = roots(fx(:,colz(jj)));
            for kk = 1:length(rr)
                r(num+1,:) = [rr(kk) colz(jj)]; num=num+1;
            end
        end
    end
    
else
    % Use contourc to get roots to 3-4 digits.
    N = 400; 
    NewtonFail = 0;
    rfx = roots(fx); 
    rfy = roots(fy);
    x = linspace(-1, 1, N);  % this is on [-1,1] because all zero contours are represented by chebfuns on the default interval.
    r = zeros(0,2); 
    num = 0;
    
    Rx = rfx(x(:), :); 
    Ry = rfy(x(:), :);
    if ( any(size(Rx)==[1 1]) )
        Rx = Rx( : ); 
    end
    if ( any(size(Ry)==[1 1]) )
        Ry = Ry( : ); 
    end
    
    for jj = 1:size(rfx, 2)
        rx = Rx(:, jj); 
        rlx = real( rx ); 
        imx = imag( rx );
        for kk = 1:size(rfy,2)
            ry = Ry(:,kk);
            [x0, y0]=intersections(rlx, imx, real(ry), imag(ry));
            if ( ~isempty(x0) )
                r(num+1:num+length(x0),1)=x0;
                r(num+1:num+length(x0),2)=y0;
                num=num+length(x0);
            end
        end
    end
    % Use a few iterations of Newton to get 14-15 digits.
    f = fx; 
    g = fy;
    fx = diff(f,1,2); 
    fy = diff(f); 
    gx = diff(g,1,2); 
    gy = diff(g);            % derivatives.
    J = @(x,y) [feval(fx,x,y) feval(fy,x,y);...
                                    feval(gx,x,y) feval(gy,x,y)];  % Jacobian
    
    warnstate = warning('off','CHEBFUN:CHEBFUN2:NEWTON');   % turn warnings off, and capture Newton failure instead.
    for kk = 1:size(r,1)
        x0 = [r(kk,1), r(kk,2)].';
        dx = 1; 
        iter = 1;
        while ( norm(dx) > 10*tol && iter < 15 )
            dx = J(x0(1),x0(2)) \ -[feval(f,x0(1),x0(2));feval(g,x0(1),x0(2))];    % update
            x0 = dx + x0; iter = iter + 1;
        end
        if ( norm(dx) < 10*sqrt(tol) ) % we may have diverged so don't always update.
            r(kk,:) = x0;
        else
            NewtonFail = 1;
        end
    end
    warning(warnstate); % turn them back on.
    
    
    % If all the Newton iterations failed then some roots may be
    % inaccurate.
    if ( NewtonFail )
        warning('CHEBFUN:CHEBFUN2V:roots:newtonFail', ...
            'Iterates may have diverged some of the computed roots may be not be accurate.')
    end
end

%%
% Remove the roots which lie outside of the domain.
if ( ~isempty(r) )
    r = r( (r(:,1) <= dom(2)+tol &...
        r(:,1) >= dom(1)-tol & ...
        r(:,2) <= dom(4)+tol & ...
        r(:,2) >= dom(3)-tol ), :);
end

if num==0
    xroots=[]; yroots=[];
    return
end

xroots = r(:,1); yroots = r(:,2);


end


function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)
%INTERSECTIONS Intersections of curves.
%   Computes the (x,y) locations where two curves intersect.  The curves
%   can be broken with NaNs or have vertical segments.
%
% Example:
%   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% where X1 and Y1 are equal-length vectors of at least two points and
% represent curve 1.  Similarly, X2 and Y2 represent curve 2.
% X0 and Y0 are column vectors containing the points at which the two
% curves intersect.
%
% ROBUST (optional) set to 1 or true means to use a slight variation of the
% algorithm that might return duplicates of some intersection points, and
% then remove those duplicates.  The default is true, but since the
% algorithm is slightly slower you can set it to false if you know that
% your curves don't intersect at any segment boundaries.  Also, the robust
% version properly handles parallel and overlapping segments.
%
% The algorithm can return two additional vectors that indicate which
% segment pairs contain intersections and where they are:
%
%   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% For each element of the vector I, I(k) = (segment number of (X1,Y1)) +
% (how far along this segment the intersection is).  For example, if I(k) =
% 45.25 then the intersection lies a quarter of the way between the line
% segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for
% the vector J and the segments in (X2,Y2).
%
% You can also get intersections of a curve with itself.  Simply pass in
% only one curve, i.e.,
%
%   [X0,Y0] = intersections(X1,Y1,ROBUST);
%
% where, as before, ROBUST is optional.

% Version: 1.12, 27 January 2010
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})

% License:
%
% Copyright (c) 2008, Douglas M. Schwarz
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% Theory of operation:
%
% Given two line segments, L1 and L2,
%
%   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
%   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
%
% we can write four equations with four unknowns and then solve them.  The
% four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
% L1 and L2, t1 is the distance from the starting point of L1 to the
% intersection relative to the length of L1 and t2 is the distance from the
% starting point of L2 to the intersection relative to the length of L2.
%
% So, the four equations are
%
%    (x1(2) - x1(1))*t1 = x0 - x1(1)
%    (x2(2) - x2(1))*t2 = x0 - x2(1)
%    (y1(2) - y1(1))*t1 = y0 - y1(1)
%    (y2(2) - y2(1))*t2 = y0 - y2(1)
%
% Rearranging and writing in matrix form,
%
%  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
%        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
%   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
%        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
%
% Let's call that A*T = B.  We can solve for T with T = A\B.
%
% Once we have our solution we just have to look at t1 and t2 to determine
% whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
% line segments cross and we can include (x0,y0) in the output.
%
% In principle, we have to perform this computation on every pair of line
% segments in the input data.  This can be quite a large number of pairs so
% we will reduce it by doing a simple preliminary check to eliminate line
% segment pairs that could not possibly cross.  The check is to look at the
% smallest enclosing rectangles (with sides parallel to the axes) for each
% line segment pair and see if they overlap.  If they do then we have to
% compute t1 and t2 (via the A\B computation) to see if the line segments
% cross, but if they don't then the line segments cannot cross.  In a
% typical application, this technique will eliminate most of the potential
% line segment pairs.


% Input checks.
narginchk(2,5)

% Adjustments when fewer than five arguments are supplied.
switch nargin
    case 2
        robust = true;
        x2 = x1;
        y2 = y1;
        self_intersect = true;
    case 3
        robust = x2;
        x2 = x1;
        y2 = y1;
        self_intersect = true;
    case 4
        robust = true;
        self_intersect = false;
    case 5
        self_intersect = false;
end

% x1 and y1 must be vectors with same number of points (at least 2).
if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
        length(x1) ~= length(y1)
    error('CHEBFUN:CHEBFUN2:roots:intersections:badInputs1', ...
        'X1 and Y1 must be equal-length vectors of at least 2 points.')
end
% x2 and y2 must be vectors with same number of points (at least 2).
if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
        length(x2) ~= length(y2)
    error('CHEBFUN:CHEBFUN2:roots:intersections:badInputs2', ...
        'X2 and Y2 must be equal-length vectors of at least 2 points.')
end


% Force all inputs to be column vectors.
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);

% Compute number of line segments in each curve and some differences we'll
% need later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;
xy1 = [x1 y1];
xy2 = [x2 y2];
dxy1 = diff(xy1);
dxy2 = diff(xy2);

% Determine the combinations of i and j where the rectangle enclosing the
% i'th line segment of curve 1 overlaps with the rectangle enclosing the
% j'th line segment of curve 2.
[i,j] = find(repmat(min(x1(1:end-1),x1(2:end)),1,n2) <= ...
    repmat(max(x2(1:end-1),x2(2:end)).',n1,1) & ...
    repmat(max(x1(1:end-1),x1(2:end)),1,n2) >= ...
    repmat(min(x2(1:end-1),x2(2:end)).',n1,1) & ...
    repmat(min(y1(1:end-1),y1(2:end)),1,n2) <= ...
    repmat(max(y2(1:end-1),y2(2:end)).',n1,1) & ...
    repmat(max(y1(1:end-1),y1(2:end)),1,n2) >= ...
    repmat(min(y2(1:end-1),y2(2:end)).',n1,1));

% Force i and j to be column vectors, even when their length is zero, i.e.,
% we want them to be 0-by-1 instead of 0-by-0.
i = reshape(i,[],1);
j = reshape(j,[],1);

% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
% At the same time we can remove redundant combinations of i and j in the
% case of finding intersections of a line with itself.
if self_intersect
    remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
else
    remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
end
i(remove) = [];
j(remove) = [];

% Initialize matrices.  We'll put the T's and B's in matrices and use them
% one column at a time.  AA is a 3-D extension of A where we'll use one
% plane at a time.
n = length(i);
T = zeros(4,n);
AA = zeros(4,4,n);
AA([1 2],3,:) = -1;
AA([3 4],4,:) = -1;
AA([1 3],1,:) = dxy1(i,:).';
AA([2 4],2,:) = dxy2(j,:).';
B = -[x1(i) x2(j) y1(i) y2(j)].';

% Loop through possibilities.  Trap singularity warning and then use
% lastwarn to see if that plane of AA is near singular.  Process any such
% segment pairs to determine if they are colinear (overlap) or merely
% parallel.  That test consists of checking to see if one of the endpoints
% of the curve 2 segment lies on the curve 1 segment.  This is done by
% checking the cross product
%
%   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
%
% If this is close to zero then the segments overlap.

% If the robust option is false then we assume no two segment pairs are
% parallel and just go ahead and do the computation.  If A is ever singular
% a warning will appear.  This is faster and obviously you should use it
% only when you know you will never have overlapping or parallel segment
% pairs.

if robust
    overlap = false(n,1);
    warning_state = warning('off','MATLAB:singularMatrix');
    % Use try-catch to guarantee original warning state is restored.
    try
        lastwarn('')
        for k = 1:n
            T(:,k) = AA(:,:,k)\B(:,k);
            [unused,last_warn] = lastwarn;
            lastwarn('')
            if strcmp(last_warn,'MATLAB:singularMatrix')
                % Force in_range(k) to be false.
                T(1,k) = NaN;
                % Determine if these segments overlap or are just parallel.
                overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
            end
        end
        warning(warning_state)
    catch err
        warning(warning_state)
        rethrow(err)
    end
    % Find where t1 and t2 are between 0 and 1 and return the corresponding
    % x0 and y0 values.
    in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
    % For overlapping segment pairs the algorithm will return an
    % intersection point that is at the center of the overlapping region.
    if any(overlap)
        ia = i(overlap);
        ja = j(overlap);
        % set x0 and y0 to middle of overlapping region.
        T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
            min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
        T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
            min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
        selected = in_range | overlap;
    else
        selected = in_range;
    end
    xy0 = T(3:4,selected).';
    
    % Remove duplicate intersection points.
    [xy0,index] = unique(xy0,'rows');
    x0 = xy0(:,1);
    y0 = xy0(:,2);
    
    % Compute how far along each line segment the intersections are.
    if nargout > 2
        sel_index = find(selected);
        sel = sel_index(index);
        iout = i(sel) + T(1,sel).';
        jout = j(sel) + T(2,sel).';
    end
else % non-robust option
    for k = 1:n
        [L,U] = lu(AA(:,:,k));
        T(:,k) = U\(L\B(:,k));
    end
    
    % Find where t1 and t2 are between 0 and 1 and return the corresponding
    % x0 and y0 values.
    in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
    x0 = T(3,in_range).';
    y0 = T(4,in_range).';
    
    % Compute how far along each line segment the intersections are.
    if nargout > 2
        iout = i(in_range) + T(1,in_range).';
        jout = j(in_range) + T(2,in_range).';
    end
end

% Plot the results (useful for debugging).
% plot(x1,y1,x2,y2,x0,y0,'ok');
end


function x = mypoints(n, dom)
% Get the sample points that correspond to the right grid for a particular
% technology.

% What tech am I based on?:
tech = chebfunpref().tech();
if ( ischar(tech) )
    tech = eval(tech);
end

if ( isa(tech, 'chebtech2') )
    x = chebpts( n, dom, 2 );   % x grid.
elseif ( isa(tech, 'chebtech1') )
    x = chebpts( n, dom, 1 );   % x grid.
elseif ( isa(tech, 'trigtech') )
    x = trigpts( n, dom );   % x grid.
else
    error('CHEBFUN:CHEBFUN2V:roots:techType', 'Unrecognized technology');
end

end
