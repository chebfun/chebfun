function r = roots(f, g, h, varargin)
%roots   Poor man's method for finding common roots of three 3D functions:
%   1) Find some initial guesses from the tensor of values of 
%   obj_fun = f.^2 + g.^2 + h.^2.
%   2) Apply 10 steps of Newton's method.
%   3) Filter inaccurate results.
%
% See also CHEBFUN3/ROOT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

dom = f.domain;
f = f./f.vscale; g = g./g.vscale; h = h./h.vscale; 

[diffF1, diffF2, diffF3] = grad(f);
[diffG1, diffG2, diffG3] = grad(g);
[diffH1, diffH2, diffH3] = grad(h);

Jac = @(x,y,z) [feval(diffF1,x,y,z),  feval(diffF2,x,y,z),  feval(diffF3,x,y,z)
                feval(diffG1,x,y,z),  feval(diffG2,x,y,z),  feval(diffG3,x,y,z)
                feval(diffH1,x,y,z),  feval(diffH2,x,y,z),  feval(diffH3,x,y,z)];

[mF, nF, pF] = length(f);
[mG, nG, pG] = length(g);
[mH, nH, pH] = length(h);
len=max([mF, nF, pF; mG, nG, pG; mH, nH, pH]);

%% A naiive 3D subdivision
% r = [];
% if max(len) > 50
%     levels = 2;
%     domTensor = subdivide3D(dom,levels);
%     %[~,f_fiberDim] = min(size(f.core)); 
%     %[~,f_fiberDim] = min(rank(f)); 
%     [r1, r2, r3] = rank(f); [~, f_fiberDim] = min( [r1, r2, r3] ); f_fiberDim = f_fiberDim(1);
%     %[~,g_fiberDim] = min(rank(g)); 
%     [r1, r2, r3] = rank(g); [~, g_fiberDim] = min( [r1, r2, r3] ); g_fiberDim = g_fiberDim(1);
%     %[~,h_fiberDim] = min(rank(h)); 
%     [r1, r2, r3] = rank(h); [~, h_fiberDim] = min( [r1, r2, r3] ); h_fiberDim = h_fiberDim(1);
%     for i = 1:levels
%         for j = 1:levels
%             for k = 1:levels
%                 f1 = chebfun3(@(x,y,z) feval(f,x,y,z), squeeze(domTensor(i,j,k,:))','fiberDim',f_fiberDim);
%                 g1 = chebfun3(@(x,y,z) feval(g,x,y,z), squeeze(domTensor(i,j,k,:))','fiberDim',g_fiberDim);
%                 h1 = chebfun3(@(x,y,z) feval(h,x,y,z), squeeze(domTensor(i,j,k,:))','fiberDim',h_fiberDim);
%                 %r = {r, roots(f1,g1,h1,varargin)}
%                 rr = roots(f1,g1,h1,varargin);
%                 if ~isempty(rr)
%                     r = [r; rr];
%                 end
%             end
%         end
%     end
% end
%%

% bndNumRoots = (sum(length(f))-3)*(sum(length(g))-3)*(sum(length(h))-3); % Bezout's classical bound 
% bndNumRoots = (lenF(1)+lenG(1)+lenH(1)-3)*(lenF(2)+lenG(2)+lenH(2)-3)*(lenF(3)+lenG(3)+lenH(3)-3)...
% - (lenF(1)+lenG(1)-2)*(lenF(2)+lenG(2)-2)*(lenF(3)+lenG(3)-2)...
% - (lenF(1)+lenH(1)-2)*(lenF(2)+lenH(2)-2)*(lenF(3)+lenH(3)-2)...
% - (lenH(1)+lenG(1)-2)*(lenH(2)+lenG(2)-2)*(lenH(3)+lenG(3)-2)...
% + (lenF(1)-1)*(lenF(2)-1)*(lenF(3)-1)...
% + (lenG(1)-1)*(lenG(2)-1)*(lenG(3)-1)...
% + (lenH(1)-1)*(lenH(2)-1)*(lenH(3)-1);% Bernstein's bound 
% % [TODO: Alex Townsend suggests that he feels it to be better to use the bounds that are employed in continuous homotopy methods]
% %bndNumRoots = 50;

factor = 2;
% %factor = 10;
% maxNdfs = max([ndf(f) ndf(g) ndf(h)]);
% if  maxNdfs < 500
%     factor = 8;
% elseif max(maxNdfs) > 1e3
%     factor = 2;
% end
m = factor*len(1); 
n = factor*len(2); 
p = factor*len(3); 

xx = chebpts(m,dom(1:2)); 
yy = chebpts(n,dom(3:4)); 
zz = chebpts(p,dom(5:6));
[xx,yy,zz] = ndgrid(xx,yy,zz);

f2 = f./f.vscale; g2 = g./g.vscale; h2 = h./h.vscale; 

%T = chebpolyval3(obj_fun,m,n,p);
T = chebpolyval3(f2,m,n,p).^2 + chebpolyval3(g2,m,n,p).^2 + chebpolyval3(h2,m,n,p).^2;
%T = abs(chebpolyval3(f2,m,n,p)) + abs(chebpolyval3(g2,m,n,p)) + abs(chebpolyval3(h2,m,n,p));
% Common roots of f,g, and h are the same as common roots of f2, g2 and h2.


% % Find the smallest len entries in |T|.
% [TSorted TIdx] = sort(abs(T));
% %smallestNElements = TSorted(1:sum(len)-1);
% allCandInd = TIdx(1:bndNumRoots);
% %allCandInd = TIdx(1:max(len)-1);
% [allColind,allRowInd,allTubeInd] = ind2sub(size(T),allCandInd);
% allCandLoc = zeros(numel(allCandInd),3);
% for i=1:numel(allCandInd)
%     allCandLoc(i,:) = [xx(allColind(i),allRowInd(i),allTubeInd(i)),yy(allColind(i),allRowInd(i),allTubeInd(i)),zz(allColind(i),allRowInd(i),allTubeInd(i))];
% end

pref = chebfunpref; 
% Find all the small entries of |T|.
%tol1 = pref.eps^(1/6);
tol1 = pref.chebfuneps^(1/6);
%tol1 = pref.eps^(1/4);
%if min(abs(T(:)))<tol
if min(T(:))<tol1
    allSmallInd = find(T<tol1);
    %allSmallInd = find(T<1e-2);
    [allColind,allRowInd,allTubeInd] = ind2sub(size(T),allSmallInd);
    allCandLoc = zeros(size(allColind,1),3);
    
    for i=1:size(allColind,1)
        allCandLoc(i,:) = [xx(allColind(i),allRowInd(i),allTubeInd(i)),yy(allColind(i),allRowInd(i),allTubeInd(i)),zz(allColind(i),allRowInd(i),allTubeInd(i))];
    end
    
else
    warning('There might be no common roots at all!')
    r = [];
    return;
end
r = allCandLoc;

% for k = 1:size(r,1)
%     rr = r(k,:);
%     upd_norm = norm(Jac(rr(1),rr(2),rr(3))\ [feval(f,rr(1),rr(2),rr(3)); feval(g,rr(1),rr(2),rr(3)); feval(h,rr(1),rr(2),rr(3))]);
%     %if upd_norm > tol1
%         
%         
% end


NewtonFail = 0;
warning_state = warning('off', 'MATLAB:nearlySingularMatrix'); % turn warnings off, and capture Newton failure instead.
warning_state = warning('off', 'MATLAB:singularMatrix'); % turn warnings off, and capture Newton failure instead.
warning_state = warning('off', 'MATLAB:illConditionedMatrix');
%tol2 = pref.eps;
tol2 = pref.chebfuneps;
for k = 1:size(r,1)
    update = 1; iter = 1;
    r1 = [r(k,1), r(k,2), r(k,3)];
    while ( norm(update) > 10*tol2 && iter <10 )
        update = Jac(r1(1),r1(2),r1(3))\ [feval(f,r1(1),r1(2),r1(3)); feval(g,r1(1),r1(2),r1(3)); feval(h,r1(1),r1(2),r1(3))];
        r1 = r1 - update.'; 
        iter = iter + 1;
    end
    
    if ( norm(update) < 10*sqrt(tol2) ) % we may have diverged so don't always update.
        r(k,:) = r1;
    else
        NewtonFail = 1;
    end
    
end

warning(warning_state); % turn them back on.
% If all the Newton iterations failed then some roots may be
% inaccurate.
if ( NewtonFail )
    warning('CHEBFUN:CHEBFUN3:roots:newtonFail', ...
        'Iterates may have diverged some of the computed roots may be not be accurate.')
end

% Remove the roots which lie outside of the domain.
if ( ~isempty(r) )
    r = r( (r(:,1) <= dom(2)+tol2 &...
        r(:,1) >= dom(1)-tol2 & ...
        r(:,2) <= dom(4)+tol2 & ...
        r(:,2) >= dom(3)-tol2 & ... 
        r(:,3) <= dom(6)+tol2 & ...
        r(:,3) >= dom(5)-tol2),: );
    % Remove the roots which are inacuurate.
    r = r( abs(feval(f,r(:,1), r(:,2),r(:,3))) <tol2^(2/3) & ...
           abs(feval(g,r(:,1), r(:,2),r(:,3))) <tol2^(2/3) & ...
           abs(feval(h,r(:,1), r(:,2),r(:,3))) <tol2^(2/3),:);
end

indices = [];
for k=1:size(r,1)-1
    rr = r(k,:);
    for kk=k+1:size(r,1)        
        if ( norm(rr-r(kk,:)) < tol2^(2/3) )
            indices = [indices; kk];
        end
    end
end
indices = unique(indices);

r(indices,:) = [];
if size(r,1)>1
    if norm(r(1,:)-r(end,:)) < tol2^(2/3)
        r(1,:) = [];
    end
end

end


function domTensor = subdivide3D(dom,levels)
domTensor = zeros(levels,levels,levels,6); % Each subvolume will be a 6*1 column vector.

xDiv = linspace(dom(1),dom(2),levels+1);
yDiv = linspace(dom(3),dom(4),levels+1);
zDiv = linspace(dom(5),dom(6),levels+1);

for j = 1:length(xDiv)-1
    for k = 1:length(yDiv)-1
        for l = 1:length(zDiv)-1
            domTensor(j,k,l,:) = [xDiv(j:j+1) yDiv(k:k+1) zDiv(l:l+1)];
        end
    end
end

end


function deg = topDeg(f,g,h,dom)
temp1 = [feval(f,dom(1),dom(3),dom(5)), feval(f,dom(1),dom(4),dom(5)), feval(f,dom(1),dom(3),dom(6)),feval(f,dom(1),dom(4),dom(6))];

f2 = chebfun2();
f2.domain = dom(3:6);
[U,S,V] = svd(squeeze(chebfun3.txm(f.core,f.cols(dom(1)),1)));
f2.rows = f.rows*U;
f2.cols = f.tubes*V;
f2.pivotValues = 1./diag(S);
[val,~] = minandmax2(f2);
temp1 = [temp1, val(1), val(2)];
fLs1 = min(temp1);
fUs1 = max(temp1);


temp2 = [feval(f,dom(2),dom(3),dom(5)), feval(f,dom(2),dom(4),dom(5)), feval(f,dom(2),dom(3),dom(6)),feval(f,dom(2),dom(4),dom(6))];
fLs2 = min(temp2);
fUs2 = max(temp2);

temp3 = [feval(f,dom(1),dom(3),dom(5)), feval(f,dom(1),dom(3),dom(6)), feval(f,dom(2),dom(3),dom(5)),feval(f,dom(2),dom(3),dom(6))];
fLs3 = min(temp3);
fUs3 = max(temp3);

temp4 = [feval(f,dom(1),dom(4),dom(5)), feval(f,dom(1),dom(4),dom(6)), feval(f,dom(2),dom(4),dom(5)),feval(f,dom(2),dom(4),dom(6))];
fLs4 = min(temp4);
fUs4 = max(temp4);

temp5 = [feval(f,dom(1),dom(3),dom(5)), feval(f,dom(1),dom(4),dom(5)), feval(f,dom(2),dom(3),dom(5)),feval(f,dom(2),dom(4),dom(5))];
fLs5 = min(temp5);
fUs5 = max(temp5);

temp6 = [feval(f,dom(1),dom(3),dom(6)), feval(f,dom(1),dom(4),dom(6)), feval(f,dom(2),dom(3),dom(6)),feval(f,dom(2),dom(4),dom(6))];
fLs6 = min(temp6);
fUs6 = max(temp6);

end