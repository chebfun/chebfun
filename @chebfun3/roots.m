function r = roots(f, g, h)
%ROOTS   Finding common roots of three CHEBFUN3 object.
%   Here are the 3 steps in this algorithm:
%   1) Find some initial guesses for the roots. This is done using a tensor
%      of values of the objective function obj_fun = f.^2 + g.^2 + h.^2. 
%      Specifiaclly, our initial guesses are local minima (valleys) in 
%      that discrete tensor found e.g. by "findpeaks".
%   2) Apply some steps of Newton's method to improve initial guesses.
%   3) Filter out inaccurate results.
%
% See also CHEBFUN3/ROOT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Scale the functions, as common roots of (f, g, h) are the same as 
% (f2, g2, h2).
f = f./f.vscale;
g = g./g.vscale; 
h = h./h.vscale; 

[diffF1, diffF2, diffF3] = grad(f);
[diffG1, diffG2, diffG3] = grad(g);
[diffH1, diffH2, diffH3] = grad(h);

Jac = @(x,y,z) [...
    feval(diffF1, x, y, z),  feval(diffF2, x, y, z),  feval(diffF3, x, y, z)
    feval(diffG1, x, y, z),  feval(diffG2, x, y, z),  feval(diffG3, x, y, z)
    feval(diffH1, x, y, z),  feval(diffH2, x, y, z),  feval(diffH3, x, y, z)];

[mF, nF, pF] = length(f);
[mG, nG, pG] = length(g);
[mH, nH, pH] = length(h);
len = max([mF, nF, pF; mG, nG, pG; mH, nH, pH]);

%% % Bezout's classical bound 
% bndNumRoots = (sum(length(f))-3)*(sum(length(g))-3)*(sum(length(h))-3); 
% bndNumRoots = (lenF(1)+lenG(1)+lenH(1)-3)*(lenF(2)+lenG(2)+lenH(2)-3)*(lenF(3)+lenG(3)+lenH(3)-3)...
% - (lenF(1)+lenG(1)-2)*(lenF(2)+lenG(2)-2)*(lenF(3)+lenG(3)-2)...
% - (lenF(1)+lenH(1)-2)*(lenF(2)+lenH(2)-2)*(lenF(3)+lenH(3)-2)...
% - (lenH(1)+lenG(1)-2)*(lenH(2)+lenG(2)-2)*(lenH(3)+lenG(3)-2)...
% + (lenF(1)-1)*(lenF(2)-1)*(lenF(3)-1)...
% + (lenG(1)-1)*(lenG(2)-1)*(lenG(3)-1)...
% + (lenH(1)-1)*(lenH(2)-1)*(lenH(3)-1);% Bernstein's bound 
% % [TODO: Alex Townsend suggests that he feels it to be better to use the bounds that are employed in continuous homotopy methods]
% %bndNumRoots = 50;
%% Sample f^2 + g^2 + h^2 on an apropraite grid:
factor = 2;
m = factor*len(1);
n = factor*len(2);
p = factor*len(3);

% If the function is easy enough, increase size of the sampling grid as 
% much as you tolerable. See issue #1900.
% if ( max([ndf(f), ndf(g), ndf(h)]) < 5e4 )
%     m = max(m, 51);
%     n = max(n, 51);
%     p = max(p, 51);
% end

dom = f.domain;
xx = chebpts(m, dom(1:2));
yy = chebpts(n, dom(3:4));
zz = chebpts(p, dom(5:6));
[xx,yy,zz] = ndgrid(xx,yy,zz);

%T = chebpolyval3(obj_fun,m,n,p);
T = chebpolyval3(f,m,n,p).^2 + chebpolyval3(g,m,n,p).^2 + chebpolyval3(h,m,n,p).^2;
%T = abs(chebpolyval3(f,m,n,p)) + abs(chebpolyval3(g,m,n,p)) + abs(chebpolyval3(h,m,n,p));

%% Find initial guesses: 
%[Maxima,MaxIdx] = findpeaks(T(:));
%[Maxima,MaxIdx] = findpeaks(T(:),'MinPeakDistance',6);
%DataInv = 1.01*max(T(:)) - T(:); % Keep all the entries of the inverted vector positive.
%[ignored,allSmallInd] = findpeaks(DataInv);
[allSmallInd,ignored]=findpeaks2(T(:),'v'); % Find valleys in T(:)
[v,x,t,m,ze]=quadpeak(T);

%[Minima,allSmallInd] = findpeaks(DataInv, 'MinPeakDistance',6);
[allColind,allRowInd,allTubeInd] = ind2sub(size(T), allSmallInd);
allCandLoc = zeros(size(allColind,1),3);
for i=1:size(allColind,1)
    allCandLoc(i,:) = [xx(allColind(i), allRowInd(i), allTubeInd(i)), ...
                       yy(allColind(i), allRowInd(i), allTubeInd(i)), ...
                       zz(allColind(i), allRowInd(i), allTubeInd(i))];
end
r = allCandLoc;

NewtonFail = 0;
warning_state = warning('off', 'MATLAB:nearlySingularMatrix'); % turn warnings off, and capture Newton failure instead.
warning_state = warning('off', 'MATLAB:singularMatrix'); % turn warnings off, and capture Newton failure instead.
warning_state = warning('off', 'MATLAB:illConditionedMatrix');
tol2 = chebfunpref().chebfuneps;
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
    %TODO: A root with 1 correct digit is still better to keep than to remove...
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