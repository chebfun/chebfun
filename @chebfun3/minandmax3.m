function [vals, locs] = minandmax3(f)
%MINANDMAX3   returns the global minimum and maximum value of a CHEBFUN3.
%   VALS = minandmax3(F) returns the global minimum and maximum value of a 
%   CHEBFUN3 object F over its domain. VALS is a vector of length 2 such 
%   that Y(1) = min3(F(x,y,z)) and Y(2) = max3(F(x,y,z)).
%
%   [VALS, LOCS] = minandmax3(F) also returns the position of the global 
%   minimum and maximum.
%
% See also CHEBFUN3/MAX, CHEBFUN3/MAX2, CHEBFUN3/NORM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% check for empty CHEBFUN3.
if ( isempty(f) )
    vals = []; 
    locs = [];
    return
end

doNewton = 1;   % Newton polishing?

% Extract low rank representation of f:
[fCore, fCols, fRows, fTubes] = tucker(f);
dom = f.domain;

% Is f the zero function?
if ( iszero(f) ) 
    locs = [(dom(1) + dom(2))/2  (dom(3) + dom(4))/2  (dom(5) + dom(6))/2];
    locs = [locs; locs];
    vals = [0; 0];
end

[m, n, p] = length(f);
if ( ndf(f) > 5e4 )
    m = min(m, 129);
    n = min(n, 129);
    p = min(p, 129);
end

% We seek a fast initial guess. So we first discretize the object.
xpts = chebpts(m, fCols.domain);
ypts = chebpts(n, fRows.domain);
zpts = chebpts(p, fTubes.domain);
colVals = feval(fCols, xpts); 
rowVals = feval(fRows, ypts); 
tubeVals = feval(fTubes, zpts); 
T = chebfun3.txm(chebfun3.txm(chebfun3.txm(fCore,colVals,1), rowVals,2), ...
    tubeVals,3);
    
% Minimum entry in the sample:
[ignored, ind] = min(T(:));
[col, row, tub] = ind2sub(size(T), ind);
loc(1,1) = xpts(col);
loc(1,2) = ypts(row);
loc(1,3) = zpts(tub);
vals(1) = feval(f, loc(1, 1), loc(1, 2), loc(1, 3));

% Maximum entry in the sample:
[ignored, ind] = max(T(:)); 
[col, row, tub] = ind2sub(size(T), ind);
loc(2,1) = xpts(col);
loc(2,2) = ypts(row);
loc(2,3) = zpts(tub);
vals(2) = feval(f, loc(2, 1), loc(2, 2), loc(2, 3));

% Get more digits with optimisation algorithms.
lb = [dom(1); dom(3); dom(5)];
ub = [dom(2); dom(4); dom(6)];

if ( ~isempty(ver('optim')) )
% Matlab's Optimization Toolbox is available. So, use fmincon command.
    options = optimset('Display', 'none', 'TolFun', eps, 'TolX', eps, ...
        'algorithm', 'active-set');
    [minLoc, vals(1)] = fmincon(@(x) feval(f, x(1), x(2), x(3)), ...
        loc(1, :), [], [], [], [], lb, ub, [], options);
    
    [maxLoc, vals(2)] = fmincon(@(x) -feval(f, x(1), x(2), x(3)), ...
        loc(2,:), [], [], [], [], lb, ub, [], options);
    vals(2) = -vals(2);
    loc(1,:) = minLoc;
    loc(2,:) = maxLoc;
    locs = loc;
else
    try
        % Use core Matlab fminsearch command by converting to an 
        % unconstrained optimization problem.
        % Maps from [-1, 1] to dom(1:2), dom(3:4), and dom(5:6) 
        % respectively.
        map1 = bndfun.createMap(dom(1:2));
        map2 = bndfun.createMap(dom(3:4));
        map3 = bndfun.createMap(dom(5:6));
        
        % Unconstrained initial guesses:
        loc(:, 1) = asin(map1.Inv(loc(:, 1)));
        loc(:, 2) = asin(map2.Inv(loc(:, 2)));
        loc(:, 3) = asin(map3.Inv(loc(:, 3)));
        
        % Maps from R to dom(1:2), dom(3:4), and dom(5:6) respectively.
        map1 = @(x) map1.For(sin(x));
        map2 = @(x) map2.For(sin(x));
        map3 = @(x) map3.For(sin(x));
        
        % Set options:
        options = optimset('Display', 'off', 'TolFun', eps, 'TolX', eps);
        warnstate = warning;
        warning('off'); % Disable verbose warnings from fminsearch.
        f_mapped = @(x) feval(f, map1(x(1)), map2(x(2)),map3(x(3)));
        
        [minLoc, vals(1)] = fminsearch(@(x) f_mapped(x), loc(1, :), options);
        [maxLoc, vals(2)] = fminsearch(@(x) -f_mapped(x), loc(2, :), options);
        
        vals(2) = -vals(2);
        loc(1:2, 1) = map1([minLoc(1); maxLoc(1)]);
        loc(1:2, 2) = map2([minLoc(2); maxLoc(2)]);
        loc(1:2, 3) = map3([minLoc(3); maxLoc(3)]);
        warning(warnstate);
        locs = loc;
        
    catch
        % Nothing is going to work. So we will have to go with initial 
        % guesses.
    end
end
    
if ( doNewton )
    % Store values before applying any Newton iterations to restore if 
    % Newton was diverging.
    locsOld = locs;
    valsOld = vals;
    
    % If the global max or min is already on the edge or out of domain, 
    % do NOT apply Newton iteration:
    if ( isOutOfDomain(locs, dom) )
        return
    end    
    
    % A few steps of Newton for optimization (involves computing Hessian)
    H = Hessian(f);
    gradF = grad(f);
    gradF = @(x,y,z) gradF(x,y,z);
    
    % Disable verbose warnings from Newton step if e.g., the Hessian 
    % matrix is singular which happens e.g., if a 2D function is passed to
    % minandmax3.    
    warning_state = warning('off','MATLAB:singularMatrix');
	% Use try-catch to guarantee original warning state is restored.
    try
        lastwarn('')
        k = 1; % k Newton iterations:
        for iter = 1:k
            locs(1,:) = newton_opt(locs(1,:), H, gradF);
            vals(:,1) = feval(f, locs(1, 1), locs(1, 2), locs(1, 3));
            
            locs(2,:) = newton_opt(locs(2, :), H, gradF);
            vals(:,2) = feval(f,locs(2, 1), locs(2, 2), locs(2, 3));
        end
        [ignored, last_warn] = lastwarn;
        if strcmp(last_warn,'MATLAB:singularMatrix')
            % Return the old value if the Hessian matrix was singular.
            locs = locsOld;
            vals = valsOld;
            return
        end
        warning(warning_state)
    catch err
        warning(warning_state)
        rethrow(err)
    end
   
    % If Newton's method moved us outside of the domain, then return old 
    % value.
    if ( isOutOfDomain(locs, dom) )
        locs = locsOld;
        vals = valsOld;
    end
    
end
   
end

%%
function H = Hessian(f)
% Forms the Hessian matrix of a CHEBFUN3 object f

H11 = diff(f, 2, 1);
H12 = diff(diff(f, 1, 2), 1, 1); 
H13 = diff(diff(f, 1, 3), 1, 1);

H21 = diff(diff(f, 1, 1), 1, 2);
H22 = diff(f, 2, 2);
H23 = diff(diff(f, 1, 3), 1, 2);

H31 = diff(diff(f, 1, 1), 1, 3); 
H32 = diff(diff(f, 1, 2), 1, 3);
H33 = diff(f, 2, 3);
H = @(x,y,z) [feval(H11, x, y, z)  feval(H12, x, y, z)  feval(H13, x, y, z);
              feval(H21, x, y, z)  feval(H22, x, y, z)  feval(H23, x, y, z);
              feval(H31, x, y, z)  feval(H32, x, y, z)  feval(H33, x, y, z)];
end

%%
function sol = newton_opt(init, H, gradF)
% Perform one step of multivariate optimization using Newton's method.
x = init(1); 
y = init(2); 
z = init(3);
rhs = -gradF(x, y, z);
h = H(x,y,z)\rhs;
sol = init.' + h;
end

%%
function check = isOutOfDomain(locs, dom)
% Check whether both points in the rows of LOCS are inside the cuboid DOM.
if ( locs(1, 1) <= dom(1) || locs(1, 1) >=dom(2) || ...
         locs(1, 2) <= dom(3) || locs(1, 2)>=dom(4) || ...
         locs(1, 3) <= dom(5) || locs(1, 3)>=dom(6) )
     check = true;
elseif ( locs(2, 1) <= dom(1) || locs(2, 1) >= dom(2) || ...
             locs(2, 2) <= dom(3) || locs(2, 2) >= dom(4) || ...
             locs(2, 3) <= dom(5) || locs(2, 3) >= dom(6) )
         check = true;
else
    check = false;
end
end
%%