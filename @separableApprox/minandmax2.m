function [Y, X] = minandmax2( f )
%MINANDMAX2   Find global minimum and maximum of a SEPARABLEAPPROX.
%   Y = minandmax2(F) returns the minimum and maximum value of a SEPARABLEAPPROX over
%   its domain. Y is a vector of length 2 such that Y(1) = min(f(x,y)) and Y(2)
%   = max(f(x,y)).
%
%   [Y, X] = minandmax2(F) also returns the position of the minimum and 
%   maximum. For example,
%
%       F(X(1,1),X(1,2)) = Y(1)     and      F(X(2,1),X(2,2)) = Y(2)
%
% See also MAX2, MIN2, NORM.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% check for empty SEPARABLEAPPROX.
if ( isempty( f ) )
    Y = []; 
    X = [];
    return
end

% This is a default implementation for separableApprox children that do not
% have know how to properly implement this method.

% Maximum possible sample matrix size:
maxsize = 4e3; 

% Is the function the zero function?
if ( iszero( f )  ) 
    dom = f.domain;
    X = [ (dom(2) + dom(1))/2 (dom(4) + dom(3))/2 ];
    X = [ X ; X ];
    Y = [0 ; 0];
end

% Extract low rank representation:
frows = f.rows;
fcols = f.cols;
piv = f.pivotValues;
dom = f.domain;

% Convert rows and columns to chebfuns.
frows = chebfun(frows,dom(1:2));
fcols = chebfun(fcols,dom(3:4));

% Share out scaling:
sgn = sign( piv ).';
sq = 1 ./ sqrt( abs( piv ) );
frows = frows * diag( sq.'.*sgn );
fcols = fcols * diag( sq );


if ( length(f) == 1 ) % Rank-1 is easy:
    % We can find it from taking maximum and minimum in x and y direction.
    
    % Find minandmax of rows and columns:
    [yr, xr] = minandmax( frows );
    [yc, xc] = minandmax( fcols );
    % All possible combinations:
    vv = [yr(1)*yc(1), yr(1)*yc(2), yr(2)*yc(1), yr(2)*yc(2)];
    
    [Y(2), indmx] = max( vv );
    [Y(1), indmn] = min( vv );
    
    % Work out the location of the maximum.
    X = zeros(2);
    X(1,1) = xr(2);
    X(1,2) = xc(2);
    
    if ( indmn <= 2 )
        X(1,1) = xr(1);
    end
    if ( mod(indmn,2) == 1 )
        X(1,2) = xc(1);
    end
    X(2,1) = xr(2);
    X(2,2) = xc(2);
    
    if ( indmx <= 2 )
        X(2,1) = xr(1);
    end
    
    if ( mod(indmx,2) == 1 ),
        X(2,2) = xc(1);
    end
    
elseif ( length(f) <= maxsize )
    
    % We seek a fast initial guess. So we first truncate the SEPARABLEAPPROX.
    ypts = chebpts(length(fcols), fcols.domain); 
    xpts = chebpts(length(frows), frows.domain);
    cvals = feval(fcols, ypts); 
    rvals = feval(frows, xpts); 

    A = cvals*rvals.';
    % Maximum entry in discretisation.
    [ignored, ind] = min( A(:) ); %#ok<ASGLU>
    [row, col] = ind2sub(size(A), ind);
    X(1,1) = xpts( col );
    X(1,2) = ypts( row );
    Y(1) = feval( f, X(1,1), X(1,2) );
    % Minimum entry in discretisation.
    [ignored, ind] = max(A(:)); %#ok<ASGLU>
    [row, col] = ind2sub(size(A), ind);
    X(2,1) = xpts(col);
    X(2,2) = ypts(row);
    Y(2) = feval(f, X(2,1), X(2,2));
    
    % Get more digits with optimisation algorithms.
    lb = [ dom(1) ; dom(3) ];
    ub = [ dom(2) ; dom(4) ];
    
    try
        
        % If the optimization toolbox is available then use it to get a better maximum.

        warnstate = warning;
        warning('off'); %#ok<WNOFF> % Disable verbose warnings from fmincon.
        options = optimset('Display', 'none', 'TolFun', eps, 'TolX', eps, ...
            'algorithm', 'active-set');
        [mn, Y(1)] = fmincon(@(x,y) feval(f, x(1), x(2)), X(1, :), ...
            [], [], [], [], lb, ub, [], options);
        [mx, Y(2)] = fmincon(@(x) -feval(f, x(1), x(2)), X(2,:), ...
            [], [], [], [], lb, ub, [], options);
        Y(2) = -Y(2);
        X(1,:) = mn;
        X(2,:) = mx;
        warning(warnstate);
        
    catch
        
        try
            % Try converting to an unconstrained problem and using built-in solver.
            
            % Maps from [-1, 1] to [dom(1:2)] and [dom(3:4)], respectively.
            map1 = bndfun.createMap(dom(1:2));
            map2 = bndfun.createMap(dom(3:4));
            % Unconstrained initial guesses:
            Z(:,1) = asin(map1.Inv(X(:,1)));
            Z(:,2) = asin(map2.Inv(X(:,2)));
            % Maps from R to [dom(1), dom(2)] and [dom(3), dom(4)], respectively.
            map1 = @(x) map1.For(sin(x));
            map2 = @(x) map2.For(sin(x));
            % Set options:
            options = optimset('Display', 'off', 'TolFun', eps, 'TolX', eps);
            warnstate = warning;
            warning('off'); %#ok<WNOFF> % Disable verbose warnings from fminsearch.
            f_mapped = @(x) feval(f, map1(x(1)), map2(x(2)));
            [mn, Y(1)] = fminsearch(@(x) f_mapped(x), Z(1, :), options);
            [mx, Y(2)] = fminsearch(@(x) -f_mapped(x), Z(2, :), options);
            Y(2) = -Y(2);            
            X(1:2,1) = map1([mn(1) ; mx(1)]);
            X(1:2,2) = map2([mn(2) ; mx(2)]);
            warning(warnstate);
            
        catch
            
            % Nothing is going to work so initial guesses will have to do.
            
        end
        
    end
    
    
elseif ( length(f) >= maxsize )
    
    error('CHEBFUN:SEPARABLEAPPROX:minandmax2:length', 'Rank is too large.');
    
end


end

%%% 
% Use the approach below when bivariate rootfinding is fully implemented. 
%%%
% % Use bivariate rootfinding to find all the local extrema:
% F = gradient( f );
% r = roots( F );
% if ( ~isempty( r ) )
%     [inMax, idInMax] = max( feval(f, r(:,1), r(:,2) ) );
%     [inMin, idInMin] = min( feval(f, r(:,2), r(:,2) ) );
% else
%     inMax = inf;   % max and min must occur on boundary.
%     inMin = inf;
% end
%
% % Search along boundary:
% dom = f.domain;
% left = feval(f, dom(1), ':');
% right = feval(f, dom(2), ':');
% down = feval(f, ':', dom(3));
% up = feval(f, ':', dom(4));
% [Yleft, Xleft] = minandmax( left );
% [Yright, Xright] = minandmax( right );
% [Yup, Xup] = minandmax( up );
% [Ydown, Xdown] = minandmax( down );
%
% % Store Min/Max location for later:
% BcMinLocations = [ Xleft(1,:) ; Xright(1,:) ; Xup(1,:) ; Xdown(1,:) ];
% BcMaxLocations = [ Xleft(2,:) ; Xright(2,:) ; Xup(2,:) ; Xdown(2,:) ];
%
% [ BcMax, idBcMax ] = max( [ Yleft(2), Yright(2), Yup(2), Ydown(2) ].' );
% [ BcMin, idBcMin ] = min( [ Yleft(1), Yright(1), Yup(1), Ydown(1) ].' );
%
% % What is the global min and max?:
% [Ymax, inOrOutMax] = max( [inMax, BcMax].' );
% [Ymin, inOrOutMin] = min( [inMin, BcMin].' );
% Y = [Ymin, Ymax];
% X = zeros(2);
%
% % Unravel to find locations:
% if ( inOrOutMin == 1 )
%     X(1,:) = [ r(idInMin,1), r(idInMin,2) ];
% else
%     X(1,:) = BcMinLocations(idBcMin, :);
% end
% if ( inOrOutMax == 1 )
%     X(2,:) = [ r(idInMax,1), r(idInMax,2) ];
% else
%     X(2,:) = BcMaxLocations(idBcMax, :);
% end
%
% end
