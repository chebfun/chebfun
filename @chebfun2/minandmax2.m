function [Y, X] = minandmax2( f )
%MINANDMAX2, find global minimum and maximum of a chebfun2.
%
% Y = minandmax2(F) returns the minimum and maximum value of a chebfun2 over
% its domain. Y is a vector of length 2 such that Y(1) = min(f(x,y)) and
% Y(2) = max(f(x,y)).
%
% [Y X] = minandmax2(F) also returns the position of the minimum and maximum.
% That is,
%
%  F(X(1,1),X(1,2)) = Y(1)     and      F(X(2,1),X(2,2)) = Y(2)
%
% See also MAX2, MIN2, NORM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty(f) ) % check for empty chebfun2.
    Y = []; X = []; 
    return
end

% Use bivariate rootfinding to find all the local extrema: 
F = gradient( f ); 
r = roots( F ); 
if ( ~isempty( r ) ) 
    [inMax, idInMax] = max( f( r(:,1), r(:,2) ) );
    [inMin, idInMin] = min( f( r(:,2), r(:,2) ) ); 
else
    inMax = inf;   % max and min must occur on boundary. 
    inMin = inf; 
end

% Search along boundary: 
dom = f.domain;
left = f(dom(1), :); 
right = f(dom(2), :);
down = f(:, dom(3)); 
up = f(:, dom(4)); 
[Yleft, Xleft] = minandmax2( left ); 
[Yright, Xright] = minandmax2( right ); 
[Yup, Xup] = minandmax2( up ); 
[Ydown, Xdown] = minandmax2( down );

% Store Min/Max location for later: 
BcMinLocations = [ Xleft(1,:) ; Xright(1,:) ; Xup(1,:) ; Xdown(1,:) ]; 
BcMaxLocations = [ Xleft(2,:) ; Xright(2,:) ; Xup(2,:) ; Xdown(2,:) ]; 

[ BcMax, idBcMax ] = max( [ Yleft, Yright, Yup, Ydown ] );
[ BcMin, idBcMin ] = min( [ Yleft, Yright, Yup, Ydown ] );

% What is the global min and max?: 
[Ymax, inOrOutMax] = max( inMax, BcMax ); 
[Ymin, inOrOutMin] = min( inMin, BcMin );
Y = [Ymin, Ymax]; 
X = zeros(2);

% Unravel to find locations:
if ( inOrOutMin == 1 )
    X(1,:) = [ r(idInMin,1), r(idInMin,2) ];
else
    X(1,:) = BcMinLocations(idBcMin, :);
end
if ( inOrOutMax == 1 )
    X(2,:) = [ r(idInMax,1), r(idInMax,2) ];
else
    X(2,:) = BcMaxLocations(idBcMax, :);
end

end