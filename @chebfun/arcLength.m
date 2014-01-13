function out = arcLength(f, a, b)
%ARCLENGTH	compute the length of the arc defined by a CHEBFUN.
%   ARCLENGTH(F) returns the arc length of the curve defined by CHEBFUN F in the
%   x-y plane over the interval where it is defined. If F is a CHEBFUN of 
%   complex values, the output is the arc length of the curve in the complex 
%   plane.
%
%   ARCLENGTH(F,A,B) returns the arc length of F over the interval [A B], where 
%   [A, B] is a subinterval of the domain in which F is defined. In the case of
%   complex-valued F, ARCLENGTH(F,A,B) computes the length of the arc whose ends
%   correspond to A and B.   
%
%   If F is a quasimatrix, the arc length of each CHEBFUN in F will be
%   computed and a vector is returned.
%
%   Examples:
%   f = chebfun(@(x) sin(x), [0 1]);
%   L = arclength(f);

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% TODO: Handle array-valued chebfuns/quasimatrices.

% Empty CHEBFUN:
if ( isempty(f) )
    out = [];
    return
end

% Check the first input for its type:
if ( ~isa(f,'chebfun') )
    error('CHEBFUN:arclength:Input', ...
        'The first argument must be a chebfun object.')
end

if ( nargin == 3 )
    
    % Full arguments:
    dom(1) = a;
    dom(2) = b;
    
elseif ( nargin == 2 )
    
    % Two arguments: the second argument is a vector:
    if ( max( size(a) ) ~= 2 )
        error('CHEBFUN:arclength:Input', ...
            'The second argument must be a 1x2 vector.')
    end
    dom = a;
    
else
    % Single argument:
    dom = [f.domain(1) f.domain(end)];
end

% numCols = min(size(f));
% out = zeros(1, numCols);
    
% Loop over each column or row:
fPrime = diff(f);
if ( isreal(f) )
    out = sum(sqrt(1 + fPrime.^2), dom(1), dom(2));
else
    out = sum(abs(fPrime), dom(1), dom(2));
end

% Reform the output to accommodate the transposedness with the input:
if ( f.isTransposed )
    out = out.';
end

end