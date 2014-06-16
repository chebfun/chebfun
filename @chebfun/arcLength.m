function out = arcLength(f, a, b)
%ARCLENGTH	Compute the length of the arc defined by a CHEBFUN.
%   ARCLENGTH(F) returns the arc length of the curve defined by CHEBFUN F in the
%   x-y plane over the interval where it is defined. If F is a CHEBFUN of
%   complex values, the output is the arc length of the curve in the complex
%   plane.
%
%   ARCLENGTH(F, A, B) returns the arc length of F over the interval [A, B],
%   where [A, B] is a subinterval of the domain in which F is defined. In the
%   case of complex-valued F, ARCLENGTH(F, A, B) computes the length of the arc
%   whose ends correspond to A and B.
%
%   If F is a quasimatrix, the arc length of each CHEBFUN in F will be computed
%   and a vector is returned.
%
% Examples:
%   f = chebfun(@(x) sin(x), [0 1]);
%   L = arcLength(f);

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty CHEBFUN:
if ( isempty(f) )
    out = [];
    return
end

% Check the first input for its type:
if ( ~isa(f, 'chebfun') )
    error('CHEBFUN:CHEBFUN:arcLength:Input', ...
        'The first argument must be a chebfun object.')
end

% Get the interval corresponding to the arc:
if ( nargin == 3 )
    % Full arguments:
    dom(1) = a;
    dom(2) = b;
elseif ( nargin == 2 )
    % Two arguments: the second argument is a vector:
    if ( max( size(a) ) ~= 2 )
        error('CHEBFUN:CHEBFUN:arcLength:Input', ...
            'The second argument must be a 1x2 vector.')
    end
    dom = a;
else
    % Single argument:
    dom = [f(1).domain(1) f(1).domain(end)];
end

% Loop over each column or row:
fPrime = f;
for i = 1:numel(f.funs)
    % This makes sure that no delta functions are 
    % generated due to jump discontinuities. pointValus of fPrime are not
    % updated since they are not needed.
    fPrime.funs{i} = diff(f.funs{i});
end

if ( isreal(f) )
    g = sqrt(1 + fPrime.^2);
    out = sum(g, dom(1), dom(2));
    
else
    out = sum(abs(fPrime), dom(1), dom(2));
end

% Reform the output to accommodate the transposedness with the input:
if ( f(1).isTransposed )
    out = out.';
end

end