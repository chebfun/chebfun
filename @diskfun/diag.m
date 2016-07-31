function F = diag(f,varargin)
%DIAG(F)   Diagonal of a DISKFUN.
%   G = DIAG(F) returns the CHEBFUN representing g(x) = f(pi/4, r).
%
%   G = diag(F,T) returns the CHEBFUN representing g(x) = f(T, r), where 
%  -pi < T < pi. 
%
% See also DISKFUN/TRACE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%[varargout{1:nargout}] = diag@separableApprox(varargin{:});

if ( isempty( f ) ) 
    F = chebfun;
    return
end 

f.coords = 'polar';

%when c = 0, choose the diagonal radial slice for t = pi/4
if varargin 
    F = f(pi/4,:); 
else 
c = varargin{1};
    if abs(c) > pi
        error('CHEBFUN:DISKFUN:diag: diagonal parameter must be between -pi and pi.')
    else    
    F = f(varargin{1},:); 
    end
end

 
end
