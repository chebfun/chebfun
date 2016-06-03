function g = std2(f, varargin)
%STD2   Standard deviation of a CHEBFUN3 along two variables.
%   G = STD2(F) returns the standard deviation of F in the x and y 
%   variables (default). If F is defined on the cuboid 
%   [a, b] x [c, d] x [e, g] then
%
%                            b d
%                           / /
%     G.^2 = 1/((b-a)(d-c)) | | (F(x,y,z) - mean2(F, [1, 2]))^2 dx dy
%                           / /
%                           a c
%
%   The output G is a CHEBFUN object over the interval [e, g].
%
%   G = STD2(F, FLAG, DIMS) takes the standard deviation along the
%   variables x and y, if DIMS = [1, 2] or DIMS = [2, 1], along x and z if 
%   DIMS = [1, 3] or DIMS = [3, 1] and along y and z if DIMS = [2, 3] or 
%   DIMS = [3, 2]. The FLAG is ignored and kept in this function so the 
%   syntax agrees with the Matlab STD command.
%
% See also CHEBFUN/STD, CHEBFUN2/STD, CHEBFUN2/STD2, CHEBFUN3/STD and 
% CHEBFUB3/STD3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) )
    g = chebfun();
    return
end

dom = f.domain; 
if ( nargin < 3 )
    dim1 = 1;   % default to std over x and y.
    dim2 = 2;
elseif ( nargin == 3 )
    dims = varargin{2}; 
    dim1 = dims(1);
    dim2 = dims(2);
else
    error( 'CHEBFUN:CHEBFUN3:std2:nargin', 'Too many input arguments.' ); 
end

%   Recall that std(X) = sqrt( mean( (X - mean(X) )^2 ) ).
if ( (dim1 == 1 && dim2 == 2) || (dim1 == 2 && dim2 == 1) )
    % std2 over x and y: 
    mxy = chebfun3(@(x,y,z) feval(mean2(f, [1 2]), z), dom);
    g = sqrt(1/(diff(dom(1:2))*diff(dom(3:4))) * sum2((f - mxy).^2, [1 2]));
    
elseif ( (dim1 == 1 && dim2 == 3) || (dim1 == 3 && dim2 == 1) )
    % mean2 over x and z: 
    mxz = chebfun3(@(x,y,z) feval(mean2(f, [1 3]), y), dom);
    g = sqrt(1/(diff(dom(1:2))*diff(dom(5:6))) * sum2((f - mxz).^2, [1 3]));
    
elseif ( (dim1 == 2 && dim2 == 3) || (dim1 == 3 && dim2 == 2) )
    % mean2 over y and z: 
    myz = chebfun3(@(x,y,z) feval(mean2(f, [2 3]), x), dom);
    g = sqrt(1/(diff(dom(3:4))*diff(dom(5:6))) * sum2((f - myz).^2, [2 3]));
else
    error('CHEBFUN:CHEBFUN3:std2:dims', ...
        'Third argument should contain two distinct values from 1, 2 or 3.');
end

end