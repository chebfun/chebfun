function r = range(f, varargin)
%RANGE   Range of CHEBFUN.
%   R = RANGE(F) returns the range R = max(F) - min(F) of the CHEBFUN F.
%
%   R = RANGE(F, DIM) operates along the dimension DIM of the quasimatrix F. If
%   DIM represents the continuous variable, then R is a vector. If DIM
%   represents the discrete dimension, then R is a quasimatrix. The default for
%   DIM is 1, unless F has a singleton dimension, in which case DIM is the
%   continuous variable.
%
% See also CHEBFUN/MINANDMAX.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www..chebfun.org/ for Chebfun information.

% TODO: This method needs a test.

% Compute the MIN and MAX:
r = minandmax(f, varargin{:});

if ( isnumeric(r) ) 
    % Range along continuous dimension.
    
    if ( numel(f) == 1 )
        % Scalar f:
        r = diff(r.');  
    elseif ( ~f(1).isTransposed )
        % Column quasimatrix:
        r = diff(r);    
    else   
        % Row quasimatrix:
        r = diff(r, 1, 2); 
    end
    
else
    % Range across discrete dimension.
    
    if ( ~r(1).isTransposed )
        % Column CHEBFUN:
        r = r(:,2) - r(:,1);
    
    else
        % Row CHEBFUN:
        r = r(2,:) - r(1,:);
        
    end
    
end

end
