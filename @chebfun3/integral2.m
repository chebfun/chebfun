function I = integral2(f, C)
%INTEGRAL2  Surface integral: Integral of the CHEBFUN3 object F over the 
%           surface defined by a CHEBFUN2 object C. By default C represents
%           z = f(x,y). [TODO: Implement x = f(y,z) and y = f(x,z)].
%
%   I = INTEGRAL2(F) returns a value representing the double integral of a
%   CHEBFUN3.
%
%   I = INTEGRAL2(F, [a b c d e g]) integrate F over the cube [a b] x [c
%   d] x [e g] provided that this cube is in the domain of F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) ) 
    I = 0;
    return
end

%if ( nargin == 1 )    
%     % Double definite integral:
%     I = sum( cols ) * D * sum( rows ).';
%    elseif ( nargin == 2 )
if ( nargin == 2 )        
    if ( isa( C, 'chebfun2' ) ) 
        % Integral over a surface represented by a CHEBFUN2.
        % Get surface: 
        %C = varargin{1};

        % Surface integral: 
        diffCx = diff(C, 1, 2);
        diffCy = diff(C, 1, 1);
        %ds_squared = 1+ diffCx.^2 + diffCy.^2;
        op = @(x,y) feval(f, x, y, feval(C, x, y)) .* ...
            sqrt(1 + feval(diffCx, x, y).^2 + feval(diffCy, x, y).^2);
        %I = sum2( chebfun2(op, C.domain));
        ICheb = chebfun2(op, C.domain);
        I = integral(ICheb);
    end
else
    error('CHEBFUN:CHEBFUN3:integral2:nargin', 'Incorrect number of input arguments.');
end

end