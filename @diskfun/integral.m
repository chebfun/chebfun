function I = integral( f, varargin )
%INTEGRAL   Complete definite integral of DISKFUN.
%
%   I = INTEGRAL(F), returns the definite integral of a DISKFUN integrated
%   over its domain of definition.
%
%   I = INTEGRAL(F, g), returns the integral of F along the
%   curve defined by the complex-valued CHEBFUN g.
%
%   I = INTEGRAL(F, 'unitcircle') returns the integral of F along the
%   unit circle.
% See also DISKFUN/SUM2, DISKFUN/INTEGRAL2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%check for empty
if isempty(f)
    I = 0;
    return
end


% Another way to do sum2(f)
if ( nargin == 1 )                         
    I = integral2( f );
else
    if ( ischar(varargin{1}) )
        if ( strcmpi(varargin{1},'unitcircle') )
            c = chebfun(@(t) exp(1i*t), [-pi, pi]);
        else
            error('CHEBFUN:DISKFUN:INTEGRAL:unrecognizedType',...
                'Unrecognized line integral type.  Did you mean "unitcircle"');
        end
    else
        c = varargin{1};
    end
    if ( ~isa( c, 'chebfun' ) )
        I = integral2( f, varargin{ : } );
    else                            % Line integral over a CHEBFUN
        % Make complex:
        c = c + realmin*1i;
        % Line integral:
        I = sum( feval(f, c ) .* abs( diff( c ) ) );
    end    
end

end
