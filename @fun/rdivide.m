function f = rdivide(f, g, varargin)
%./   Right array divide for a FUN.
%   F ./ C divides the FUN F by an array C. If F is a vectorised FUN with M
%   columns, then C must be either a scalar or a 1xM array.
%
%   Alternatively if C is a FUN with the same number of columns as F, or if F is
%   a scalar, the resulting division is returned if C is found to have no roots
%   in [a, b], where [a, b] is the domain of F and C. The division is performed
%   column-wise.
%
% See also MRDIVIDE, TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% FUN ./ [] = []:
if ( isempty(f) || isempty(g) )
    f = [];
    return
end

% Look at different cases
if ( isa(g, 'double') )         % FUN ./ DOUBLE
    f.onefun = rdivide(f.onefun, g, varargin{:});
elseif isa(f, 'double')         % DOUBLE ./ FUN
    % In this case, g is the FUN. But if g has roots in [a, b], we're going to
    % get an error, so enclose in a try-catch in order to be able to throw a
    % reasonable error message
    try
        g.onefun = rdivide(f, g.onefun, varargin{:});
    catch ME
        if strcmp(ME.identifier,'CHEBFUN:CHEBTECH:rdivide:DivideByZeros')
            % The right argument has roots on [a ,b].
            error('CHEBFUN:FUN:rdivide:DivideByZeros', ...
                'Cannot divide by a FUN with roots in [%.5g, %.5g].', ...
                g.domain(1), g.domain(2));
        else
            % Unexpected error, rethrow it.
            rethrow(ME)
        end
    end
    % Assign g to be the output argument
    f = g;
else                            % FUN ./ FUN
    % Both f and g are FUN objects. Similarly as above, enclose in a try-catch
    % loop in case g has roots on [a, b].
    try
        f.onefun = rdivide(f.onefun, g.onefun, varargin{:});
    catch ME
        if strcmp(ME.identifier,'CHEBFUN:CHEBTECH:rdivide:DivideByZeros')
            % The right argument has roots on [a ,b].
            error('CHEBFUN:FUN:rdivide:DivideByZeros', ...
                'Cannot divide by a FUN with roots in [%.5g, %.5g].', ...
               g.domain(1), g.domain(2));
        else
            % Unexpected error, rethrow it.
            rethrow(ME)
        end
    end
end

end
