function f = rdivide(f, c, pref)
%./   Right array divide for a FUNCHEB1
%   F ./ C Divides the FUNCHEB1 F by an array C. If F is a vectorised
%   FUNCHEB1 with M columns, then C must be either a scalar or a 1xM array.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(c, 'double') )
    % Dividing by a constant is easy.
    
    % This can never work (as size(f, 1) == inf):
    if ( size(c, 1) > 1 || ...
           ( numel(c) > 1 && size(f.values, 2) ~= size(c, 2) ) )
        error('CHEBUFN:FUNCHEB1:rdivide:size', ...
            'Matrix dimensions must agree.');
    end
    
    if ( ~any(c(:)) )  
        f = funcheb1(NaN(1, size(f, 2)));
        
    elseif ( numel(c) == 1 )
        % Scalar
        f.values = f.values/c;      % Divide values
        f.coeffs = f.coeffs/c;      % Divide coeffs
        f.vscale = f.vscale/abs(c); % Divide vscale
        
    else
        % Vectorised 
        n = size(f.values, 1);   
        f.values = f.values./repmat(c, n, 1);   % Divide values
        f.coeffs = f.coeffs./repmat(c, n, 1);   % Divide coeffs
        f.vscale = f.vscale/norm(c, inf);       % Divide vscale
        
        f.values(:, c == 0) = NaN;
        f.coeffs(:, c == 0) = NaN;
        
    end
    
else
    % Dividing by another FUNCHEB1 is harder. Call compose.
    
    % Obtain preferencess:
    if ( nargin < 3 )
        pref = funcheb1.pref;
    end
    
    % Check if c (the denominator) has any roots.
    if ( ~isempty(roots(c)) )
        error('CHEBFUN:FUNCHEB1:rdivide:DivideByZeros', ...
        'Division by a FUNCHEB1 with roots in [-1, 1].');
    else
%         warning('CHEBFUN:FUNCHEB1:RDIVIDE:DivideByFunCHEB1', ...
%         'Division by a FUNCHEB1');
    end
    
    % Call COMPOSE.
    if ( isa(f, 'funcheb1') )  % funcheb1 / funcheb1
        f = compose(f, @rdivide, c, pref);
    else                   % double / funcheb1
        op = @(x) f./x;
        f = compose(c, op, pref);
    end
end

end
