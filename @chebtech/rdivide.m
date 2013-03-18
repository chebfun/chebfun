function f = rdivide(f, c, pref)
%./   Right array divide for a CHEBTECH.
%   F ./ C Divides the CHEBTECH F by a an array C. If F is a vectorised CHEBTECH
%   with M columns, then C must be either a scalar or a 1xM array. 
%
%   Alternatively if C is a CHEBTECH with the same number of columns as F, or if
%   F is a scalar, the resulting division is returned if C is found to have no
%   roots in [-1,1]. The division is performed columnwise.
%
% See also MRDIVIDE, TIMES.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(c, 'double') )
    % Dividing by a constant is easy.
    
    % This can never work (as size(f, 1) == inf):
    if ( size(c, 1) > 1 || ...
           ( numel(c) > 1 && size(f.values, 2) ~= size(c, 2) ) )
        error('CHEBUFN:CHEBTECH:rdivide:size', ...
            'Matrix dimensions must agree.');
    end
    
    if ( ~any(c(:)) )  
        % Division by zero produces a NaN chebtech:
        f = f.make(NaN(1, size(f, 2)));
        
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
        f.vscale = f.vscale./abs(c);            % Divide vscale
        
        f.values(:, c == 0) = NaN;
        f.coeffs(:, c == 0) = NaN;
        f.vscale(:, c == 0) = NaN;
        
    end
    
else
    % Dividing by another CHEBTECH is harder. Call COMPOSE.
    
    % Obtain preferencess:
    if ( nargin < 3 )
        pref = chebtech.pref; % c is a CHEBTECH.
    end
    
    % Check if c (the denominator) has any roots.
    if ( ~isempty(roots(c)) )
        error('CHEBFUN:CHEBTECH:rdivide:DivideByZeros', ...
        'Cannot divide by a CHEBTECH with roots in [-1, 1].');
    else
%         warning('CHEBFUN:CHEBTECH:RDIVIDE:DivideByFun2', ...
%         'Division by a CHEBTECH');
    end
    
    % Call COMPOSE.
    if ( isa(f, 'chebtech') )   % chebtech / chebtech
        f = compose(f, @rdivide, c, pref);
    else                       % double / chebtech
        op = @(x) f./x;
        f = compose(c, op, pref);
    end
end

end
