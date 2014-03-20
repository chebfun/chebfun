function f = rdivide(f, c, pref)
%./   Right array divide for a FOURIERTECH.
%   F ./ C divides a FOURIERTECH F by an array C. If F is an array-valued FOURIERTECH
%   with M columns, then C must be either a scalar or a 1xM array.
%
%   Alternatively C can be a FOURIERTECH and F can either be a FOURIERTECH with the
%   same number of columns as C or a scalar.  In this case, C must have no
%   roots in [-pi, pi], or else F ./ C may return garbage with no warning.  The
%   division is performed column-wise.
%
% See also MRDIVIDE, TIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(c, 'double') )
    % Dividing by a constant is easy.
    
    % This can never work (as size(f, 1) == inf):
    if ( (size(c, 1) > 1) || ...
           ( (numel(c) > 1) && (size(f.values, 2) ~= size(c, 2)) ) )
        error('CHEBFUN:FOURIERTECH:rdivide:size', ...
            'Matrix dimensions must agree.');
    end
    
    if ( ~any(c(:)) )  
        % Division by zero produces a NaN FOURIERTECH:
        f = f.make(NaN(1, size(f, 2)));
    elseif ( numel(c) == 1 )
        % Scalar
        f.values = f.values/c;      % Divide values
        f.coeffs = f.coeffs/c;      % Divide coeffs
        f.vscale = f.vscale/abs(c); % Divide vscale
    else
        % Array-valued FOURIERTECH
        n = size(f.values, 1);   
        f.values = f.values./repmat(c, n, 1);   % Divide values
        f.coeffs = f.coeffs./repmat(c, n, 1);   % Divide coeffs
        f.vscale = f.vscale./abs(c);            % Divide vscale
        
        f.values(:, c == 0) = NaN;
        f.coeffs(:, c == 0) = NaN;
        f.vscale(:, c == 0) = NaN;
    end
else
    % Dividing by another fouriertech is harder. Call COMPOSE.
    
    % Obtain preferences:
    if ( nargin < 3 )
        pref = fouriertech.techPref(); % c is a fouriertech.
    end

    % Call COMPOSE.
    if ( isa(f, 'fouriertech') )   % FOURIERTECH / FOURIERTECH
        f = compose(f, @rdivide, c, pref);
    else                       % DOUBLE / FOURIERTECH
        op = @(x) f./x;
        f = compose(c, op, [], pref);
    end
end

end
