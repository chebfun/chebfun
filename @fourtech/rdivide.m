function f = rdivide(f, c, pref)
%./   Right array divide for a FOURTECH.
%   F ./ C divides a FOURTECH F by an array C. If F is an array-valued FOURTECH
%   with M columns, then C must be either a scalar or a 1xM array.
%
%   Alternatively C can be a FOURTECH and F can either be a FOURTECH with the
%   same number of columns as C or a scalar.  In this case, C must have no
%   roots in [-1, 1], or else F ./ C may return garbage with no warning.  The
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
        error('CHEBFUN:FOURTECH:rdivide:size', ...
            'Matrix dimensions must agree.');
    end
    
    if ( ~any(c(:)) )  
        % Division by zero produces a NaN FOURTECH:
        f = f.make(NaN(1, size(f, 2)));
    elseif ( numel(c) == 1 )
        % Scalar.
        f.values = f.values/c;      % Divide values.
        f.coeffs = f.coeffs/c;      % Divide coeffs.
        f.vscale = f.vscale/abs(c); % Divide vscale.       
    else
        % Array-valued FOURTECH.
        n = size(f.values, 1);   
        f.values = f.values./repmat(c, n, 1);   % Divide values.
        f.coeffs = f.coeffs./repmat(c, n, 1);   % Divide coeffs.
        f.vscale = f.vscale./abs(c);            % Divide vscale.
        
        f.values(:, c == 0) = NaN;
        f.coeffs(:, c == 0) = NaN;
        f.vscale(:, c == 0) = NaN;
    end
    f.isReal = f.isReal & isreal(c);
    
else
    % Dividing by another fourtech is harder. Call COMPOSE.
    
    % Obtain preferences:
    if ( nargin < 3 )
        pref = fourtech.techPref(); % c is a FOURTECH.
    end
    
    % Call COMPOSE.
    if ( isa(f, 'fourtech') )   % Possibly FOURTECH / FOURTECH.
        if ( isa(c, 'fourtech') )
            f = compose(f, @rdivide, c, pref);
        else
            error('CHEBFUN:FOURTECH:rdivide:fourtechRdivideUnknown',...
                'rdivide does not know how to divide a FOURTECH and a %s.', class(c));
        end
    else                       % DOUBLE / FOURTECH.
        op = @(x) f./x;
        f = compose(c, op, [], pref);
    end
    
end

end