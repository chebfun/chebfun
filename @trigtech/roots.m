function out = roots(f, varargin)
%ROOTS   Roots of a TRIGTECH in the interval [-1,1].
%   ROOTS(F) returns the real roots of the TRIGTECH F in the interval [-1,1].
%
%   ROOTS(F, PROP1, VAL1, PROP2, VAL2, ...) modifies the default ROOTS
%   properties. The PROPs (strings) and VALs may be any of the following:
%
%   ALL: 
%       [0] - Return only real-valued roots in [-1,1].
%        1  - Return roots outside of [-1,1] (including complex roots).
%
%   RECURSE:
%        0  - Compute roots without interval subdivision (slower).
%       [1] - Subdivide until length(F) < 50. (causes additional complex roots).
%
%   If F is an array-valued TRIGTECH then there is no reason to expect each
%   column to have the same number of roots. In order to return a useful output,
%   the roots of each column are computed and then padded with NaNs so that a
%   matrix may be returned. The columns of R = ROOTS(F) correspond to the
%   columns of F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    out = [];
    return
end

% Here is the present strategy:
% 
% 1. The default option is to simply create a chebtech of the TRIGTECH and
%    call roots on it. 
%
% 2. If the option is to find 'all' or 'complex' roots then we simply use
%    matlab's built-in roots routine.
%
% [TODO]: Figure out a more elegant way to find roots of a TRIGTECH.

useMatlabsRootsCommand = false;
pruneRoots = false;

if ( nargin > 1 )
    if isa(varargin{1}, 'struct')
        rootsPref = varargin{1};
        useMatlabsRootsCommand = rootsPref.all;
        pruneRoots = rootsPref.prune;
    else
        j = 1;
        while ( j <= length(varargin) )
            % Determine if the 'all' or 'complex' flag was passed in.
            if ( any(strcmp(lower(varargin{j}), 'complex')) ) %#ok<STCI>
                useMatlabsRootsCommand = varargin{j+1};
                pruneRoots = true;
                break;
            elseif ( any(strcmp(lower(varargin{j}), 'all')) ) %#ok<STCI>
                useMatlabsRootsCommand = varargin{j+1};                
                break;
            end
            j = j+2;
        end
    end
end

if ( useMatlabsRootsCommand )
    numCols = size(f.coeffs, 2);
    r = cell(1,numCols);
    for j = 1:numCols
        % Simplify the current column to get the minimal number of
        % roots.
        fj = simplify(extractColumns(f,j));
        % Flip coeffs to match Matlab roots:
        rTemp = roots(fj.coeffs(end:-1:1,:));
        % Roots here finds the roots in the transformed variable z=exp(i*pi*x)
        % so we need to take the log (and scale it) to get back the roots in x.
        rTemp = -1i/pi*log(rTemp);
       
        % Prune the roots if required.  
        if ( pruneRoots )
            % For trigtech pruning means only keep the roots that fit in
            % the strip of analyticity (estimated using 10*machine eps).
            % This is the analog to what chebtech does by pruning to the
            % Bernstein ellipse.  The distance from the real-axis of the
            % strip is given by a = 1i/N/pi*log(4/e + 1), where e is the
            % defined max precision (10*eps)) and N is the highest non-zero
            % Fourier term.
            N = ceil(length(fj.coeffs)/2) - 1;
            a = 1/N/pi*log(4/(10*eps) + 1);
            rTemp = rTemp(abs(imag(rTemp)) <= a); 
        end
        
        r{j} = rTemp;
    end
    % Find the max length of r:
    mlr = max(cellfun(@length, r)); 

    % Pad the columns in r with NaNs:
    r = cellfun(@(x) [x ; NaN(mlr - length(x), 1)], r, 'UniformOutput', false);

    % Convert to an array for output:
    out = cell2mat(r);
else
    % An arbitrary decision was made to use CHEBTECH1.
    g = chebtech1(@(x) f.feval(x));
    out = roots(g,varargin{:});
end

end
