function f = abs(f, pref)
%ABS    Absolute value of a chebfun.
%   ABS(F) is the absolute value of the chebfun F.
%
% See also SIGN, ANGLE, UNWRAP, HYPOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Trivial case: (f is empty)
if ( isempty(f) )
    return
end

% Grab some preferences if none are given:
if ( nargin < 2 )
    pref = chebfun.pref();
end

% Three cases: Real, imaginary, and complex.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REAL CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( isreal(f) )

    % Abs is singular at roots, so locate these:
    r = roots(f, 'nozerofun');
    
    % Since each column of an array-valued CHEBFUN must have the same breakpoints,
    % we simply take unique(r(:)) and remove any remaining NaNs.
    r = unique(r(:));
    r(isnan(r)) = [];

    % Also discard any roots that are closer than the accuracy of the CHEBFUN:
    el = epslevel(f);
    vs = vscale(f);
    hs = hscale(f);
    r([false ; diff(r) < el*hs*vs]) = [];

    % Avoid adding new breaks where not needed:
    if ( ~isempty(r) )
        
        % Choose a tolerance:
        rtol = el*vs*max(min(diff(f.domain)), 1);

        % Remove if sufficiently close to an existing break points:
        [rem, ignored] = find(abs(bsxfun(@minus, r , f.domain)) < rtol);
        r(rem) = [];

    end

    % Get the domain with the new breakpoints: (union is not required, by above)
    dom = sort([f.domain, r']);

    % Introduce these breakpoints to f:
    f = restrict(f, dom);
    f = simplify(f);

    % Call ABS on each of the funs: (result will be smooth)
    for k = 1:numel(f.funs)
        f.funs{k} = abs(f.funs{k});
    end
    
    % Impulses are dealt with below.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGINARY CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( isreal(1i*f) )

    % 1i*f will be real and abs(1i*f) = abs(f), so call ABS again with 1i*f.
    f = abs(1i*f, pref);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPLEX CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else

    % SQRT will deal with introducing breakpoints.
    f = sqrt(conj(f).*f, pref);

end

% Take the absolute value of the impulses in the first row:
f.impulses(:,:,1) = abs(f.impulses(:,:,1));

end