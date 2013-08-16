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

% Three cases: real, imaginary, and complex.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REAL CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( isreal(f) )

    %% %%%%%%%%% LOCATE NEW BREAKS %%%%%%%%%
    % Abs is singular at roots, so locate these:
    r = roots(f, 'nozerofun', 'nojump', 'noimps');
    
    % Since each column of an array-valued CHEBFUN must have the same breakpoints,
    % we simply take unique(r(:)) and remove any remaining NaNs.
    r = unique(r(:));
    r(isnan(r)) = [];
    
    % Discard any roots which are closer than the accuracy of the CHEBFUN:
    el = epslevel(f);
    vs = vscale(f);
    hs = hscale(f);
    rtol1 = el*hs*vs;
    r([false ; diff(r) < rtol1]) = [];
    
    % Avoid introducing new breakpoints close to an existing ones:
    rtol2 = el*vs*max(min(diff(f.domain)), 1);
    r(any(abs(bsxfun(@minus, r, f.domain)) < rtol2, 2)) = [];

    %% %%%%%%%%%% CREATE NEW FUNS %%%%%%%%%%
    % Get the domain with the new breakpoints: (union is not required, by above)
    oldDom = f.domain;
    newDom = sort([f.domain, r.']);

    % Introduce these breakpoints to f:
    f = restrict(f, newDom);

    % Call ABS on each of the FUNs: (result will be smooth)
    for k = 1:numel(f.funs)
        f.funs{k} = abs(f.funs{k});
    end
    
    %% %%%%%%%%%% TIDY THE OUTPUT %%%%%%%%%%
    
    % Take the absolute value of the impulses in the first row:
    f.impulses(:,:,1) = abs(f.impulses(:,:,1));
    
    % Simplify the output:
    f = simplify(f);
    
    % [TODO]: Do we want to do this?
    [ignored, idx] = setdiff(newDom, oldDom);
    f = merge(f, idx.'); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGINARY CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( isreal(1i*f) )

    % 1i*f will be real and abs(1i*f) = abs(f), so call ABS again with 1i*f.
    f = abs(1i*f, pref);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPLEX CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else

    % SQRT will deal with introducing breakpoints.
    imps = f.impulses(:,:,1);
    f = sqrt(real(conj(f).*f), pref);
    
    % Take the absolute value of the impulses in the first row:
    f.impulses = abs(imp);

end




end