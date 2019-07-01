function out = matrix(L, varargin)
%   MATRIX(L, DIM) returns an DIM by (DIM+D) discretization of the linear
%   chebop L, where D is the differential order of L. If L is not linear
%   then an error is thrown. 
%
%   If L has assigned boundary conditions then these are included in the
%   discretization, increasing the number of rows. If the domain of L has
%   breakpoints then DIM must be a vector specifying the discretization size
%   on each interval. In this case, appropriate continuity constraints are
%   included in the discretization, further increasing the number of rows
%   (by D times the number of breakpoints).
%
%   MATRIX(L, DIM, PREFS) allows additional preferences to be passed 
%   via the CHEBOPPREF PREFS.
%
%   MATRIX(L, DIM, 'oldschool') forces the returned differentiation
%   matrices to be square rather than rectangular. See LINOP/MATRIX and
%   LINOP/FEVAL for further details.
%
%   MATRIX(L, DIM, 'oldschool', PREFS) allows, again, additional 
%   preferences to be passed via the CHEBOPPREF PREFS.
% 
%   For mathematical details of these matrices see Aurentz and Trefethen,
%   "Block operators and spectral discretizations," SIAM Review, 2017.
%   For the 'oldschool' discretizations, see Trefethen, Spectral Methods
%   in MATLAB, SIAM, 2000.
%
% See also FEVAL, LINOP/MATRIX, LINOP/FEVAL.

% Get the preferences if given:
isPrefGiven = 0;
for j = 1:nargin-1
    item = varargin{j};
    if ( isa(item,'cheboppref') )
        prefs = item;
        isPrefGiven = 1;
    end
end

% Otherwise, use default prefs
if ( ~isPrefGiven )
    prefs = cheboppref();
end

% Linearize:
[L_linop, ~, fail] = linop(L);
if ( fail )
    error('CHEBFUN:CHEBOP:matrix:nonlinear',...
        'Matrix expansion is only allowed for linear CHEBOP objects.')
end

% Determine the discretization:
prefs = determineDiscretization(L, L_linop, prefs);

% Derive continuity conditions in the case of breakpoints:
L_linop = deriveContinuity(L_linop);

% Call LINOP/FEVAL or LINOP/MATRIX:
if ( (numel(varargin) > 1) && strcmpi(varargin{2}, 'oldschool') )
    warnState = warning('off', 'CHEBFUN:LINOP:feval:deprecated');
    out = feval(L_linop, varargin{:});
    warning(warnState);
else
    % Add the preferences in vargarin to pass them to LINOP/MATRIX:
    if ( isPrefGiven )
        % If a CHEBOPPREF was passed to the method, it will have been at the
        % last position of varargin, indexed at nargin-1. Overwrite it with the
        % current PREFS.DISCRETIZATION, as the discretization might have changed
        % in the periodic case:
        varargin{nargin-1} = prefs;
    else
        % Otherwise, add the PREFS.DISCRETIZATION to VARARGIN, so that it can be
        % passed to the call to LINOP/MATRIX below.
        varargin{nargin} = prefs;
    end
    out = matrix(L_linop, varargin{:});
end

end
