function out = matrix(N, varargin)
%   OUT = MATRIX(N, DIM) returns an DIM-point discretization of the linear
%   operator N. If N is not linear an error is thrown. 
%
%   OUT = MATRIX(N, DIM, PREFS) allows additional preferences to be passed 
%   via the CHEBOPPREF, PREFS.
%
%   OUT = MATRIX(N, DIM, 'oldschool') forces the returned differentiation
%   matrices to be square, rather than rectangular. See LINOP/MATRIX and
%   LINOP/FEVAL for further details.
%
%   OUT = MATRIX(N, DIM, 'oldschool', PREFS) allows, again, additional 
%   preferences to be passed via the CHEBOPPREF, PREFS.
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
[L, ~, fail] = linop(N);
if ( fail )
    error('CHEBFUN:CHEBOP:matrix:nonlinear',...
        'Matrix expansion is only allowed for linear CHEBOP objects.')
end

% Determine the discretization:
prefs = determineDiscretization(N, L, prefs);


% Call LINOP/FEVAL or LINOP/MATRIX:
if ( (numel(varargin) > 1) && strcmpi(varargin{2}, 'oldschool') )
    warnState = warning('off', 'CHEBFUN:LINOP:feval:deprecated');
    out = feval(L, varargin{:});
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
    out = matrix(L, varargin{:});
end

end
