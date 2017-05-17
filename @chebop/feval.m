function out = feval(N, varargin)
%FEVAL   Evaluate the operator of the CHEBOP at a CHEBFUN or CHEBMATRIX.
%   OUT = FEVAL(N, U) for a CHEBFUN or CHEBMATRIX U applies the CHEBOP N to U,
%   i.e., it returns N(U). Here, N.OP should be of the form @(u) diff(u,2) + ...
%   If N.op is of the form @(x, u) diff(u,2) + ... then an x variable is
%   instantiated internally and included automatically, however, this is not
%   the preferred syntax and may not be supported in future releases.
%
%   OUT = FEVAL(N, X, U) for the CHEBFUN X and CHEBFUN or CHEBMATRIX U applies
%   the CHEBOP N to X and  U, i.e., it returns N(X, U) where N.OP has the form
%   @(x, u) diff(u,2) + .... Here, X should be the dependent variable on
%   N.DOMAIN.
%
%   OUT = FEVAL(N, X, U1, U2, ..., UM) for a CHEBFUN X and CHEBFUN or CHEMBATRIX
%   objects U1, ..., UM applies the CHEBOP N to the functions Uk; i.e., it
%   returns N(X, U1, U2, ..., UM) where N.OP has the form @(x, u1, u2, ..., um).
%   Note that for systems of equations, X _must_ be included in N.OP.
%
%   OUT = FEVAL(N, X, U) where U is a CHEBMATRIX of M entries and N.OP has the
%   form @(X, U1, U2, ..., UM) is equivalent FEVAL(N, X, U{1}, ..., U{M}).
%   Again, OUT = FEVAL(N, U) will also work in this situation, but this is not
%   the preferred syntax.
%
%   OUT = FEVAL(N, DIM) returns an DIM-point discretization of the linear
%   operator N. If N is not linear an error is thrown. OUT = FEVAL(N, DIM,
%   'oldschool') uses boundary bordering to deal wiuth the boundary conditions,
%   rather than the rectangular projection approach of hale and Driscoll. Note
%   that this syntax exists only the support ATAP and doesn't necessarily give a
%   clear picture of the discretizations now being used in the Chebfun release.
%   CHEBOP/MATRIX, which is the preferred syntax for this functionality,
%   provides further details
%
% See also CHEBOP/SUBSREF, LINOP/MTIMES, CHEBOP/MATRIX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Support for calling a linear CHEBOP with a numerical input to get its
% discretization. This is deprecated
if ( (nargin == 2) && isnumeric(varargin{1}) )
    warning('CHEBFUN:CHEBOP:feval:deprecated', ...
        ['FEVAL(N, DIM) or N(DIM) exists only to provide backwards\n', ...
        'compatibility with ATAP. The preferred method for visualizing a\n', ...
        'discretization of a linear CHEBOP is MATRIX(N, DIM). Note, however,\n', ...
        'that these may not give the same result due to changes in how\n', ...
        'CHEBOP discretizes differential operators.'])
    warning('off', 'CHEBFUN:CHEBOP:feval:deprecated');

    % We cannot support boundary conditions in this way.
    throwBCwarning = false;
    if ( ~isempty(N.lbc) )
        N.lbc = [];
        throwBCwarning = true;
    end
    if ( ~isempty(N.rbc) )
        N.rbc = [];
        throwBCwarning = true;
    end
    if ( ~isempty(N.bc) )
        N.bc = [];
        throwBCwarning = true;
    end
    if ( throwBCwarning )
        warning('CHEBFUN:CHEBOP:feval:BCs', ...
            'Boundary conditions are not supported in FEVAL(N, DIM).')
    end

    % Because the native V5 matrix is rectangular, we have to do some resizing
    % to recapture V4 behavior. Note: this may break for systems, piecewise
    % cases.
    n = varargin{1};
    L = linop(N);
    pref = cheboppref();
    disc = pref.discretization;
    if ( isa(disc, 'char') && strcmpi(disc, 'values') )
        pref.discretization = @chebcolloc2;
    elseif ( isa(disc, 'function_handle') && ~isequal(disc, @chebcolloc2) ) 
        % We only support CHEBCOLLOC2 discretizations!
        error('CHEBFUN:CHEBOP:feval:notColloc2', ...
            ['FEVAL(N, DIM) only supports CHEBCOLLOC2 discretizations.\n', ...
             'Use MATRIX(N, DIM) or change the discretization in CHEBOPPREF.']);
    end
    A = matrix(L, n, pref); % has n rows but maybe more columns

    % We want an n by n result. We have that A is (n-d) x n where d =
    % L.diffOrder. Also, the output of A is at 1st kind points. We need to do:
    %  (map from n 1st kind to n 2nd kind) * A * (map from n 2nd kind to
    %  n+d 2nd kind).

    % Note that we don't have to translate/scale to the domain for barymat
    % matrices that go between grids.
    d = L.diffOrder;
    [x1,~,w1] = chebtech1.chebpts( n );
    x2 = chebtech2.chebpts( n );
    x3 = chebtech2.chebpts( n + d );

    out = barymat(x2, x1, w1) * A * barymat(x3, x2);
    return

elseif ((nargin == 3) && ischar(varargin{2}) ...
        && strcmpi(varargin{2}, 'oldschool'))

    % TODO: Is this still correct?
    out = matrix(N, varargin{:});
    return
end

% We must expand CHEBMATRIX entries out to a cell for {:} to work below.
isChebMat = cellfun(@(u) isa(u, 'chebmatrix'), varargin);
if ( any(isChebMat) )
    args = {};
    for k = 1:numel(varargin)
        % Append variables from the kth input:
        if ( isChebMat(k) )
            args = [args , varargin{k}.blocks.']; %#ok<AGROW>
        else
            args = [args , varargin(k)];          %#ok<AGROW>
        end
    end
else
    args = varargin;
end

% How many input arguments are there to N.op?
numberOfInputs = nargin(N);

% If no arguments, return empty:
if ( numberOfInputs == 0 )
    out = [];
    return
end

if ( numberOfInputs == 1)
    % If N has one input arguments, either we have a scalar problem, or the
    % problem is specified on chebmatrix format, e.g.,
    %   N.op = @(u) [ diff(u{1},2) + u{2}; u{1} + diff(u{2}];
    % Here, importantly, x does not appear in the argument list for N.op.
    u = varargin{1};

    % If we have a scalar problem, but U is still a CHEBMATRIX, we need to
    % extract the BLOCK of U in order to be evaluate N. If on the other hand, U
    % has more than one block, but NUMBEROFINPUTS is still equal to 1 (which got
    % us here in the first place), we must be dealing with a CHEBOP N, whose OP
    % is specified on CHEBMATRIX format, e.g.
    %   N.op = @(u) [diff(u{1}) + u{2}; u{1} + diff(u{2})];
    if ( isa(u, 'chebmatrix') && max(size(u)) == 1 )
        u = u.blocks{1};
    end

    out = N.op(u);

elseif ( numberOfInputs == 2 )
    % If N has two input arguments, either we have a scalar problem, or the
    % problem is specified on chebmatrix format, e.g.,
    %   N.op = @(x, u) [ diff(u{1},2) + u{2}; u{1} + diff(u{2}];
    % Here, importantly, x must appear in the argument list for N.op.

    % Did we not get the x variable passed in as argument?
    if ( numel(varargin) == 1 )
        u = varargin{1};

        % Construct the independent variable X.
        x = chebfun(@(x) x, N.domain);
        % If the CHEBFUN U passed in is a quasimatrix, we need to tile the
        % independent variable X to have the same dimensions as U:
        x = repmat(x, 1, size(u,2));

    elseif ( numel(varargin) == numberOfInputs )
        % Got passed both X and U.
        x = varargin{1};
        u = varargin{2};

    else
        error('CHEBFUN:CHEBOP:feval:numInputs', ...
            'Unexpected number of input arguments.')
    end

    % If we have a scalar problem, but U is still a CHEBMATRIX, we need to
    % extract the BLOCK of U in order to be evaluate N. If on the other
    % hand, U has more than one block, but NUMBEROFINPUTS is still less than
    % or equal to 2 (which got us here in the first place), we must be
    % dealing with a CHEBOP N, whose OP is specified on CHEBMATRIX format,
    % e.g.,
    %   N.op = @(x,u) [diff(u{1}) + u{2}; u{1} + diff(u{2})];
    if ( isa(u, 'chebmatrix') && max(size(u)) == 1 )
        u = u.blocks{1};
    end

    % Evaluate the operator!
    out = N.op(x, u);
    
else
    % The operator is specified on the form
    %   N.op = @(x, u, v) = [diff(u,2) + v; u + diff(v)]
    
    % Count the number of RHSs:
    numCols = max(max(cellfun(@(v) size(v, 2), varargin)));

    % We must expand CHEBMATRIX entries out to a cell for {:} to work below.
    isChebMat = cellfun(@(u) isa(u, 'chebmatrix'), varargin);
    if ( any(isChebMat) )
        args = {};
        for k = 1:numel(varargin)
            % Append variables from the kth input:
            if ( isChebMat(k) )
                args = [args , varargin{k}.blocks.']; %#ok<AGROW>
            else
                args = [args , varargin(k)];          %#ok<AGROW>
            end
        end
        % ARGS need to be vertically concatenated for operator to be evaluated
        % correctly.
        args = args.';
    else
        % ARGS need to be vertically concatenated for operator to be evaluated
        % correctly.
        args = varargin.';
    end
    
    if ( numCols > 1 )
        % Deal with multiple RHS
        out = cell(1, numCols);
        for k = 1:numCols
            out{k} = doEval(N, args(:,k));
        end
        out = horzcat(out{:});
    else
        out = doEval(N, args);
    end
    
end

% If the operator is written in a function file or nested function, then it
% is natural for it to be returned in cell form. Convert it to a chebmatrix:
if ( iscell(out) )
    out = vertcat(out{:}); % TODO: is this correct for multiple RHSs?
end

end

function out = doEval(N, args)

    numberOfInputs = nargin(N);
    if ( numel(args) == numberOfInputs - 1 )
        % Check if we need to include an x (independent variable):
        x = chebfun(@(x) x, N.domain);
        args = [ {x} ; args ];
    end
    % Evaluate the operator:
    out = N.op(args{:});
    
end
