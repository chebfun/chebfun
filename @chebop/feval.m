function out = feval(N, varargin)
%FEVAL   Evaluate the operator of the CHEBOP at a CHEBFUN or CHEBMATRIX.
%   OUT = FEVAL(N, U) for a CHEBFUN or CHEBMATRIX U applies the CHEBOP N to U,
%   i.e., it returns N(U). Here, N.OP should be of the form @(u) diff(u,2) + ...
%   If N.op is of the form @(x, u) diff(u,2) + ... then an x variable is
%   instantiated internally and included automatically, however, this is not
%   the prefered syntax and may not be supported in future releases.
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
%   the prefered syntax.
%
%   OUT = FEVAL(N, DIM) returns an DIM-point discretization of the linear
%   operator N. If N is not linear an error is thrown. OUT = FEVAL(N, DIM,
%   'oldschool') forces the returned differentiation matrices to be square,
%   rather than rectangular. See CHEBOP/MATRIX for further details (which is the
%   prefered syntax for this functionality)
%
% See also CHEBOP/SUBSREF, LINOP/MTIMES, CHEBOP/MATRIX.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Support for calling a linear CHEBOP with a numerical input to get its
% discretization. Here we simply call the MATRIX() method.

if ( ((nargin == 2) && isnumeric(varargin{1})) || ...
   ((nargin == 3) && ischar(varargin{2}) && strcmpi(varargin{2}, 'oldschool')) )
    % TODO: Throw a deprecated warning?
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
        % Construct the independent variable X.
        x = chebfun(@(x) x, N.domain);
        u = varargin{1};
        
    elseif ( numel(varargin) == numberOfInputs )
        % Got passed both x and u.
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
    
    if ( numel(args) == numberOfInputs - 1 )
        % Check if we need to include an x (independent variable):
        x = chebfun(@(x) x, N.domain);
        args = [ {x} , args ];
    end
    
    % Evaluate the operator:
    out = N.op(args{:});
end

end
