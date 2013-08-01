function fx = feval(f, x, varargin)
%FEVAL   Evaluate a chebfun at one or more points.
%
%   FEVAL(F, X) ...
%
%   FEVAL(F, X, 'left') or FEVAL(F, X, 'right')

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% [TODO]: Add support for multi-rowed impulses and deltas.

% If F or x are empty, there's nothing to do.
if ( isempty(f) || isempty(x) )
    fx = [];
    return
end

% Support for feval(f, 'left') and feval(f, 'end'), etc.
if ( ischar(x) )
    if ( any(strcmpi(x, {'left', 'start' ,'-'})) )
        fx = get(f, 'lval');
    elseif ( any(strcmpi(x, {'right', 'end', '+'})) )
        fx = get(f, 'rval');
    else
        error('CHEBFUN:feval:strInput', 'Unknown input argument "%s".', x);
    end
    return
end

% Reshape x to be a column vector.
sizex = size(x);
ndimsx = ndims(x);
x = x(:);

% Initialise:
[numFuns, numCols] = size(f);
fx = zeros(size(x, 1), numCols);
funs = f.funs;
dom = f.domain;

% Deal with feval(f, x, 'left') and feval(f, x, 'right'):
if ( nargin > 2 )
    lr = varargin{1};
    parse = strcmpi(lr, {'left', 'right', '-', '+'});
    if ( ~any(parse) )
        if ( ischar(lr) )
            error('CHEBFUN:feval:leftrightchar',...
                'Unknown input argument "%s".', lr);
        else
            error('CHEBFUN:feval:leftright', 'Unknown input argument.');
        end
    end
    % We deal with this by reassigning imps to be left/right values.
    if ( parse(1) || parse(3) ) % left
        f.impulses(:,:,2:end) = []; % Level 2 imps are not needed here
        for j = 1:numFuns
            f.impulses(j+1,1,1) = get(funs{j}, 'rval');
        end
    elseif ( parse(2) || parse(4) ) % right
        f.impulses(:,:,2:end) = []; % Level 2 imps are not needed here
        for j = 1:numFuns
            f.impulses(j,1,1) = get(funs{j}, 'lval');
        end
    end
end

% Points to the left of the domain:
xReal = real(x);
I = xReal < dom(1);
if ( any(I(:)) )
    fx(I,:) = feval(funs{1}, x(I));
end
% Points within the domain:
for k = 1:numFuns
    I = ( xReal >= dom(k) ) & ( xReal < dom(k+1) );
    if ( any(I(:)) )
        % Evaluate the appropriate fun:
        fx(I,:) = feval(funs{k}, x(I));
    end
end
% Points to the right of the domain:
I = ( xReal >= dom(end) );
if ( any(I(:)) )
    fx(I,:) =  feval(funs{end}, x(I));
end

%% DEALING WITH IMPS:
% If the evaluation point corresponds to a breakpoint, we get the value from
% imps. If there is only one row, the value is given by the corresponding entry
% in that row. If the second row is nonzero the value is -inf or inf
% corresponding to the sign of the entry in the 2nd row. If the entry in the
% corresponding 3rd or higher rows is nonzero, we return NaN.

higherImpulses = f.impulses(:,:,2:end);
% Only one row:
if ( (size(f.impulses, 3) == 1) || ~any(higherImpulses(:)) )
    % Loop over the funs:
    for k = 1:numFuns + 1
        idx = x == dom(k);
        if ( any(idx) )
            fx((x == dom(k)),:) = f.impulses(k,:,1);
        end
    end
    
else
    % [TODO]: Multiple imps rows:    
end

%% Reshape if possible:

% [TODO]: Document what we return in these cases.
if ( (ndimsx == 2) && (sizex(1) == 1) )
    fx = fx.';
elseif ( ((ndimsx > 2) || (sizex(2) > 1)) && (numCols == 1) )
    fx = reshape(fx, sizex);
elseif ( ((ndimsx == 2) || (sizex(2) > 1)) && (numCols > 1))
    fx = reshape(fx, sizex(1), numCols*numel(x)/sizex(1));
end
