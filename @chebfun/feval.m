function fx = feval(f, x, varargin)
%FEVAL   Evaluate a chebfun at one or more points.
%
%   FEVAL(F, X) ...
%
%   FEVAL(F, X, 'left') or FEVAL(F, X, 'right')

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% [TODO]: Add support for multi-valued functions.
% [TODO]: Add support for multi-rowed impulses and deltas.

% If F or x are empty, there's nothing to do. (fx defined above)
if ( isempty(f) || isempty(x) )
    fx = [];
    return
end

% Support for feval(f,'left') and feval(f,'end'), etc.
if ( ischar(x) )
    if ( any(strcmpi(x, {'left', 'start' ,'-'})) )
        fx = get(f, 'lval');
    elseif ( any(strcmpi(x, {'right', 'end', '+'})) )
        fx = get(f, 'rval');
    else
        error('CHEBFUN:feval:strinput', 'Unknown input argument "%s".', x);
    end
    return
end

% Reshape x to be a column vector.
sizex = size(x);
ndimsx = ndims(x);
m = size(f.impulses, 3);
x = x(:);

% Initialise:
fx = zeros(size(x, 1), m);
funs = f.funs;
nfuns = numel(funs);
dom = f.domain;

% Deal with feval(f,x,'left') and feval(f,x,'right'):
if ( nargin > 2 )
    lr = varargin{1};
    parse = strcmpi(lr, {'left', 'right', '-', '+'});
    if ( ~any(parse) )
        if ( ischar(lr) )
            error('CHEBFUN:feval:leftrightchar',...
                'Unknown input argument "%s".',lr);
        else
            error('CHEBFUN:feval:leftright','Unknown input argument.');
        end
    end
    % We deal with this by reassigning imps to be left/right values.
    if ( parse(1) || parse(3) ) % left
        f.impulses(2:end,:) = []; % Level 2 imps are not needed here
        for j = 1:nfuns
            f.impulses(1,j+1) = get(funs{j}, 'rval');
        end
    elseif parse(2) || parse(4) % right
        f.impulses(2:end,:) = []; % Level 2 imps are not needed here
        for j = 1:nfuns
            f.impulses(1,j) = get(funs{j}, 'lval');
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
for k = 1:nfuns
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

% Reshape if possible:
if ( ((ndimsx > 2) || (sizex(2) > 1)) && (m == 1) )
    fx = reshape(fx, sizex);

elseif ( ((ndimsx == 2) || (sizex(2) > 1)) && (m > 1))
    fx = reshape(fx, sizex(1), m*numel(x)/sizex(1));
end

if ( m > 1 )
    % [TODO]: Impulses won't work here!
   return
end

%% DEALING WITH IMPS:
% If the evaluation point corresponds to a breakpoint, we get the value from
% imps. If there is only one row, the value is given by the corresponding entry
% in that row. If the second row is nonzero the value is -inf or inf
% corresponding to the sign of the entry in the 2nd row. If the entry in the
% corresponding 3rd or higher rows is nonzero, we return NaN.

% Only one row:
if ( size(f.impulses,1) == 1 || ~any(any(f.impulses(2:end,:))) )
    % Loop over the funs:
    for k = 1:nfuns + 1
        fx( x == dom(k) ) = f.impulses(1,k);
    end
end


