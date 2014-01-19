function f = cumsum(f, k, dim, shift)
%CUMSUM   Indefinite integral of an UNBNDFUN.
%   CUMSUM(F) is the indefinite integral of the UNBNDFUN F on an interval [a,b],
%   with the constant of integration chosen so that F(a) = 0.
%
%   CUMSUM(F, K) will compute the Kth indefinite integral with the constant of
%   integration chosen so that each intermediate integral evaluates to 0 at x=a.
%   Thus CUMSUM(F, 2) is equivalent to CUMSUM(CUMSUM(F)).
%
%   CUMSUM(F, K, 2) will take the Kth cumulative sum over the columns F an
%   array-valued UNBNDFUN.
%
%   CUMSUM(F, K, DIM, S) will shift F up by S. Note that this could be useful at
%   the CHEBFUN level to concatenate different pieces forming a countinuous
%   object.
%
% See also DIFF, SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%
% Trivial case of an empty UNBNDFUN:
if ( isempty(f) )
    return
end

% Parse inputs:
if ( nargin == 1 || isempty(k) )
    % Compute first indefinite intergral by default
    k = 1;
end

if ( nargin < 3 )
    % Assume dim = 1 by default
    dim = 1;
end

if ( nargin < 4 )
    % Assume no need to shift:
    shift = 0;
end

if ( dim == 1 )
    
    % Rescaling factor is the derivative of the forward map:
    pref = chebpref();
    pref.singPrefs.exponents = f.mapping.forderExps;
    rescaleFactor = onefun.constructor(@(x) f.mapping.forder(x), ...
        [], [], pref);
    
    % Check if we need to add breakpoints:
    f = addBreaksForCumSum(f);
    
    if ( iscell(f) )
        
        % If f is a cell of two UNBNDFUNs which are obtained by adding a new
        % breakpoint at the midpoint of the original domain due to non-trivial
        % exponents at both endpoints, we integrate each UNBNDFUN individually.
        
        for j = 1:2
            
            % Loop K times for the Kth integral:
            for l = 1:k
                
                % The ONEFUN of the integral of F should be the integral of the
                % ONEFUN of the F multiplied by the derivative of the forward
                % map:
                f{j}.onefun = f{j}.onefun.*rescaleFactor;
                f{j}.onefun = cumsum(f{j}.onefun);
            end
            
        end
        
    else
        
        % Loop K times for the Kth integral:
        for l = 1:k
            
            % The ONEFUN of the integral of F should be the integral of the
            % ONEFUN of the F multiplied by the derivative of the forward
            % map:
            f.onefun = f.onefun.*rescaleFactor;
            f.onefun = cumsum(f.onefun);
        end
        
        
        % Shift F up or down. This is useful at the chebfun level to concatenate the
        % piece making the entire function as continuous as possible.
        if ( ~any( issing(f) ) )
            
            % Grab the indice correspond to infinite shift:
            ind = isinf(shift);
            
            % Zero the infinite shift:
            shift( ind ) = 0;
            
            % Shift:
            f = f + shift - get(f, 'lval');
        end
        
    end
    
elseif ( dim == 2 )
    
    % When the third argument is 2, i.e. dim = 2, we compute the cumlative sum
    % over columns, in which case, no rescale is needed.
    
    f.onefun = cumsum(f.onefun, k, dim);
    
else
    error('CHEBFUN:UNBNDFUN:cumsum:input', ...
        'The third argument is unrecognizable.');
end

end