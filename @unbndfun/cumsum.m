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
    % Compute first indefinite intergral by default.
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
    
    % Make a copy of F:
    g = f;
    
    % Grab the preference:
    pref = chebpref();
    
    % Loop K times:
    for j = 1:k
        
        if ( iscell(g) )
            
            for l = 1:2
                % Rescaling factor is the derivative of the forward map:
                pref.singPrefs.exponents = g{l}.mapping.forderExps;
                rescaleFactor = onefun.constructor( ...
                    @(x) g{l}.mapping.forder(x), [], [], pref);
                numRoots = -repmat(pref.singPrefs.exponents, 1, size(g{l}, 2));
                
                % Try to see if we can extract boundary roots:
                [h, rootsLeft, rootsRight] = extractBoundaryRoots(g.onefun, ...
                    numRoots);
                
                if ( all( rootsLeft == numRoots(1,:) ) && ...
                        all( rootsRight == numRoots(2,:) ) )
                    
                    g{l}.onefun = h;
                    % The ONEFUN of the integral of F should be the integral of the
                    % ONEFUN of the F multiplied by the derivative of the forward
                    % map. Here the singularities of the RESCALEFACTOR is cancelled off by
                    % the boundary roots of H. Therefore, only the smoothPart of
                    % RESCALEFACTOR is involved.
                    g{l}.onefun = cumsum(g{l}.onefun.*rescaleFactor.smoothPart);
                    
                else
                    
                    % The ONEFUN of the integral of F should be the integral of the
                    % ONEFUN of the F multiplied by the derivative of the forward
                    % map.
                    g{l}.onefun = g{l}.onefun.*rescaleFactor;
                    g{l}.onefun = cumsum(g{l}.onefun);
                    
                end
                
            end
            
        else
            
            % Rescaling factor is the derivative of the forward map:
            pref.singPrefs.exponents = g.mapping.forderExps;
            rescaleFactor = onefun.constructor(@(x) g.mapping.forder(x), ...
                [], [], pref);
            numRoots = -repmat(pref.singPrefs.exponents.', 1, size(g, 2));
            
            % Try to see if we can extract boundary roots:
            [h, rootsLeft, rootsRight] = extractBoundaryRoots(g.onefun, ...
                numRoots);
            
            if ( all( rootsLeft == numRoots(1,:) ) && ...
                    all( rootsRight == numRoots(2,:) ) )
                
                % The ONEFUN of the integral of F should be the integral of the
                % ONEFUN of the F multiplied by the derivative of the forward
                % map. Here the singularities of the RESCALEFACTOR is cancelled off by
                % the boundary roots of H. Therefore, only the smoothPart of
                % RESCALEFACTOR is involved.
                g.onefun = h;
                g.onefun = cumsum(g.onefun.*rescaleFactor.smoothPart);
                
            else
                
                % The ONEFUN of the integral of F should be the integral of the
                % ONEFUN of the F multiplied by the derivative of the forward
                % map.
                g.onefun = g.onefun.*rescaleFactor;
                
                % Check if we need to add breakpoints:
                g = addBreaksForCumSum(g);
                
                if ( iscell(g) )
                    g{1}.onefun = cumsum(g{1}.onefun);
                    g{2}.onefun = cumsum(g{2}.onefun);
                else
                    g.onefun = cumsum(g.onefun);
                end
                
            end
        end
        
    end
    
    f = g;
    
    % Shift F up or down. This is useful at the chebfun level to concatenate the
    % piece making the entire function as continuous as possible.
    if ( ~any( issing(f) ) && ( ~iscell(f) ) )
        
        % Grab the indice correspond to infinite shift:
        ind = isinf(shift);
        
        % Zero the infinite shift:
        shift( ind ) = 0;
        
        % Shift:
        f = f + shift - get(f, 'lval');
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