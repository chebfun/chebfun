function g = restrict(f, s)
%RESTRICT   Restrict a SINGFUN to subinterval(s).
%   RESTRICT(F, S) returns a ONEFUN that is restricted to the subinterval
%   [S(1),S(2)] of [-1,1]. Note that since SINGFUN objects only live on
%   [-1,1], a linear change of variables is implicitly applied.
%
%   If length(S) > 2, i.e., S = [S1, S2, S3, ...], then RESTRICT(F, S) returns
%   an cell of ONEFUN objects, where the entries hold F restricted to each of
%   the subintervals defined by S.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Check if s is actually a subinterval:
if ( (s(1) < -1) || (s(end) > 1) || (any(diff(s) <= 0)) )
    error('CHEBFUN:SINGFUN:restrict:badInterval', 'Not a valid interval.')
elseif ( (numel(s) == 2) && all(s == [-1, 1]) )
    % Nothing to do here!
    return
end

% Number of subintervals:
numInts = numel(s) - 1;

% Preallocate a cell:
g = cell(1, numInts);

% In the following construction process, we try to make the best of the
% information stored in the original SINGFUN f.

for j = 1:numInts
    % Check if the first subinterval includes the left endpoint of the original
    % domain (i.e., -1).
    
    if ( ( s(j) == -1 ) && ( f.exponents(1) ) )
        
        % Create an empty SINGFUN:
        gtmp = singfun();
        
        % Inherit the exponents:
        gtmp.exponents = [f.exponents(1), 0];
        
        % Construct the smoothPart of the restricted SINGFUN:
        if ( ~f.exponents(2) )
            
            % If the exponent of the right end is zero, it's an easy task:
            gtmp.smoothPart = ((1 + s(2))/2)^f.exponents(1)*...
                restrict(f.smoothPart, s(1:2));
            
        else
            
            % Otherwise:
            gtmp.smoothPart = ((1+s(2))/2)^f.exponents(1)*...
                restrict(f.smoothPart, s(1:2)).*...
                smoothfun.constructor(@(x) (1 - s(2)*(x + 1)/2 - ...
                    (x - 1)/2).^f.exponents(2));
            
        end
        
        % Put in cell:
        g{1} = gtmp;
        
        continue
    end
    
    % Check if the last subinterval includes the right endpoint of the original
    % domain (i.e., 1).
    
    if ( ( s(j+1) == 1 ) && ( f.exponents(2) ) )
        
        % Create an empty SINGFUN:
        gtmp = singfun();
        
        % Inherit the exponents:
        gtmp.exponents = [0, f.exponents(2)];
        
        % Construct the smoothPart of the restricted SINGFUN:
        if ( ~f.exponents(1) )

            % If the exponent of the left end is zero, it's an easy task:
            gtmp.smoothPart = ((1 - s(end-1))/2)^f.exponents(2)*...
                restrict(f.smoothPart, s(end-1:end));
            
        else
            
            % Otherwise:
            gtmp.smoothPart = ((1 - s(end-1))/2)^f.exponents(2)*...
                restrict(f.smoothPart, s(end-1:end)).*...
                smoothfun.constructor(@(x) (s(end-1)*(1 - x)/2 + ...
                    (x + 1)/2 + 1).^f.exponents(1));

        end
        
        % Put in cell:
        g{end} = gtmp;
        
        continue
    end
    
    % For the subintervals where the piece is just a smoothfun:
    
    % Define the new operator:
    op = @(x) feval(f, ((1 - x)*s(j) + (1 + x)*s(j+1))/2);
    
    % Call the smoothfun constructor
    data.vscale = [];
    data.hscale = 1;
    gtmp = smoothfun.constructor(op, data);
    
    % Put in cell:
    g{j} = gtmp;
    
end

% If there is only one subinterval, return the singfun, instead of a cell:
if ( numInts == 1 )
    g = g{:};
end

end
