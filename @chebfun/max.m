function [y, x] = max(f, flag)
%MAX    Maximum value of a CHEBFUN.
%   MAX(F) returns the maximum value of the CHEBFUN F.
%
%   [Y, X] = MAX(F) returns also point X such that F(X) = Y.
%
%   [Y, X] = MAX(F, 'local') returns not just the global maximum value, but all
%   of the local maxima.
%
%   If F is complex-valued, absolute values are taken to determine the maxima,
%   but the resulting values correspond to those of the original function.
%
%   If F is array-valued, then the columns of X and Y correspond to the columns
%   of F. NaNs are used to pad Y and X when the 'local' flag is used and the
%   columns are not of the same length
%
% See also MIN, MINANDMAX, ROOTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL MAXIMUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( nargin == 1 ) 
    
    % Call MINANDMAX():
    [y, x] = minandmax(f);
    
    % Extract the minimum:
    y = y(2,:);
    x = x(2,:);
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL MINIMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    % Call MINANDMAX():
    [y, x] = minandmax(f, flag);

    % Determine which are maxima.
    
    ends = f.domain([1, end]).'; % Endpoints of the domain are special.
    f = mat2cell(f); % Convert f into a cell of scalar-valued CHEBFUNs.
    
    % Loop over the FUNs:
    for k = 1:numel(f)
        % For interior extrema, look at 2nd derivative:
        idx = feval(diff(f{k}, 2), x(:,k)) < 0;
        % For end-points, look at 1st derivative:
        idx2 = feval(diff(f{k}), ends).*[1, -1]' < 0;
        idx(1) = idx2(1);
        idx(x(:,k) == ends(2)) = idx2(2);
        % Set points corresponding to local minima to NaN:
        y(~idx,k) = NaN;
        x(~idx,k) = NaN;
        % Sort the result
        [x(:,k), idx] = sort(x(:,k));
        y(:,k) = y(idx,k);
    end
    
    % Remove any rows which contain only NaNs.
    x(all(isnan(x), 2),:) = []; 
    y(all(isnan(y), 2),:) = [];

end

end