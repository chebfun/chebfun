function F = diff(F, k, dim)
%DIFF   Partial derivative of a CHEBFUN3.
%   DIFF(F) is the 1st derivative of F along the 1st input variable.
%
%   DIFF(F, K) is the K-th derivative of F along the 1st input variable.
%  
%   DIFF(F, K, DIM) is the K-th derivative of F along the dimension DIM.
%   DIM = 1 (default) returns the derivative in the 1st input variable.
%   DIM = 2 returns the derivative in the 2nd variable.
%   DIM = 3 returns the derivative in the 3rd varialbe.
%
%   DIFF(F, [K1 K2], [DIM1 DIM2]) means K1-th derivative of F in dimension
%   DIM1 and K2-th derivative in dimension DIM2. DIM1 and DIM2 can be 1, 2 
%   or 3 in any order.
%
%   DIFF(F, [K1 K2 K3]) is the K1-th partial derivative of F in the first 
%   variable, K2-th partial derivative of F in the second variable and 
%   K3-th partial derivative of F in the third variable. 
%   For example, DIFF(F, [1 2 3]) is d^6F/(dx d^2y d^3z).
%    
% See also CHEBFUN3/GRAD, CHEBFUN3/LAP, and CHEBFUN3/BIHARM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) )
    return
end

% Default to first derivative in the 1st variable x:
if ( nargin == 1 )    
    k = 1;
    dim = 1;
end

% Default to k-th derivative in the 1st variable:
if ( nargin == 2 )
    if ( numel(k) == 1 )    
        dim = 1;
    end
elseif ( nargin == 3 && numel(dim) ~= 1 && numel(dim) ~= 2 )
    error('CHEBFUN:CHEBFUN3:diff:dim', 'Dim should have either 1 or 2 entries.');
    % Developer Note: If K has three entries, we always assume that dim is 
    % [1 2 3]. So, we don't accept dim with 3 indices.
end

% Diff the individual column, row, and tubes:
if ( (numel(k) == 3) && (nargin < 3) )
    F.cols =  diff(F.cols, k(1));
    F.rows =  diff(F.rows, k(2));
    F.tubes =  diff(F.tubes, k(3));
    
elseif ( (numel(k) == 2) && (nargin == 3) && (numel(k) == numel(dim)) )
    [is1,loc1] = ismember(1, dim);
    [is2,loc2] = ismember(2, dim);
    [is3,loc3] = ismember(3, dim);
    if ( is1 )
        % k1-th diff in the 1st variable
        k1 = k(loc1);
        F.cols = diff(F.cols, k1);
    end
    if ( is2 ) 
        % k2-th diff in the 2nd variable.
        k2 = k(loc2);
        F.rows = diff(F.rows, k2);
    end
    if ( is3 )
        % k3-th diff in the 3rd variable.
        k3 = k(loc3);
        F.tubes = diff(F.tubes, k3);
    end
    
elseif ( dim == 1 ) 
    % numel(k) = 1 & dim = 1. So, we want k-th diff in the 
    % 1st variable.
    F.cols = diff(F.cols, k);
    
elseif ( dim == 2 ) 
    % numel(k)=1 & dim = 2. So, k-th diff in the 
    % y-direction. Call Chebfun2/diff.
    F.rows = diff(F.rows, k);
    
elseif ( dim == 3 )
    % numel(k)=1 & dim=3. So, k-th diff in the z-direction. 
    % Call Chebfun2/diff.
    F.tubes = diff(F.tubes, k);

else
    % dim is not 1, 2 or 3. Or, numel(k) ~= numel(dim)
    error('CHEBFUN:CHEBFUN3s:diff:dim', ...
          'Can compute derivative in x, y or z only.');
end

end