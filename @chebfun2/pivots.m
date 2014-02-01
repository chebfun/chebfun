function p = pivots(f, varargin)
%PIVOTS  pivot values of a chebfun2
% 
% PIVOTS(F) returns the pivot values taken during in the constructor by the 
%  GE algorithm. 
% 
% PIVOTS(F, 'normalize'), returns the normalised pivot values.  These
% numbers are scaled so that the columns and rows have unit 2-norm.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin < 2 )
    
    p = f.pivotValues;  % Return the pivot values. 
    p = p(:);            % make column vector
    
elseif ( strcmpi(varargin{1},'normalise') ... 
                               || strcmpi(varargin{1},'normalize') )
                           
    cscl = norm( f.cols );   
    rscl = norm( f.rows );
    p = f.pivotValues; 
    p = p.*cscl.*rscl; % normalized pivots. 
    p = p(:);   % make column vector.
    
else
    
    error('CHEBFUN2:PIVOTS:InPuts','Unrecognised second argument');
    
end
end