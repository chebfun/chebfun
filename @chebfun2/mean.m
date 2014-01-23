function s = mean(f, varargin)
%MEAN   Average or mean value of a chebfun2. 
% 
%  MEAN(F) takes the mean in the y-direction (default), i.e., 
% 
%          MEAN(F) = 1/(ymax-ymin) sum(F).
%
%  MEAN(F,DIM) takes the mean along the direction DIM. If DIM = 1 it is the
%  y-direction and if DIM = 2 then it is the x-direction. 
%
% See also MEAN2, STD2.

if ( isempty( f ) ) % check for empty chebfun2.
    s = chebfun;
    return
end 

if ( nargin == 1) 
    dim = 1;    % default to the y-direction
else
    dim = varargin{ 1 };  
end

dom = f.domain;

s = sum( f, dim ); 
if ( dim == 1 )
    s = s / diff( dom(3:4) );            % mean in the y direction (default)
elseif ( dim == 2 )
    s = s / diff( dom(1:2) );            % mean in the x direction
else
    error('CHEBFUN2:MEAN:DIM','Mean not in x or y direction')
end

end
