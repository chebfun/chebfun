function vals = feval(f,x,y,varargin)
%FEVAL pointwise evaluate a chebfun2v.
%
%  F(X,Y) returns the evaluation of F at the coordinate (X,Y).
%
% See also SUBSREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.


% Should this be length(x) by 2 or 2 by length(x)

% Make it length(x) by 2 from now.

if isempty(f.zcheb)
    vals = nan(2,length(x));
    vals(1,:) = feval(f.xcheb,x,y,varargin);
    vals(2,:) = feval(f.ycheb,x,y,varargin);
else
    vals = nan(3,length(x));
    vals(1,:) = feval(f.xcheb,x,y,varargin);
    vals(2,:) = feval(f.ycheb,x,y,varargin);
    vals(3,:) = feval(f.zcheb,x,y,varargin);
end

end