function v = integral( F, c )
%INTEGRAL line integration of a chebfun2v
%
% INTEGRAL(F,C) computes the line integral of F along the curve C, that is
%
%                   
%                  /
% INTEGRAL(F,C) =  |  < F(r), dr > 
%                 /
%                 C 
%
% where the curve C is parameterised by the complex curve r(t).  

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

if ( F.nComponent == 3 )
    warning('CHEBFUN2V:INTEGRAL','Ignoring third component of chebfun2v.')
end

% restrict to chebfun: 
F1 = restrict(F.component{1}, c); 
F2 = restrict(F.component{2}, c); 
 
% Line integral: 
dc = diff(c); 
v = sum(F1 .* real(dc) + F2 .* imag(dc) );  

end