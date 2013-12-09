function f = sqrt( f )
%SQRT   Square root of a chebfun2

if ( isempty( f ) ) % check for empty chebfun2.
    f = chebfun2;
    return
end 

%  positive/negative test. 
% [bol wzero] = singlesigntest(f); 
%
% if bol == 0 || wzero == 1
%    error('CHEBFUN2:SQRT','A change of sign/zero has been detected, unable to represent the result.'); 
% end

% Still call the constructor in case we missed a change of sign. 

op = @(x,y) sqrt( feval( f, x, y ) ); % resample.  
f = chebfun2( op, f.domain );         % Call constructor.

end