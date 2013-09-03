function f = classifyExponents(f)
%CLASSIFYEXPONENTS   Function to assign types of exponents in a SINGFUN object.
%   Based on the values in F.EXPONENTS, this functions decides the type that 
%   should be assigned to F.SINGTYPE. The valid types can be 'sing', 'pole', 
%   'root' or 'none'.
%
% See also CHECKSINGTYPES

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Get the SINGFUN tolearnce
tol = singfun.pref.singfun.eps;

%%
% Loop on the left and right end point of the domain
for k = 1:2
    % if positive exponent
    if ( f.exponents(k) >= 0 )
        if ( abs(f.exponents(k) - round(f.exponents(k))) < 100*tol )
            % positive integer exponent, i.e. no singularity
            f.singType{k} = 'none';
        else
            % the function is bounded but there is a root of fractional order.
            f.singType{k} = 'root';
        end
    else
        % the exponents are negative
        if ( f.exponents(k) > -100*tol )
            % the exponent is negative but almost zero,
            % remove the singularity
            f.singType{k} = 'none';
        else
            % non-trivial negative exponent
            if ( abs(f.exponents(k) - round(f.exponents(k))) < 100*tol )
                % pole if integer valued exponent
                f.singType{k} = 'pole';
            else
                % a fractional pole, which we call 'sing'
                f.singType{k} = 'sing';
            end
        end
    end
end

end
