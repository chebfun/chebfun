function f = classifyExponents(f)
%CLASSIFYEXPONENTS   Function to assign types of exponents in a SINGFUN object.
%   Based on the values in F.EXPONENTS, this functions decides the type that 
%   should be assigned to F.SINGTYPE. The valid types can be 'sing', 'pole', 
%   'root' or 'none'. F.SINGTYPE is a 1X2 cell array and the pair of string
%   contained in this field describes the types of singularities at each end,
%   -1 or 1 of the SINGFUN F. These types have the following meaning:
%    
%      'pole' - A pole, i.e. a negative integer exponent at the 
%               corresponding end.
%      'sing' - A negative real exponent at the corresponding end. 
%               Can be an integer as well.
%      'root' - A root of fractional order at the corresponding end point.
%      'nont' - No singularity at the end point.

%
% See also CHECKSINGTYPES

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Get the SINGFUN tolerance
% [TODO]: This should depend on scales, but what are the scales?
%         This TODO will be setteled once scales are finalised.
tol = singfun.pref.singfun.eps;

%%
% Store the exponents in exps (for brevity):
exps = f.exponents;
% Loop on the left and right end point of the domain
for k = 1:2
    % If positive exponent
    if ( exps(k) >= 0 )
        if ( abs(exps(k) - round(exps(k))) < 100*tol )
            % Positive integer exponent, i.e. no singularity
            f.singType{k} = 'none';
        else
            % The function is bounded but there is a root of fractional order.
            f.singType{k} = 'root';
        end
    else        
        % Negative exponents        
        if ( exps(k) > -100*tol )
            % The exponent is negative but almost zero,
            % remove the singularity
            f.singType{k} = 'none';
        else
            % Non-trivial negative exponent
            if ( abs(exps(k) - round(exps(k))) < 100*tol )
                % Pole if integer valued exponent
                f.singType{k} = 'pole';
            else
                % A fractional pole, which we call 'sing'
                f.singType{k} = 'sing';
            end
        end
    end
end
end