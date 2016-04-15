function v = mean3(f)
%MEAN3   Mean of a CHEBFUN3.
%   V = MEAN3(F) returns the mean of a CHEBFUN3, i.e., sum3(F) / V,
% 	where V is the volume of the domain of F.
%
%   See also MEAN, MEAN2 and STD3.

% Empty check:
if ( isempty(f) ) 
    return
end

% Apply the formula: 
v = sum3(f) / domainvolume(f);  

end