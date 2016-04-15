function v = std3(f)
%STD3   Standard deviation of a CHEBFUN3.
%   V = STD3(F) computes the standard deviation of a CHEBFUN3, i.e., 
%
%     STD3(F)^2 = 1/V*sum3(|F(x,y,z) - m|^2)
%
%   where V is the volume of the domain of F and m is the mean of F.
%
% See also MEAN, MEAN2, MEAN3, STD and STD2.

h = f - mean3(f);
v = sqrt(mean3(h .* conj(h)));

end