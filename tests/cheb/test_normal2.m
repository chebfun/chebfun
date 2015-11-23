function pass = test_normal2(~)

% Test that an indefinite covariance matrix gives an error.
Sigma = [1 2; 2 1];
pass(1) = false;
try
    cheb.normal2([0,0], Sigma);
catch ME
    if strcmpi(ME.identifier, 'CHEB:NORMAL2:covariance:nonSymPosDef')
        pass(1) = true;
    end
end

% Test that a nonsymmetric matrix gives an error.
Sigma = [3 1; 2 3];
pass(2) = false;
try
    cheb.normal2([0,0], Sigma);
catch ME
    if strcmpi(ME.identifier, 'CHEB:NORMAL2:covariance:nonSymPosDef')
        pass(2) = true;
    end
end

% Test that for this arbitrary bivariate gaussian the integral is close to 1.
Sigma = [4 -1.2; -1.2 4];
p = cheb.normal2([2 6], Sigma);
I = sum2(p);
pass(3) = abs(I-1) < 1e-13;

end
