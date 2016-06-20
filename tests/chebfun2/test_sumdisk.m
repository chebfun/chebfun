function pass = test_sumdisk( pref ) 
% Test for sumdisk command of a chebfun2 object. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 1000*eps;
j = 1; 
pass = [];

% Load set of test functions

ListLength = 18;

TestFn_Cell = cell(ListLength,3);           % 1st column function name, 2nd colnum function itself, 3rd column is a non-defult domain


TestFn_Cell{1,1} = 'Exp fn';
TestFn_Cell{1,2} = @(x,y) exp(-30.*(x.^2+y.^2));
TestFn_Cell{1,3} = [0 3 0 3];

TestFn_Cell{2,1} = 'Bessel fn';
TestFn_Cell{2,2} = @(x,y) besselj(0,10*sqrt(x.^2+y.^2));
TestFn_Cell{2,3} = [-2 0 -2 0];

TestFn_Cell{3,1} = 'Cos fn 1';
TestFn_Cell{3,2} = @(x,y) .3*cos(2*(x.^2+y.^2)/sqrt(2));
TestFn_Cell{3,3} = [-2 2 -2 2];

TestFn_Cell{4,1} = 'Runge fn';
TestFn_Cell{4,2} = @(x,y) 1./(0.1+4*(x.^2+y.^2));
TestFn_Cell{4,3} = [1 2 3 4];

TestFn_Cell{5,1} = 'Cos fn 2';
TestFn_Cell{5,2} = @(x,y) (1.25+cos(5*y))./(6+6*(3*x-1).^2);
TestFn_Cell{5,3} = [-1 2 -2 1];

TestFn_Cell{6,1} = 'test fn 1';
TestFn_Cell{6,2} = @(x,y) (1/3)*cos(x.*2 - y.^2);   
TestFn_Cell{6,3} = [1 2 1 2];

TestFn_Cell{7,1} = 'test fn 2';
TestFn_Cell{7,2} = @(x,y) (1/3)*sin(x - y) .* cos(x + y);
TestFn_Cell{7,3} = [1 3 0.5 1];

TestFn_Cell{8,1} = 'test fn 3';
TestFn_Cell{8,2} = @(x,y) x.^2 + y - 0.5;
TestFn_Cell{8,3} = [-0.5 0 -3 0];

TestFn_Cell{9,1} = 'test fn 4';
TestFn_Cell{9,2} = @(x,y) 2*cos(10*x).*sin(10*y) + sin(10*x.*y);
TestFn_Cell{9,3} = [2 4 2 4];

TestFn_Cell{10,1} = 'even fn';
TestFn_Cell{10,2} = @(x,y) exp(-50*x.^2) + 0.75*exp(-50*y.^2) +0.75*exp(-50*x.^2).*exp(-50*y.^2);
TestFn_Cell{10,3} = [-2 4 -2 4];

TestFn_Cell{11,1} = 'odd fn';
TestFn_Cell{11,2} = @(x,y) (sin(x) + sin(y)).*(exp(x.^2));
TestFn_Cell{11,3} = [0 1 -1 5];

TestFn_Cell{12,1} = 'const fn';
TestFn_Cell{12,2} = @(x,y) 1;
TestFn_Cell{12,3} = [-3 3 -1 4];

TestFn_Cell{13,1} = 'linear fn';
TestFn_Cell{13,2} = @(x,y) x + y;
TestFn_Cell{13,3} = [-3 0 0 3];

TestFn_Cell{14,1} = 'quadratic fn 3';
TestFn_Cell{14,2} = @(x,y) x.^2 + y.^2;
TestFn_Cell{14,3} = [-0.5 0 -3 0];

TestFn_Cell{15,1} = 'trigfun fn 1';
TestFn_Cell{15,2} = @(x,y) cos(pi*x) + sin(pi*y);
TestFn_Cell{15,3} = [1 3 1 3];

TestFn_Cell{16,1} = 'trigfun fn 2';
TestFn_Cell{16,2} = @(x,y) cos(pi*x).*sin(pi*y);
TestFn_Cell{16,3} = [0 6 0 4];

TestFn_Cell{17,1} = 'trigfun fn 3';
TestFn_Cell{17,2} = @(x,y) (cos(2*pi*x)).^2 + sin(pi*y);
TestFn_Cell{17,3} = [1.5 3.5 -1 3];

TestFn_Cell{18,1} = 'trigfun fn 4';
TestFn_Cell{18,2} = @(x,y) 1./(sin(2*pi*x) + 2);
TestFn_Cell{18,3} = [10 12 -2 2];


% display the test functions
figure(1)
for I = 1:ListLength
    disp(['Display Test function' num2str(I)])
    f = TestFn_Cell{I,2};
    f_Chebfun2 = chebfun2(f,'vectorize');
    subplot(3,6,I)
    plot(f_Chebfun2)
    title(TestFn_Cell{I,1},'FontSize',10)
end

% Test sumdisk using chebfun2 test functions defined on the unit disk
disp('Test sumdisk using chebfun2 test functions defined on the unit disk')

Int_Err_Vec = zeros(14,1);

for I = 1:14

    disp(['Test fn ' num2str(I)])
    f = TestFn_Cell{I,2};
    polarfun = @(theta,r) f(r.*cos(theta),r.*sin(theta)).*r;
    Integral_Disc_integral2 = integral2(polarfun,0,2*pi,0,1,'AbsTol',1e-15,'RelTol',1e-15);
    disp(['integral2 gives I = ' num2str(Integral_Disc_integral2,20)])
    
    f_Chebfun2 = chebfun2(f, [-1 1 -1 1], 'vectorize');
    Integral_Disc_sumdisk = sumdisk(f_Chebfun2);
    disp(['sumdisk gives I = ' num2str(Integral_Disc_sumdisk,20)])
    
    Int_Err_Vec(I) = Integral_Disc_integral2 - Integral_Disc_sumdisk;
    disp(['Err = ' num2str(Int_Err_Vec(I))])
    pass(j) = abs(Int_Err_Vec(I) < tol);
    j = j+1;
end

Max_Err = max(abs(Int_Err_Vec));
disp(['Maximum error among all test functions (trigfun2 on unit disk) = ' num2str(Max_Err,20)])

% Test sumdisk using trigfun2 test functions defined on the unit disk
disp('Test sumdisk using trigfun2 test functions defined on the unit disk')

Int_Err_Vec = zeros(14,1);

for I = 15:18
    disp(['Test fn ' num2str(I)])
    f = TestFn_Cell{I,2};
    polarfun = @(theta,r) f(r.*cos(theta),r.*sin(theta)).*r;
    Integral_Disc_integral2 = integral2(polarfun,0,2*pi,0,1,'AbsTol',1e-15,'RelTol',1e-15);
    disp(['integral2 gives I = ' num2str(Integral_Disc_integral2,20)])
    
    f_Chebfun2 = chebfun2(f, [-1 1 -1 1], 'vectorize','trig');
    Integral_Disc_sumdisk = sumdisk(f_Chebfun2);
    disp(['sumdisk gives I = ' num2str(Integral_Disc_sumdisk,20)])
    
    Int_Err_Vec(I) = Integral_Disc_integral2 - Integral_Disc_sumdisk;
    disp(['Err = ' num2str(Int_Err_Vec(I))])
    pass(j) = abs(Int_Err_Vec(I) < tol);
    j = j+1;
end

Max_Err = max(abs(Int_Err_Vec));
disp(['Maximum error among all test functions (trigfun2 on unit disk) = ' num2str(Max_Err,20)])


% Test sumdisk using chebfun2 test functions defined on the non-default
% domain
disp('Test sumdisk using chebfun2 test functions defined on non-default domain')

Int_Err_Vec = zeros(14,1);

for I = 1:14

    disp(['Test fn ' num2str(I)])
    fdomain = TestFn_Cell{I,3};
    disp(['Domain [' num2str(fdomain) ']'])

    f = TestFn_Cell{I,2};
    g = @(x,y) f(((fdomain(2) - fdomain(1))/2)*x + (fdomain(2) + fdomain(1))/2,((fdomain(4) - fdomain(3))/2)*y + (fdomain(4) + fdomain(3))/2);
    polarfun = @(theta,r) g(r.*cos(theta),r.*sin(theta)).*r;
    Integral_Disc_integral2 = integral2(polarfun,0,2*pi,0,1,'AbsTol',1e-15,'RelTol',1e-15);
    Integral_Disc_integral2 = Integral_Disc_integral2 * (fdomain(2) - fdomain(1)) * (fdomain(4) - fdomain(3))/4;
    disp(['integral2 gives I = ' num2str(Integral_Disc_integral2,20)])
    
    f_Chebfun2 = chebfun2(f, fdomain, 'vectorize');
    Integral_Disc_sumdisk = sumdisk(f_Chebfun2);
    disp(['sumdisk gives I = ' num2str(Integral_Disc_sumdisk,20)])
    
    Int_Err_Vec(I) = Integral_Disc_integral2 - Integral_Disc_sumdisk;
    disp(['Err = ' num2str(Int_Err_Vec(I))])
    pass(j) = abs(Int_Err_Vec(I) < tol);
    j = j+1;
end

Max_Err = max(abs(Int_Err_Vec));
disp(['Maximum error among all test functions (trigfun2 on unit disk) = ' num2str(Max_Err,20)])
 
% Test sumdisk using trigfun2 test functions defined on the non-default
% domain
disp('Test sumdisk using trigfun2 test functions defined on non-default domain')

Int_Err_Vec = zeros(14,1);

for I = 15:18

    disp(['Test fn ' num2str(I)])
    fdomain = TestFn_Cell{I,3};
    disp(['Domain [' num2str(fdomain) ']'])
    
    f = TestFn_Cell{I,2};
    g = @(x,y) f(((fdomain(2) - fdomain(1))/2)*x + (fdomain(2) + fdomain(1))/2,((fdomain(4) - fdomain(3))/2)*y + (fdomain(4) + fdomain(3))/2);
    polarfun = @(theta,r) g(r.*cos(theta),r.*sin(theta)).*r;
    Integral_Disc_integral2 = integral2(polarfun,0,2*pi,0,1,'AbsTol',1e-15,'RelTol',1e-15);
    Integral_Disc_integral2 = Integral_Disc_integral2 * (fdomain(2) - fdomain(1)) * (fdomain(4) - fdomain(3))/4;
    disp(['integral2 gives I = ' num2str(Integral_Disc_integral2,20)])
    
    f_Chebfun2 = chebfun2(f, fdomain, 'trig', 'vectorize');
    Integral_Disc_sumdisk = sumdisk(f_Chebfun2);
    disp(['sumdisk gives I = ' num2str(Integral_Disc_sumdisk,20)])
    
    Int_Err_Vec(I) = Integral_Disc_integral2 - Integral_Disc_sumdisk;
    disp(['Err = ' num2str(Int_Err_Vec(I))])
    pass(j) = abs(Int_Err_Vec(I) < tol);
    j = j+1;
end

Max_Err = max(abs(Int_Err_Vec));
disp(['Maximum error among all test functions (trigfun2 on non-defualt domain) = ' num2str(Max_Err,20)])
 


end