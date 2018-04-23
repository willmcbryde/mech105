
function[root, fx, ea, iter] = FalsePosition(func, xl, xu, es, maxiter)
tic
format long
%falseposition: This function uses the false position algorithm to find the roots of any
%function that the user inputs
%func = the function being evaluated 
%xl = the lower bracket
%xu = the upper bracket of the root
%es = the desired relative error (will default to 0.0001)
%maxiter = number of iterations (will default to 200)
%root = estimated root location
%fx = function value at the root
%ea = approximate relative error
%iter = number of iterations performed 

if nargin < 3
    error('At least 3 input arguments are required');
    %Make sure enough arguments are input
end 
if isempty('es') == 0
    es = 0.0001;
    %Default es to 0.0001 if user does not provide input
end
if isempty('maxiter') == 0
    maxiter = 200;
    %Default maxiter to 200 if user does not provide input
end 
%Define all known starting variables
root = 0;
ea = 100;
iter = 100;
xr = 0;
%The false position function is now run 
while ea > es && iter < maxiter
    xr = root;
    if (func(xl)*func(xu)) > 0
        error('Values for xl and xu are out of the bounds of the root; choose new bounds')
    end 
    root = xu - ((func(xu)*(xl-xu))/(func(xl)-func(xu)));
    if (func(root) == inf)
        fprintf('The root for this function is %f\n', root)
        fprintf('The function evaluated at root is %f\n', 0)
        fprintf('The approximate relative error is %f\n', 0)
        fprintf('The number of iterations performed was %0.0f\n',1)
        return 
        %If the right root is found, end the function
    end 
    if (func(xl)*func(root)) < 0
        xu = root;
    else 
        xl = root;
    end 
    if iter > 0
        ea = (abs((xr-root)/(root)))*100;
    end 
    iter = iter + 1;
end 
fx = func(root);
fprintf('The root is %f\n', root)
fprintf('The function evaluated at the root is %f\n', fx)
fprintf('The approximate relative error is %f\n', ea)
fprintf('The number of iterations performed was %0.0f\n', iter)
toc
end 
