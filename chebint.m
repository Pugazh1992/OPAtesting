function fout = chebint(coeff)
% function that integrates chebyshev series and provides integrated series
% as output. be sure to scale the output properly when evaluating ( (tf-t0)/2 ),
% and to use appropriate boundary condition in the constant.
N = length(coeff);
fout = zeros(N,1);
n = N - 1;
fout(1) = 0;
% for k=2:n
%     fout(k) = ( coeff(k-1) - coeff(k+1) ) / ( 2.0 * (k-1) );
% end

for i=1:n-1
    fout(i+1) = ( coeff(i+1-1) - coeff(i+1+1) ) / ( 2.0 * (i) );
end

fout(n+1) = coeff(n+1-1) / ( 2.0 * n );

end