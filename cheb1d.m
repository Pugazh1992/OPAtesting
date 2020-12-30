function fout = cheb1d(coeff,x,xspan)
% function that fits chebyshev polynomial with given coefficients to
% input independent variable x (vector)

M = length(x);
N = length(coeff);
fout = zeros(1,M);


if nargin == 2
    for j=1:M
        t = 2*(x(j) - x(1)) / (x(end) - x(1)) - 1;
        %t = x(j);
        cm2 = 1; % the cm# corresponds to t = t-#. this is the first polynomial
        cm1 = t;
        fout(j) = coeff(2)*cm1 + 0.5*coeff(1)*cm2;

        for i=3:N
            c0 = 2*t*cm1 - cm2;
            fout(j) = fout(j) + coeff(i)*c0;
            cm2 = cm1;
            cm1 = c0;        
        end
    end
else
    for j=1:M
        t = 2*(x(j) - xspan(1)) / (xspan(end) - xspan(1)) - 1;
        %t = x(j);
        cm2 = 1; % the cm# corresponds to t = t-#. this is the first polynomial
        cm1 = t;
        fout(j) = coeff(2)*cm1 + 0.5*coeff(1)*cm2;

        for i=3:N
            c0 = 2*t*cm1 - cm2;
            fout(j) = fout(j) + coeff(i)*c0;
            cm2 = cm1;
            cm1 = c0;        
        end
    end
end

end