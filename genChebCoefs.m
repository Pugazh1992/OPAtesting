function coeff = genChebCoefs(input,tspan,N,M,dis)
%   function to get coefficients for chebyshev polynomial approximation.
%   this matches format of junkins paper on orthogonal polynomial
%   approximation page 4
%   f - either a function pointer or sampled function
%   tspan - sample times
%   N - Order of polynomial to fit
%   M - number of input points
%   dis - 1 for discrete samples, 0 for functional input

% for dicrete input use cosine samples

    y = zeros(M+1,1); % function evaluated at cosine nodes
    x = y;  % for debugging
    coeff = zeros(N+1,1);

    tmult = 0.5*(tspan(end)-tspan(1)); % ?

    Cmat = zeros(N+1,M+1);

    % populate "A" matrix as in junkins paper. Cmat == A in paper
    for j = 1 : (M+1)
        % convert x,y input to chebyshev space
        if dis==0 % functional input
            %xj = cos((M+1-j)*pi/M);
            %y(j) = f(tmult*(1+xj)+tspan(1)); % convert back and evaluate
        else % discrete sample input
            xj = (tspan(j)-tspan(1))/tmult -1; % convert to cheby space
%             y(j) = input(j);
%               xj = tspan(j); % these inputs are already converted to -1 to 1
              y(j) = input(j);
        end

        % initialize first two polynomials in table. we also move the 2/M
        % scaling to here.
        Cmat(1,j) = 1 * (2/M);
        Cmat(2,j) = xj * (2/M);

        % find remaining polynomials using recurrence relation
        for i = 3:(N+1)
            Cmat(i,j) = 2*xj*Cmat(i-1,j) - Cmat(i-2,j);
        end
        x(j) = xj;
    end

    % clean up scaling
    % for j=1:M+1
    %     %Cmat(1,j) = 0.5*Cmat(1,j); % not everything needed to be multiplied by 2. Some points will end up being multiplied by 1/2
    % end

    for i=1:N+1
        Cmat(i,1) = 0.5*Cmat(i,1);
        Cmat(i,M+1) = 0.5*Cmat(i,M+1);
    end

    % find coefficients
    for i=1:N+1
        for j=1:M+1
            coeff(i) = coeff(i) + Cmat(i,j)*y(j);
        end
    end

    % this is the only difference in the coefficients if interpolation v.
    % least squares
    if M>N
        coeff(N+1) = 0.5*coeff(N+1);
    end

    %coeff = Cmat*y;

end
