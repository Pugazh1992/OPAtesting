function fapprox = cheb2d( coeff , x , y , Nx , Ny , xbound , ybound )
   

xi  = 2*(x-xbound(1))/(xbound(2)-xbound(1)) - 1 ;
eta = 2*(y-ybound(1))/(ybound(2)-ybound(1)) - 1 ;

Cpx = zeros(Nx+1,1) ;
Cpy = zeros(Ny+1,1) ;

Cpx(1) = 1 ;
Cpx(2) = xi ;
Cpy(1) = 1 ;
Cpy(2) = eta ;

for i = 3 : (Nx+1)
    Cpx(i) = 2*xi*Cpx(i-1) - Cpx(i-2) ;
end

for i = 3 : (Ny+1)
    Cpy(i) = 2*eta*Cpy(i-1) - Cpy(i-2) ;
end

fapprox = 0 ;
for i = 1 : (Nx+1)
    off = (i-1)*(Ny+1) ;
    for j = 1 : (Ny+1)
        fapprox = fapprox + coeff(off+j) * Cpx(i)*Cpy(j) ;
    end
end


end

