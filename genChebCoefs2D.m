function a = genChebCoefs2D( input , Nx , Ny , Mx , My , xspan , yspan , dis )

xi  = cos( (0:Mx).*(pi/Mx) ) ;
eta = cos( (0:My).*(pi/My) ) ;

% X
MxInv = zeros( Nx+1 , Nx+1 ) ;
CP = zeros( Mx+1 , Nx+1 ) ;
for i = 1 : (Mx+1)
    % For sample i, evaluate each CP term
    CP(i,1) = 1 ;
    CP(i,2) = xi(i) ;
    for j = 3 : (Nx+1)
        CP(i,j) = 2*CP(i,j-1)*xi(i) - CP(i,j-2) ;
    end
end

% Compute M inverse
for i = 1 : (Nx+1)
    MxInv(i,i) = 0.5*CP(1,i)^2 ;
    for j = 2 : (Mx)
        MxInv(i,i) = MxInv(i,i) + CP(j,i)^2 ;
    end
    MxInv(i,i) = MxInv(i,i) + 0.5*CP(Mx+1,i)^2 ; % Now this is M_ii
    MxInv(i,i) = 1/MxInv(i,i) ; % Now this is (M^-1)_ii
end
CP = CP';

% CP(1,1) = CP(1,1)*( 0.5*2^(0.5) ) ;
% CP(1,Mx+1) = CP(1,Mx+1)*( 0.5*2^(0.5) ) ;
% CP(2:(Nx+1),1) = CP(2:(Nx+1),1) .* (2^(0.5)) ;
% CP(2:(Nx+1),Mx+1) = CP(2:(Nx+1),Mx+1) .* (2^(0.5)) ;
% CP(2:(Nx+1),2:(Mx)) = CP(2:(Nx+1),2:(Mx)) .* 2 ;
% if Nx == Mx
%     CP(Nx+1,:) = CP(Nx+1,:) .* 0.5 ;
% end

Cxbar = MxInv * CP * diag( [0.5^(0.5) ones(1,Mx-1) 0.5^(0.5)] ) ;


% Y
MyInv = zeros( Ny+1 , Ny+1 ) ;
CP = zeros( My+1 , Ny+1 ) ;
for i = 1 : (My+1)
    % For sample i, evaluate each CP term
    CP(i,1) = 1 ;
    CP(i,2) = eta(i) ;
    for j = 3 : (Ny+1)
        CP(i,j) = 2*CP(i,j-1)*eta(i) - CP(i,j-2) ;
    end
end

% Compute M inverse
for i = 1 : (Ny+1)
    MyInv(i,i) = 0.5*CP(1,i)^2 ;
    for j = 2 : (My)
        MyInv(i,i) = MyInv(i,i) + CP(j,i)^2 ;
    end
    MyInv(i,i) = MyInv(i,i) + 0.5*CP(My+1,i)^2 ; % Now this is M_ii
    MyInv(i,i) = 1/MyInv(i,i) ; % Now this is (M^-1)_ii
end
CP = CP' ;
% CP(1,1) = CP(1,1)*( 0.5*2^(0.5) ) ;
% CP(1,My+1) = CP(1,My+1)*( 0.5*2^(0.5) ) ;
% CP(2:(Ny+1),1) = CP(2:(Ny+1),1) .* (2^(0.5)) ;
% CP(2:(Ny+1),My+1) = CP(2:(Ny+1),My+1) .* (2^(0.5)) ;
% CP(2:(Ny+1),2:(My)) = CP(2:(Ny+1),2:(My)) .* 2 ;
% if Ny == My
%     CP(Ny+1,:) = CP(Ny+1,:) .* 0.5 ;
% end



%% New

fbar = zeros( (Mx+1)*(My+1) , 1 ) ;

if dis == 1
  for i = 1 : (Mx+1) % This constructs fbar for the discrete case

      fbar( (i-1)*(My+1) + 1 ) = (0.5^.5)*input( My+1, (Mx+2)-i ) ;

      for j = 2 : (My)
          fbar( (i-1)*(My+1) + j ) = input( (Mx+2)-j,(Mx+2)-i  ) ;
      end

      fbar( (i-1)*(My+1) + My+1 ) = (0.5^.5)*input( 1, (Mx+2)-i  ) ;

  end
else
  for i = 1 : (Mx+1) % This constructs fbar for the fuctional case

      x = xspan(1) + (xi(i)+1)/2*(xspan(2)-xspan(1)) ;
      y = yspan(1) + (eta(1)+1)/2*(yspan(2)-yspan(1)) ;
      fbar( (i-1)*(My+1) + 1 ) = (0.5^.5)*input( x , y ) ;


      for j = 2 : (My)
          y = yspan(1) + (eta(j)+1)/2*(yspan(2)-yspan(1)) ;
          fbar( (i-1)*(My+1) + j ) = input( x , y ) ;
      end

      y = yspan(1) + (eta(My+1)+1)/2*(yspan(2)-yspan(1)) ;
      fbar( (i-1)*(My+1) + My+1 ) = (0.5^.5)*input( x , y ) ;

  end
end


fbar(1:(My+1)) = fbar(1:(My+1)).*(0.5^.5) ;
fbar((1+(My)*(Mx+1)):(My+1)*(Mx+1)) = fbar((1+(My)*(Mx+1)):(My+1)*(Mx+1)).*(0.5^.5) ;

% Get ybar
Cybar = MyInv * CP * diag( [0.5^(0.5) ones(1,My-1) 0.5^(0.5)] ) ;
ybar = zeros( (Mx+1)*(Ny+1) , 1 ) ;
for i = 1 : (Mx + 1)
    ybar( ( (i-1)*(Ny+1)+1) : ( (i)*(Ny+1)) ) = Cybar * fbar( ( (i-1)*(My+1)+1) : ( (i)*(My+1)) );
end

% Get coefficients
a = zeros( (Nx+1)*(Ny+1) , 1 ) ;
for i = 1 : (Nx + 1)
    for j = 1 : (Mx + 1)
        a( ( (i-1)*(Ny+1)+1) : ( (i)*(Ny+1)) ) = a( ( (i-1)*(Ny+1)+1) : ( (i)*(Ny+1)) ) + Cxbar(i,j)*ybar( ( (j-1)*(Ny+1)+1) : ( (j)*(Ny+1)) );
    end
end

end
