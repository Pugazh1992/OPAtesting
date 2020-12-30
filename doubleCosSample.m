function [ out ] = doubleCosSample( t0, tc, tf, M )
%cosSample Simple Cosine Sample Functon
%   t0 start time
%   tf end time
%   M number of samples

    out = zeros(1,2*M+1);

    j = M:-1:0; % reverse point generation to produce output from -1 to 1
    ksi = cos((j)*pi/M);
    tmp1 = 0.5*(ksi + 1)*(tc-t0) + t0; % not discrepancy
    tmp2 = 0.5*(ksi + 1)*(tf-tc) + tc; % not discrepancy
    
    out(1:M+1) =  tmp1;
    out(M+2:end) =  tmp2(2:end);
    
end


% function [ out ] = cosSample( t0, tf, M )
% %cosSample Simmple Cosine Sample Functon
% %   t0 start time
% %   tf end time
% %   M number of samples
% 
%     j = 1 : (M+1);
%     
%     out = 1-cos((j-1).*pi./M);
%     
% 
%     out = 0.5 * (tf - t0) .* (out) + t0;
% end
% 
