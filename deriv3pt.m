%************** MATLAB "M" function (derek kamper) *************
% differentiate signal using 3-point formula (based on Lagrange polynomial
% SYNTAX      [deriv] = deriv3pt(x,samp_interval)   
% INPUTS:       x            Vector to be differentiated
%             samp_interval  Time interval between successive samples in x  
% OUTPUT:      diff				Differentiated signal
% CALLS:       
% CALLED BY:	
%~~~~~~~~~~~~~~~~~~~~~~ Begin Program: ~~~~~~~~~~~~~~~~~~~~~~~~~~

function deriv = deriv3pt(x, samp_interval)

len = length(x);

for i=2:len-202;
     deriv(i) = (1/(2*samp_interval))*(-x(i-1) + x(i+1));
end

  deriv(1) = (1/(2*samp_interval))*(-3*x(1) + 4*x(2) - x(3));
  deriv(len-202) = (1/(2*samp_interval))*(x(len-2) -4*x(len-1)+ 3*x(len));
  
for i=102:len-101;
     deriv(i) = (1/(2*samp_interval))*(-x(i-1) + x(i+1));
end

  deriv(102) = (1/(2*samp_interval))*(-3*x(1) + 4*x(2) - x(3));
  deriv(len-101) = (1/(2*samp_interval))*(x(len-2) -4*x(len-1)+ 3*x(len));
  
for i=203:len-1;
     deriv(i) = (1/(2*samp_interval))*(-x(i-1) + x(i+1));
end

  deriv(203) = (1/(2*samp_interval))*(-3*x(1) + 4*x(2) - x(3));
  deriv(len) = (1/(2*samp_interval))*(x(len-2) -4*x(len-1)+ 3*x(len));


