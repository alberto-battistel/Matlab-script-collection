function [Z, varargout] = RsCPERct(V,f)
%**************************************************************************
%
%**************************************************************************
Rs = 10^V(1);
Q1 = 10^V(2);
a1 = 1./(10^V(3)+1);
Rct = 10^V(4);

w = 2*pi*f;

YCPE1 = (1i*w).^a1*Q1;


Yct = (Rct).^-1;

Z = Rs + (YCPE1 + Yct).^-1;

if nargout > 1
    varargout{1,1} = zeros(length(V),1);
    varargout{1,1}(1,1) = Rs;
    varargout{1,1}(2,1) = Q1;
    varargout{1,1}(3,1) = a1;
    varargout{1,1}(4,1) = Rct;

end
return