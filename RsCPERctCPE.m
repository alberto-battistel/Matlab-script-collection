function [Z, varargout] = RsCPERctCPE(V,f)
%**************************************************************************
%
%**************************************************************************
Rs = 10^V(1);
Q1 = 10^V(2);
a1 = 1./(10^V(3)+1);
Rct = 10^V(4);
Q2 = 10^V(5);
a2 = 1./(10^V(6)+1);

w = 2*pi*f;

YCPE1 = (1i*w).^a1*Q1;
YCPE2 = (1i*w).^a2*Q2;

Yct = (YCPE2.^-1 + Rct).^-1;

Z = Rs + (YCPE1 + Yct).^-1;

if nargout > 1
    varargout{1,1} = zeros(length(V),1);
    varargout{1,1}(1,1) = Rs;
    varargout{1,1}(2,1) = Q1;
    varargout{1,1}(3,1) = a1;
    varargout{1,1}(4,1) = Rct;
    varargout{1,1}(5,1) = Q2;
    varargout{1,1}(6,1) = a2;
end
return