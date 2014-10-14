function [c, ceq]=hestonconstraint(x)
%     x=[ mu, alpha, beta, sigma, rho]
        
c(1) = 2*x(2)*x(3)-x(4)*x(4);  % 2*alpha*beta < sigma^2
c(2) = (-1/x(3)^3)*(exp(-x(3))+x(3)-1+4*((2+x(3))*exp(-x(3))+x(3)-2)*x(5)*x(5))*x(1)*x(4)*x(4);
ceq=[];

end