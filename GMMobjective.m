function J = GMMobjective(Params, Model, W)
% =============================================================================================
% Objective function for SV models
% 
% INPUT: 
% OUTPUT: J
%
% email: hitdragon@163.com
% =============================================================================================

Data = Model.Data;
Data2= Data.^2;
Data_2F = Data(3:end);
Data_2L = Data(1:end-2);
Data_3F = Data(4:end);
Data_3L = Data(1:end-3);
Data_4F = Data(5:end);
Data_4L = Data(1:end-4);
Data_5F = Data(6:end);
Data_5L = Data(1:end-5);
Data_6F = Data(7:end);
Data_6L = Data(1:end-6);

Nobs = length(Data);

DeltaT = Model.DeltaT;
mu    = Params(1);
alpha = Params(2);
beta  = Params(3);
sigma = Params(4); 
rho   = Params(5);

% Calculate the sample moments 
switch Model.Name
    case 'Heston' 
%         g1 = mean(Data) - mu;
        g1 = sum(Data)/Nobs - mu;
        g2 = var(Data) - alpha;
        g3 = skewness(Data, 0) - 3*(exp(-beta) + beta - 1)*alpha*rho*sigma/beta^2;
        g4 = kurtosis(Data, 0) - 3*var(Data) - (3/beta^3)*(exp(-beta)+beta-1+4*((2+beta)*exp(-beta)+beta-2)*rho*rho)*alpha*sigma*sigma;
        g5 = (1/(Nobs-2))*sum((Data_2F-mu).*(Data_2F-mu).*(Data_2L-mu).*(Data_2L-mu)) - (1/(2*beta^3))*exp(-3*beta)*(exp(beta)-1)*(exp(beta)-1+4*rho*rho*(exp(beta)-beta-1))*alpha*sigma*sigma;
        g6 = (1/(Nobs-3))*sum((Data_3F-mu).*(Data_3F-mu).*(Data_3L-mu).*(Data_3L-mu)) - (1/(2*beta^3))*exp(-4*beta)*(exp(beta)-1)*(exp(beta)-1+4*rho*rho*(exp(beta)-beta-1))*alpha*sigma*sigma;
        g7 = (1/(Nobs-4))*sum((Data_4F-mu).*(Data_4F-mu).*(Data_4L-mu).*(Data_4L-mu)) - (1/(2*beta^3))*exp(-5*beta)*(exp(beta)-1)*(exp(beta)-1+4*rho*rho*(exp(beta)-beta-1))*alpha*sigma*sigma;
        g8 = (1/(Nobs-5))*sum((Data_5F-mu).*(Data_5F-mu).*(Data_5L-mu).*(Data_5L-mu)) - (1/(2*beta^3))*exp(-6*beta)*(exp(beta)-1)*(exp(beta)-1+4*rho*rho*(exp(beta)-beta-1))*alpha*sigma*sigma;
        g9 = (1/(Nobs-6))*sum((Data_6F-mu).*(Data_6F-mu).*(Data_6L-mu).*(Data_6L-mu)) - (1/(2*beta^3))*exp(-7*beta)*(exp(beta)-1)*(exp(beta)-1+4*rho*rho*(exp(beta)-beta-1))*alpha*sigma*sigma;
end
g = [g1 g2 g3 g4 g5 g6 g7 g8 g9];                                       
J = g*W*g';

end

