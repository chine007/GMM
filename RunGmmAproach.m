% =============================================================================================
%                            GMM ESTIMATION ROUTINE OF THE SV MODEL
% Jiang G. J. and J. L. Knight, 2002, "Estimation of continuous-time processes via the 
% empirical characteristic function," Journal of Business & Economic Statistics, 20(2).
%
% =============================================================================================

clc
clear variables
close all

%% INPUT
fid=fopen('../sp500.csv'); 
data = textscan(fid, '%s %f %f %f %f %f %f', -1, 'headerlines', 1, 'Delimiter',','); 
fclose(fid); 

%%
logS=log(data{5});
Model.Data = (logS(10055:1:12581)-logS(10054:1:12580)).*100;  % Year 1990-1999
Model.Name = 'Heston';        % 'Heston'
Model.DeltaT = 1/252;      % recommended: 1/12 for monthly data, 1/252 for daily data, etc
Model.MatlabDisp = 'off';   % 'off'|'iter'|'notify'|'final'  (default: off)
Model.Disp = 'y';           % 'y'|'n' (Print results in Matlab's command window, draws graphs)
Model.Iters = 1;            % # of iterations of the weighting matrix (traditionally = 1) 
Model.q = 20;               % # of lags in the weight matrix is estimated by the Barlett kernel proposed by Newey and West (1987),
                            % Model.q = 0 reduces the spectral density matrix to the sample covariance matrix   
                            

% ESTIMATION
Results = GMMestimation(Model);

