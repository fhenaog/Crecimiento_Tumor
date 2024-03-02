function [] = gsua_plot_scatter(S)
% function [] = gsua_plot_scatter(S)
%
% S                 Structure array with input and output fields
% INPUT FIELDS
% S.Ys              Scalar outputs for y = scalar_fun (Nx1 or N*Npx1)
% S.M               Sample matrix of parameters with one parameter by column and one sample by row (NxNp)
% S.sens_method     Sensitivity method: 'OAT', 'Saltelli', 'Bruteforce'
% S.factor_names    Factor names in a cell
% S.N               Sample size (number of sets of parameters to test)
% S.Np              Number of factors
%
% Global sensitivity and uncertainty analysis using GSUA Toolbox
% https://bit.ly/Matlab_GSUA
% (c) Carlos Mario Vélez S. 2022
% Universidad EAFIT, Medellin, Antioquia, Colombia
% https://sis-control.blogspot.com/

Ys = S.Ys;
M = S.M;
sens_method = S.sens_method;
factor_names = S.factor_names;
N = S.N;
Np = S.Np;
nf = floor(sqrt(Np));
nc = nf+ceil((Np-nf^2)/nf);
sgtitle('Scatter plots')
for i=1:Np
    subplot(nf,nc,i)
    if strcmp(sens_method,'OAT')
        plot(M(:,i),Ys((i-1)*N+1:i*N),'.')
    else
        plot(M(:,i),Ys,'.')
    end
    xlabel([factor_names{i}])
    ylabel('Ys')
end
end