function S = gsua_methods(S)
% S = gsua_methods(S)
%
% S                         Structure array with input and output fields
% FIELDS
% S.sens_method             Sensitivity method: 'OAT', 'Saltelli', 'Bruteforce'
% S.model                   Simulation model
% S.x                       Matrix with nominal values in the first row and uncertainty percent in the second one (2xNp)
% S.M                       Sample matrix of parameters with one parameter by column and one sample by row (NxNp)
% S.ynom                    Reference vectorial (1xNt) or scalar output: experimental or nominal (ynom = [] for nominal)
% S.scalar_characteristic   Scalar characteristic for computing scalar sensitivity indices
% S.S1                      First-order global sensitivity indices with output = y (NpxNt vectorial or Npx1 scalar output)
% S.S1s                     Scalar first-order global sensitivity indices with y = scalar_fun (Npx1)
% S.ST                      Total global sensitivity indices with output = y (NpxNt vectorial or Npx1 scalar)
% S.STs                     Scalar total global sensitivity indices with y = scalar_fun (Npx1)
% S.Y                       Outputs for every factor sample: vectorial time response (NxNt or N*NpxNt) or scalar (N*Npx1)
% S.Ys                      Scalar outputs for y = scalar_fun (Nx1 or N*Npx1)
% S.ynom                    Nominal time response (if input S.ynom = [])
% S.t                       Time simulation vector with fixed step (1xNt)
% S.Nsim                    Number of simulations
% S.tsim                    Simulation time in datatime format
%
% Global sensitivity and uncertainty analysis using GSUA Toolbox
% https://bit.ly/Matlab_GSUA
% (c) Carlos Mario Vélez S. 2022
% Universidad EAFIT, Medellin, Antioquia, Colombia
% https://sis-control.blogspot.com/


% Compute of t vector and, optionally, nominal time response
if  strcmp(S.model_type, 'Simulink model') == 1
    load_system(S.model)
    model_workspace = get_param(S.model, 'ModelWorkspace');
    assignin(model_workspace,'x', S.x(1,:));
    simout = sim(S.model,'ReturnWorkspaceOutputs','on');
    S.t = simout.yout.time';
    if strcmp(S.ynom,'Nominal')==1
        S.ynom = simout.yout.signals.values';
    end
else
    [ynom,S.t] = eval([S.model '(S.x(1,:))']);
    if strcmp(S.ynom,'Nominal')==1
        S.ynom = ynom;
    end
end
if length(S.t)>2
    S.Ts = S.t(2)-S.t(1); S.tmax = S.t(end);
end

switch S.sens_method
    case 'OAT'
        S = gsua_oat(S);
    case 'Saltelli'
        S = gsua_saltelli(S);
    case 'Bruteforce'
        S = gsua_bruteforce(S);
end