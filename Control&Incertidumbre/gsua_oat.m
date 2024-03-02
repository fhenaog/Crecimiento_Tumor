function S = gsua_oat(S)
% S = gsua_oat(S)
%
% S                         Structure array with input and output fields
% S.model                   Simulation model
% S.x                       Matrix with nominal values in the first row and uncertainty percent in the second one (2xNp)
% S.M                       Sample matrix of parameters with one parameter by column and one sample by row (NxNp)
% S.ynom                    Reference vectorial (1xNt) or scalar output: experimental or nominal (ynom = [] for nominal)
% S.scalar_characteistic    Scalar characteristic for computing scalar sensitivity indices
% S.S1                      First-order global sensitivity indices with output = y (NpxNt vectorial or Npx1 scalar output)
% S.S1s                     Scalar first-order global sensitivity indices with y = scalar_fun (Npx1)
% S.ST                      Empty matrix (for OAT method there are only first-order indices)
% S.STs                     Empty matrix (for OAT method there are only first-order indices)
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

M = S.M;
ynom = S.ynom; % It is used intrinsically in "eval(S.scalar_characteristic)"
N = S.N;
Np = S.Np;
Nt = size(S.t,2);
Y = zeros(N*Np,Nt);
Ys = zeros(N*Np,1);
V = zeros(Np,Nt);
Vs = zeros(Np,1);
if strcmp(S.model_type,'Simulink model') == 1
    load_system(S.model)
    model_workspace = get_param(S.model, 'ModelWorkspace');
end

tic
for i=1:Np
    if i==1
        tic
    end
    Yi = zeros(N,Nt); % Vectorial output
    Yis = zeros(N,1); % Scalar output
    Mi = repmat(S.x(1,:),N,1);
    Mi(:,i) = M(:,i);
    for j = 1:N
        if i==1 && j==10
            [h,m,s]=hms(seconds(toc*N*Np/10));
            [h2,m2,s2] = hms(datetime('now') + seconds(toc*N*Np/10));
            disp(['Estimated computing time: ' num2str(h) 'h:' num2str(m) 'm:' num2str(round(s)) 's'])
            disp(['Estimated stop time: ' num2str(h2,2) ':' num2str(m2,2) ':' num2str(s2,2)])
        end

        switch S.model_type
            case 'Simulink model'
                assignin(model_workspace,'x', Mi(j,:));
                simout = sim(S.model,'ReturnWorkspaceOutputs','on');
                y = simout.yout.signals.values';
            case 'M-file model'
                [y,~] = eval([S.model '(Mi(j,:))']);
        end
        Yi(j,:) = y;
        Yis(j) = eval(S.scalar_characteristic);
    end
    % Compute of local sensitivity indices
    Y((i-1)*N+1:i*N,:) = Yi;
    Ys((i-1)*N+1:i*N) = Yis;
    V(i,:) = var(Yi,1);
    Vs(i) = var(Yis);
end

% Normalization of sensitivity indices
S.S1 = V./sum(V,1);
S.S1s = Vs/sum(Vs);
S.ST = [];
S.STs = [];
S.Y = Y;
S.Ys = Ys;
S.Nsim = N*Np;
[h,m,s] = hms(seconds(toc)); S.tsim = [h m s];

end