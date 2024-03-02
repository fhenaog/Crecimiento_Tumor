function S = gsua_sample(S)
% S = gsua_sample(S)
%
% S                 Structure array with the input and output fields
% INPUT FIELDS
% S.x               Matrix with nominal values in the first row and uncertainty percent in the second one (2xNp)
% S.N               Sample size (number of sets of parameters to test)
% S.Np              Number of factors
% S.sample_method   'Uniform','LatinHypercube','Normal'
% S.rngN            Random seed (if rngN = [] or 0, then rngN = random)
% OUTPUT FIELDS
% S.M               Sample matrix of parameters with one parameter by column and one sample by row (NxNp)
%
% Global sensitivity and uncertainty analysis using GSUA Toolbox
% https://bit.ly/Matlab_GSUA
% (c) Carlos Mario Vélez S. 2022
% Universidad EAFIT, Medellin, Antioquia, Colombia
% https://sis-control.blogspot.com/


N = S.N;
x = S.x;
Np = S.Np;
sample_method = S.sample_method;

if isempty(S.rngN) || S.rngN == 0
    S.rngN = randi(1000); %  Returns a pseudorandom scalar integer between 1 and 1000.
end

% Generation of matrix M (N x Np) with all parameters for uncertanity and sensitivity analysis (one parameter by column and one sample by row)
rng(S.rngN); % Random seed
M = zeros(N,Np);
switch sample_method
    case 'Normal'
        for k=1:Np
            mu = x(1,k);
            std = x(2,k)*x(1,k)/100;
            pdfun = makedist('Normal','mu',mu,'sigma',std);
            M(:,k) = random(pdfun,1,N); % Normal distribution with two values given by the range of parameters
        end
    case 'Uniform'
        for k = 1:Np
            xmin = x(1,k) - x(2,k)*x(1,k)/100;
            xmax = x(1,k) + x(2,k)*x(1,k)/100;
            pdfun = makedist('Uniform','lower',xmin,'upper',xmax);
            M(:,k) = random(pdfun,1,N); % Uniform distribution between two values given by the range of parameters
        end        
    otherwise
        if ~strcmp(sample_method,'LatinHypercube')
            disp('Unknown sampling method. Using LatinHypercube method.')
        end
        lhs = lhsdesign(N,Np); % latin-hipercube sample in the interval [0,1]
        for k = 1:Np
            xmin = x(1,k) - x(2,k)*x(1,k)/100;
            if k==4
                xmin=x(1,k);
            end
            xmax = x(1,k) + x(2,k)*x(1,k)/100;
            M(:,k) = ones(N,1)*xmin + lhs(:,k)*(xmax - xmin); % x = xmin + random*(xmax – xmin)
        end       
end
S.M = M;
    