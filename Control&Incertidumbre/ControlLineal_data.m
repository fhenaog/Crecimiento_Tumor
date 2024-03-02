% Global sensitivity and uncertainty analysis using GSUA Toolbox
% https://bit.ly/Matlab_GSUA
% (c) Carlos Mario Vélez S. 2022
% Universidad EAFIT, Medellin, Antioquia, Colombia
% https://sis-control.blogspot.com/

% Model description:                            S.description('text...');
% Factor names:                                 S.factor_names = {'text1','text2',..};
% Nominal factors and respective uncertainties: S.x = [v1 v2,...,vn; u1 u2 ... un];
% Simulation model:                             S.model = 'text...';

S.description = 'Control del Crecimiento de un Tumor';
S.factor_names = {'a1','a2','a3','b1','alpha','c1','c2','c3','c4','d1','d2','r1','r2','s','rho','N0','T0','I0','u0'}; 
S.x = [0.2 0.3 0.1 1.0 0.3 1.0 0.5 1.0 1.0 0.2 1.0 1.5 1.0 0.33 0.01 0.4934 0.4970 0.4649 0.10
        20 20 20 20 20 20 20 20 20 20 20 19 20 20 20 20 20 20 20]; 
S.model = 'ControlLinealGSUA';
%S.model = 'mass_spring_ode';
%S.model = 'mass_spring_ss';