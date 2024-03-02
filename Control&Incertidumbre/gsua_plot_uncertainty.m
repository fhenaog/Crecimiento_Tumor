function [] = gsua_plot_uncertainty(S)
% function [] = gsua_plot_uncertainty(S)
%
% S    Structure array with all sensitivity index fields
%
% Global sensitivity and uncertainty analysis using GSUA Toolbox
% https://bit.ly/Matlab_GSUA
% (c) Carlos Mario Vélez S. 2022
% Universidad EAFIT, Medellin, Antioquia, Colombia
% https://sis-control.blogspot.com/

Y = S.Y;
t = S.t;
ynom = S.ynom;

if length(t)>1
    plot(t,Y,'b-',t,ynom,'r-')
    xlabel('t (sec)')
    ylabel('y(t)')
    title({'Uncertainty plot (y = time response)';''})
else
    histogram(Y,'Normalization','probability')
    xlabel('y (escalar output)')
    ylabel('Frequency')
    xline(ynom, 'Color', 'r', 'LineWidth', 2);
    title({'Uncertainty analysis (y = scalar)';''})
end
end