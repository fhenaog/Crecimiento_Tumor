function [] = gsua_plot_sensitivity(index_type,S,plot_type)
% function [] = gsua_plot_s(type,S,plot_type)
%
% type          Sensitivity index type: 'S1', 'S1s', 'ST', 'STs'
% S             Structure array with all sensitivity index fields
% plot_type     Plot type: 'pie', 'bar'
%
% Global sensitivity and uncertainty analysis using GSUA Toolbox
% https://bit.ly/Matlab_GSUA
% (c) Carlos Mario Vélez S. 2022
% Universidad EAFIT, Medellin, Antioquia, Colombia
% https://sis-control.blogspot.com/

if nargin==2 || isempty(plot_type)
    plot_type = 'barh';
end
t = S.t;
factor_names = S.factor_names;
fspec = '%.1f';
Np = S.Np;
Nt = length(t);

if strcmp(index_type,'S1') && Nt>1
    sens_case = 'S1 (y = time response)';
    Stype = S.S1;
elseif strcmp(index_type,'S1') && Nt==1
    sens_case = 'S1 (y = scalar function)';
    Stype = S.S1;
elseif strcmp(index_type,'ST') && Nt>1
    sens_case = 'ST (y = time response)';
    Stype = S.ST;
elseif strcmp(index_type,'ST') && Nt==1
    sens_case = 'ST (y = scalar function)';
    Stype = S.ST;
elseif strcmp(index_type,'S1s')
    sens_case = 'S1s (y = scalar characteristic)';
    Stype = S.S1s;
elseif strcmp(index_type,'STs')
    sens_case = 'STs (y = scalar characteristic)';
    Stype = S.STs;
else
    disp('Sensitivity index type is incorrect')
    return
end

switch sens_case
    case 'S1 (y = time response)'
        area(t,Stype')
        xlabel('t (sec)')
        title({'S1 (y = time response)';''})
        legend(factor_names)
    case 'S1 (y = scalar function)'
        switch plot_type
            case 'pie'
                inter = 1-sum(Stype);
                labels_pie = cell(1,Np+1);
                for i=1:Np
                    labels_pie(i) = {[factor_names{i} ' (' num2str(Stype(i)*100,fspec) '%)']};
                end
                labels_pie(Np+1) = {['Inter'  ' (' num2str(inter*100,fspec) '%)']};
                pie([Stype;inter],labels_pie)
                title({'S1 (y = scalar function)';''})
            case 'barh'
                inter = 1-sum(Stype);
                labels_bar = cell(1,Np+1);
                for i=1:Np
                    labels_bar(i) = {[num2str(Stype(i)*100,fspec) '%']};
                end
                labels_bar(Np+1) = {[num2str(inter*100,fspec) '%']};
                X = reordercats(categorical([factor_names 'Inter']),[factor_names 'Inter']);
                b = barh(X,[Stype;inter]*100,'BaseValue',0);
                xmin = min([Stype*100;0]); xmax = max(Stype*100)+20; xlim([xmin xmax])
                title({'S1 (y = scalar function)';''})
                xtips1 = b(1).YEndPoints + 1;
                ytips1 = b(1).XEndPoints;
                text(xtips1,ytips1,labels_bar,'FontSize',8)
        end
    case 'ST (y = time response)'
        area(t,Stype')
        xlabel('t (sec)')
        title({'ST (y = time response)';''})
        legend(factor_names)
    case 'ST (y = scalar function)'
        switch plot_type
            case 'pie'
                labels_pie = cell(1,Np);
                for i=1:Np
                    labels_pie(i) = {[factor_names{i} ' (' num2str(Stype(i)*100,fspec) '%)']};
                end
                pie(Stype,labels_pie(1:Np))
                title({'ST (y = scalar function)';''})
            case 'barh'
                labels_bar = cell(1,Np);
                for i=1:Np
                    labels_bar(i) = {[num2str(Stype(i)*100,fspec) '%']};
                end                
                X = reordercats(categorical(factor_names),factor_names);
                b = barh(X,Stype*100,'BaseValue',0);
                xmin = min([Stype*100;0]); xmax = max(Stype*100)+20; xlim([xmin xmax])
                title({'ST (y = scalar function)';''})
                xtips1 = b(1).YEndPoints + 0.4;
                ytips1 = b(1).XEndPoints;
                text(xtips1,ytips1,labels_bar(1:Np),'FontSize',8)
        end
    case 'S1s (y = scalar characteristic)'
        switch plot_type
            case 'pie'
                inter = 1 -sum(Stype);
                labels_pie = cell(1,Np+1);
                for i=1:Np
                    labels_pie(i) = {[factor_names{i} ' (' num2str(Stype(i)*100,fspec) '%)']};
                end
                labels_pie(Np+1) = {['Inter'  ' (' num2str(inter*100,fspec) '%)']};
                pie([Stype;inter],labels_pie)
                title({['S1s (y = ' S.scalar_characteristic ')'];''})
            case 'barh'
                inter = 1-sum(Stype);
                labels_bar = cell(1,Np+1);
                for i=1:Np
                    labels_bar(i) = {[num2str(Stype(i)*100,fspec) '%']};
                end
                labels_bar(Np+1) = {[num2str(inter*100,fspec) '%']};
                X = reordercats(categorical([factor_names 'Inter']),[factor_names 'Inter']);
                b = barh(X,[Stype;inter]*100,'BaseValue',0);
                xmin = min([Stype*100;0]); xmax = max([Stype*100;xmin])+20; xlim([xmin xmax])
                title({['S1s (y = ' S.scalar_characteristic ')'];''})
                xtips1 = max(2,b(1).YEndPoints + 2);
                ytips1 = b(1).XEndPoints;
                text(xtips1,ytips1,labels_bar,'FontSize',8)
        end
    case 'STs (y = scalar characteristic)'
        switch plot_type
            case 'pie'
                labels_pie = cell(1,Np);
                for i=1:Np
                    labels_pie(i) = {[factor_names{i} ' (' num2str(Stype(i)*100,fspec) '%)']};
                end                
                pie(Stype,labels_pie(1:Np))
                title({['STs (y = ' S.scalar_characteristic ')'];''})
            case 'barh'
                labels_bar = cell(1,Np);
                for i=1:Np
                    labels_bar(i) = {[num2str(Stype(i)*100,fspec) '%']};
                end                
                X = reordercats(categorical(factor_names),factor_names);
                b = barh(X,Stype*100,'BaseValue',0);
                xmin = min([Stype*100;0]); xmax = max(Stype*100)+20; xlim([xmin xmax])
                title({['STs (y = ' S.scalar_characteristic ')'];''})
                xtips1 = b(1).YEndPoints + 1;
                ytips1 = b(1).XEndPoints;
                text(xtips1,ytips1,labels_bar(1:Np),'FontSize',8)
        end
end
