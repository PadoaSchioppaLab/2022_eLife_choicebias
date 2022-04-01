% plot_pop_preference_bias_JCSOinTT.m

% script used for plotting in pop_preference_bias_JCSOinTT.m

%
% none zero threshold
ind_nonzero_JC1 = nonzero_inJC(:,1);
ind_nonzero_SO1 = nonzero_inSO(:,1);
ind_nonzero_SO2 = nonzero_inSO(:,2);
% ind_nonzero_JC1 = logical(ones(size(nonzero_inJC(:,1))));
% ind_nonzero_SO1 = logical(ones(size(nonzero_inSO(:,1))));
% ind_nonzero_SO2 = logical(ones(size(nonzero_inSO(:,2))));

% % %
% for OVA, OVB, OV cells
% plot slope/activity range in postoffer1 vs postoffer2
%
if ismember(classname,{'OVA';'OVB';'OV'})
    %
    figure;
    set(gcf,'position',[110 65 950 550], 'PaperPositionMode','auto')
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname]},'fontsize',11);
    axis off
    % %
    if doexamplecell
        ind_example = ismember(cellnames_iclass,examplecells);
    end
    % %
    
    ind_nonzero_JC1SO1 = ind_nonzero_JC1 & ind_nonzero_SO1;
    % ind_good1 = ind_nonzero_JC1SO1;
    ind_good1 = ind_nonzero_JC1SO1 & ind_noneout;
    ind_nonzero_JC1SO2 = ind_nonzero_JC1 & ind_nonzero_SO2;
    % ind_good2 = ind_nonzero_JC1SO2;
    ind_good2 = ind_nonzero_JC1SO2 & ind_noneout;
    ind_good = [ind_good1; ind_good2];
    try
        ind_good = ind_good  & [ind_rightslope(:,1); ind_rightslope(:,2)];
    end
    %
    ind_consistslope_JC1SO1 = slope_inJC(:,1).*slope_inSO(:,1) >0;
    ind_consistslope_JC1SO2 = slope_inJC(:,1).*slope_inSO(:,2) >0;
    ind_consistslope = [ind_consistslope_JC1SO1; ind_consistslope_JC1SO2];
    
    %
    ind_off1 = logical([ones(size(slope_inSO,1),1); zeros(size(slope_inSO,1),1)]);
    ind_off2 = logical([zeros(size(slope_inSO,1),1); ones(size(slope_inSO,1),1)]);
    
    %
    rhotype  = '_inJC';
    normalized_deltarho = (rho_inSO - rho_inJC)./(rho_inSO + rho_inJC);
    deltarho = (rho_inSO - rho_inJC);
    xlabelname = 'Task 1';
    plotypes = {'intcept', 'slope'};
    regressors = {'b0','b1'};
    nplots = length(plotypes);
    for iplot = 1:nplots
        plottype = plotypes{iplot};
        regsor = regressors{iplot};
        subplot(1,nplots,iplot);
        %
        eval(['XXX = ([rho',rhotype,'; rho',rhotype,']);'])
        % XXX = [normalized_deltarho; normalized_deltarho];
        % XXX = [deltarho; deltarho];
        eval(['YYY = [',plottype,'_inSO(:,1)-',plottype,'_inJC(:,1); ',plottype,'_inSO(:,2)-',plottype,'_inJC(:,1)];'])
        XXX = XXX(ind_good & ind_consistslope);
        YYY = YYY(ind_good & ind_consistslope);
        anacellnum =  sum([ind_good & ind_consistslope]); 
        %
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);       
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
        %
        % aa = deming(XXX,YYY);
        % % 
        aa_mdl1 = fitlm(XXX,YYY);
        aa1 = aa_mdl1.Coefficients.Estimate;
        aa_mdl2 = fitlm(YYY,XXX);
        aa2 = aa_mdl2.Coefficients.Estimate;
        aa2 = [(-aa2(1)/aa2(2));1/aa2(2)];
        % aa = (aa1+aa2)./2;
        aa = aa1;
        % 
        % XX = [(min(XXX))-((max(XXX)-min(XXX))/10),(max(XXX))+((max(XXX)-min(XXX))/10)];
        % XX = [0.25 4.5];
        XX = [0.25 4.5];
        if iplot == 1
            YY = [-9 16];
        elseif iplot == 2
            YY = [-2 2];
        else
            YY = [floor(min(YYY)),ceil(max(YYY))];
        end
        Yfit = aa(2)*XX+aa(1);
        %
        hold on; 
        plot(XX,[0 0],'--','LineWidth',0.5, 'Color', [0.2 0.2 0.2]);
        plot([1 1],YY,'--','LineWidth',0.5, 'Color', [0.2 0.2 0.2]);
        plot(XX,Yfit,'-','LineWidth',1.5, 'Color', [0.4 0.4 0.4]);
        % do not separate OFF1 and OFF2
        % plot(XXX, YYY, 'ko','MarkerSize',8);
        % separate OFF1 and OFF2
        plot(XXX(ind_off1(ind_good & ind_consistslope)), YYY(ind_off1(ind_good & ind_consistslope)), 'ko','MarkerSize',8);
        plot(XXX(ind_off2(ind_good & ind_consistslope)), YYY(ind_off2(ind_good & ind_consistslope)), 'kd','MarkerSize',8);
        if doexamplecell
            scat_col = scatcols;
            try
                sscat1 = scatter(XXX(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                                 YYY(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                                 80,'o','markerfacecolor',scat_col,'markeredgecolor','k');
                sscat1.MarkerFaceAlpha = .6;
            end
        end
        %
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        %
        xlabel(['\rho ',xlabelname]);
        % xlabel(['normalized \Delta\rho']);
        ylabel({['Difference in tuning ',plottype];['(',regsor,',Task2 - ',regsor,',Task1)']});
        text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
            {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
             }, 'fontsize', 12);
             % ['N = ',num2str(anacellnum),' responses']}, 'fontsize', 12);
        box off
        axis([XX YY])
        axis square    
        set(gca,'FontSize',14);   
    end   
    
    
    %
    if 0
    figure;
    set(gcf,'position',[110 65 1050 550], 'PaperPositionMode','auto')
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname]},'fontsize',11);
    axis off
    % %
    if doexamplecell
        ind_example = ismember(cellnames_iclass,examplecells);
    end
    % %
    
    ind_nonzero_JC1SO1 = ind_nonzero_JC1 & ind_nonzero_SO1;
    ind_good1 = ind_nonzero_JC1SO1;
    ind_nonzero_JC1SO2 = ind_nonzero_JC1 & ind_nonzero_SO2;
    ind_good2 = ind_nonzero_JC1SO2;
    ind_good = [ind_good1; ind_good2];
    try
        ind_good = ind_good & [ind_noneout; ind_noneout] & [ind_rightslope(:,1); ind_rightslope(:,2)];
    end
    %
    ind_consistslope_JC1SO1 = slope_QAB_inJC(:,1).*slope_QAB_inSO(:,1)>0;
    ind_consistslope_JC1SO2 = slope_QAB_inJC(:,1).*slope_QAB_inSO(:,2)>0;
    ind_consistslope = [ind_consistslope_JC1SO1; ind_consistslope_JC1SO2];
    
    %
    plotypes = {'intcept', 'slope'};
    nplots = length(plotypes);
    for iplot = 1:nplots
        plottype = plotypes{iplot};
        subplot(1,nplots,iplot);
        %
        XXX = [steepness_inSO-steepness_inJC; steepness_inSO-steepness_inJC];
        eval(['YYY = [',plottype,'_inSO(:,1)-',plottype,'_inJC(:,1); ',plottype,'_inSO(:,2)-',plottype,'_inJC(:,1)];'])
        XXX = XXX(ind_good & ind_consistslope);
        YYY = YYY(ind_good & ind_consistslope);
        anacellnum =  sum([ind_good & ind_consistslope]); 
        %
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);       
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
        %
        % aa = deming(XXX,YYY);
        % % 
        aa_mdl1 = fitlm(XXX,YYY);
        aa1 = aa_mdl1.Coefficients.Estimate;
        aa_mdl2 = fitlm(YYY,XXX);
        aa2 = aa_mdl2.Coefficients.Estimate;
        aa2 = [(-aa2(1)/aa2(2));1/aa2(2)];
        aa = (aa1+aa2)./2;
        % 
        XX = [floor(min(XXX)),ceil(max(XXX))];
        YY = [floor(min(YYY)),ceil(max(YYY))];
        Yfit = aa(2)*XX+aa(1);
        %
        hold on; 
        plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
        plot(XXX, YYY, 'ko','MarkerSize',8);
        if doexamplecell
            scat_col = scatcols;
            try
                sscat1 = scatter(XXX(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                                 YYY(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                                 80,'o','markerfacecolor',scat_col,'markeredgecolor','k');
                sscat1.MarkerFaceAlpha = .6;
            end
        end
        %
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        %
        xlabel(['\Delta \eta']);
        ylabel(['\Delta ',plottype]);
        text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
            {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(anacellnum),' responses']}, 'fontsize', 12);
        box off
        axis([XX YY])
        axis square    
        set(gca,'FontSize',14);   
    end   
    end
end




% % % 
function plotErrorEllipse(mu_ell, Sigma_ell, p_ell)
    s = -2 * log(1 - p_ell);
    [V, D] = eig(Sigma_ell * s);
    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    plot(a(1, :) + mu_ell(1), a(2, :) + mu_ell(2), 'LineWidth',1.5, 'Color',[0.4 0.4 0.4] );
end
    
