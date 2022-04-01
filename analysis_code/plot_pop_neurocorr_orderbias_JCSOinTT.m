% plot_pop_neurocorr_orderbias_JCSOinTT.m

% script used for plotting in pop_neurocorr_orderbias_JCSOinTT.m

%
% none zero threshold
ind_nonzero_JC1 = nonzero_inJC(:,1);
ind_nonzero_JC2 = nonzero_inJC(:,2);
ind_nonzero_JC3 = nonzero_inJC(:,3);
ind_nonzero_SO1 = nonzero_inSO(:,1);
ind_nonzero_SO2 = nonzero_inSO(:,2);
ind_nonzero_SO3 = nonzero_inSO(:,3);

% ind_nonzero_JC1 = logical(ones(size(nonzero_inJC(:,1))));
% ind_nonzero_JC2 = logical(ones(size(nonzero_inJC(:,2))));
% ind_nonzero_JC3 = logical(ones(size(nonzero_inJC(:,3))));
% ind_nonzero_SO1 = logical(ones(size(nonzero_inSO(:,1))));
% ind_nonzero_SO2 = logical(ones(size(nonzero_inSO(:,2))));
% ind_nonzero_SO3 = logical(ones(size(nonzero_inSO(:,3))));

% % %
% for OVA, OVB, OV cells
% plot slope/activity range in postoffer1 vs postoffer2
%
if ismember(classname,{'OVA';'OVB';'OV';'CV'})
    %
    figure;
    set(gcf,'position',[110 65 1550 550], 'PaperPositionMode','auto')
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname]},'fontsize',11);
    axis off
    % %
    if doexamplecell
        ind_example = ismember(cellnames_iclass,examplecells);
    end
    % %
    
    ind_nonzero_SO12 = ind_nonzero_SO1 & ind_nonzero_SO2;
    % ind_good = ind_goodorder & ind_nonzero_SO12;
    ind_good = ind_nonzero_SO12;
    cellnum = sum(ind_good);
              
    %
    plotypes = {'intcept', 'slope'};
    titletypes = {'tuning intercept', 'tuning slope \beta'};
    nplots = length(plotypes);
    for iplot = 1:nplots
        plottype = plotypes{iplot};
        titletype = titletypes{iplot};                
        subplot(1,nplots+1,iplot);
        if isequal(plottype,'Rsq')
            ind_inf = isinf(Rsq_inSO(:,1)) | isinf(Rsq_inSO(:,2));            
        else
            ind_inf = logical(zeros(size(ind_goodorder)));
        end  
        try
            ind_consislope = slope_inSO(:,1).*slope_inSO(:,2)>0 & ind_rightslope;
        catch
            ind_consislope = slope_inSO(:,1)>0 & slope_inSO(:,2)>0;
        end
        % nd_consislope = ones(size(slope_inSO(:,1)));
        eval(['XXX = ',plottype,'_inSO(ind_good & ~ind_inf & ind_consislope, 1);'])
        eval(['YYY = ',plottype,'_inSO(ind_good & ~ind_inf & ind_consislope, 2);'])        
        anacellnum =  sum([ind_good & ~ind_inf & ind_consislope]); 
        maxXY = ceil(max([XXX; YYY]));
        minXY = floor(min([XXX; YYY]));
        if isequal(plottype,'FR')
            maxXY = maxXY + 5;
            minXY = minXY - 5;
        end
        hold on
        plot([minXY maxXY], [minXY maxXY],'--','Color', [0.5 0.5 0.5],'LineWidth',1);   
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);       
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
        [~, p_ttest] = ttest(XXX, YYY);
        p_wil = signrank(XXX, YYY);
        plot(XXX, YYY, 'ko','MarkerSize',8);
        if doexamplecell
        try
            scat_col = scatcols;
            sscat1 = scatter(XXX(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                             YYY(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                             80,'o','markerfacecolor',scat_col,'markeredgecolor','k');
            sscat1.MarkerFaceAlpha = .6;
        end
        end
        % title(titletype);
        xlabel([titletype, ' in postoffer1']);
        ylabel([titletype, ' in postoffer2']);
        text(minXY+(maxXY-minXY)/10,maxXY-(maxXY-minXY)/10,...
            {['t test: p=',num2str(p_ttest,'%1.1g')],['Wilcoxon: p=',num2str(p_wil,'%1.1g')]...
             ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
        box off
        axis([minXY maxXY minXY maxXY])
        axis square     
        set(gca,'FontSize',14)
    end   
    %
    subplot(1,nplots+1,nplots+1);
    XXX = [(orderbias_inSO(ind_good & ~ind_inf & ind_consislope))];
    YYY = [(slope_inSO(ind_good & ~ind_inf & ind_consislope,2))-(slope_inSO(ind_good & ~ind_inf & ind_consislope,1))]; 
    anacellnum =  sum([ind_good & ~ind_inf & ind_consislope]); 
    % aa = deming(XXX,YYY);
    % % 
    aa_mdl1 = fitlm(XXX,YYY);
    aa1 = aa_mdl1.Coefficients.Estimate;
    aa_mdl2 = fitlm(YYY,XXX);
    aa2 = aa_mdl2.Coefficients.Estimate;
    aa2 = [(-aa2(1)/aa2(2));1/aa2(2)];
    aa = (aa1+aa2)./2;
    % % 
    XX = [floor(min(XXX)),ceil(max(XXX))];
    YY = [floor(min(YYY)),ceil(max(YYY))];
    Yfit = aa(2)*XX+aa(1);
    [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
    [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
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
    Sigma_ell = cov(XXX, YYY);
    mu_ell(1) = mean(XXX);
    mu_ell(2) = mean(YYY);  
    xlabel('order bias \epsilon (Task 2)');
    ylabel('\Delta slope (offer2-offer1)');
    text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
        {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
         ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
         ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
    box off
    axis([XX YY])
    axis square    
    set(gca,'FontSize',14);    
end


% % %
% for CV cells
% plot slope/activity range in postjuice time window (chosen value 1 v.s. chosen value 2)
%
if ismember(classname,{'CV'})
    
    %
    % ind_nonzero_SOpj = nonzeros_chV12_inSO(:,1) & nonzeros_chV12_inSO(:,2);
    ind_nonzero_SOpj =  ind_nonzero_SO3;
    ind_good = ind_nonzero_SOpj;
    % ind_good = logical(ones(size(ind_nonzero_SOpj)));
    cellnum = sum(ind_good);
              
    %
    figure;
    set(gcf,'position',[110 65 1550 550], 'PaperPositionMode','auto')
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname],['postjuice time window']},'fontsize',11);
    axis off

    %
    plotypes = {'intcepts', 'slopes'};
    titletypes = {'tuning intercept', 'tuning slope \beta'};
    nplots = length(plotypes);
    for iplot = 1:nplots
        plottype = plotypes{iplot};
        titletype = titletypes{iplot};                
        subplot(1,nplots+1,iplot);
        try
            ind_consislope = slopes_chV12_inSO(:,1).*slopes_chV12_inSO(:,2)>0 & ind_rightslope;
        catch
            ind_consislope = slopes_chV12_inSO(:,1)>0 & slopes_chV12_inSO(:,2)>0;
        end
        % ind_consislope = ones(size(slopes_chV12_inSO(:,1)));
        eval(['XXX = ',plottype,'_chV12_inSO(ind_good & ind_consislope, 1);'])
        eval(['YYY = ',plottype,'_chV12_inSO(ind_good & ind_consislope, 2);'])        
        anacellnum =  sum([ind_good & ind_consislope]); 
        maxXY = ceil(max([XXX; YYY]));
        minXY = floor(min([XXX; YYY]));
        if isequal(plottype,'FR')
            maxXY = maxXY + 5;
            minXY = minXY - 5;
        end
        hold on
        plot([minXY maxXY], [minXY maxXY],'--','Color', [0.5 0.5 0.5],'LineWidth',1);   
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);       
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
        [~, p_ttest] = ttest(XXX, YYY);
        p_wil = signrank(XXX, YYY);
        plot(XXX, YYY, 'ko','MarkerSize',8);
        if doexamplecell
        try
            scat_col = scatcols;
            sscat1 = scatter(XXX(ind_example(ind_good & ind_consislope)),...
                             YYY(ind_example(ind_good & ind_consislope)),...
                             80,'o','markerfacecolor',scat_col,'markeredgecolor','k');
            sscat1.MarkerFaceAlpha = .6;
        end
        end
        % title(titletype);
        xlabel([titletype, ' chosen value 1']);
        ylabel([titletype, ' chosen value 2']);
        text(minXY+(maxXY-minXY)/10,maxXY-(maxXY-minXY)/10,...
            {['t test: p=',num2str(p_ttest,'%1.1g')],['Wilcoxon: p=',num2str(p_wil,'%1.1g')]...
             ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
        box off
        axis([minXY maxXY minXY maxXY])
        axis square     
        set(gca,'FontSize',14)
    end   
    %
    subplot(1,nplots+1,nplots+1);
    XXX = [(orderbias_inSO(ind_good & ind_consislope))];
    YYY = [(slopes_chV12_inSO(ind_good & ind_consislope,2))-(slopes_chV12_inSO(ind_good & ind_consislope,1))]; 
    anacellnum =  sum([ind_good & ind_consislope]); 
    % aa = deming(XXX,YYY);
    % % 
    aa_mdl1 = fitlm(XXX,YYY);
    aa1 = aa_mdl1.Coefficients.Estimate;
    aa_mdl2 = fitlm(YYY,XXX);
    aa2 = aa_mdl2.Coefficients.Estimate;
    aa2 = [(-aa2(1)/aa2(2));1/aa2(2)];
    aa = (aa1+aa2)./2;
    % % 
    XX = [floor(min(XXX)),ceil(max(XXX))];
    YY = [floor(min(YYY)),ceil(max(YYY))];
    Yfit = aa(2)*XX+aa(1);
    [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
    [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
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
    Sigma_ell = cov(XXX, YYY);
    mu_ell(1) = mean(XXX);
    mu_ell(2) = mean(YYY);  
    xlabel('order bias \epsilon (Task 2)');
    ylabel('\Delta slope (chV2 - chV1)');
    text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
        {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
         ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
         ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
    box off
    axis([XX YY])
    axis square    
    set(gca,'FontSize',14);    
    
end



% % % 
function plotErrorEllipse(mu_ell, Sigma_ell, p_ell)
    s = -2 * log(1 - p_ell);
    [V, D] = eig(Sigma_ell * s);
    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    plot(a(1, :) + mu_ell(1), a(2, :) + mu_ell(2), 'LineWidth',1.5, 'Color',[0.4 0.4 0.4] );
end
    
