% plot_pop_encoding_slope_JCSOinTT.m

% script used for plotting in pop_encoding_slope_JCSOinTT.m


% %
% % none zero threshold
ind_nonzero_JC1 = nonzero_inJC(:,1);
ind_nonzero_JC2 = nonzero_inJC(:,2);
ind_nonzero_JC3 = nonzero_inJC(:,3);
ind_nonzero_SO1 = nonzero_inSO(:,1);
ind_nonzero_SO2 = nonzero_inSO(:,2);
ind_nonzero_SO3 = nonzero_inSO(:,3);

ind_G = ismember(monkeyname_list,'G');
ind_J = ismember(monkeyname_list,'J');
ind_B = ind_G | ind_J;

% % %
% ind_nonzero_SO1 = nonzero_inSO(:,1) & nonzero_inSO(:,2);
% ind_nonzero_SO2 = nonzero_inSO(:,1) & nonzero_inSO(:,2);


% none zero threshold
% ind_nonzero_JC1 = ones(size(nonzero_inJC(:,1)));
% ind_nonzero_JC2 = ones(size(nonzero_inJC(:,2)));
% ind_nonzero_JC3 = ones(size(nonzero_inJC(:,3)));
% ind_nonzero_SO1 = ones(size(nonzero_inSO(:,1)));
% ind_nonzero_SO2 = ones(size(nonzero_inSO(:,2)));
% ind_nonzero_SO3 = ones(size(nonzero_inSO(:,3)));


% % %
% all three time windows
%
JCTWtypes = {'postoffer',  'postoffer',  'postjuice'};
SOTWtypes = {'postoffer1', 'postoffer2', 'postjuice'};
scatcols = {'b','c','g'};
nTWtypes = length(JCTWtypes);

figure;
set(gcf,'position',[110 65 1550 1550], 'PaperPositionMode','auto')
axes('position',[.02 .97 .2 .05]);
text(0,0,{[classname,', ',slopesignname]},'fontsize',11);
axis off

if doexamplecell
    ind_example = ismember(cellnames_iclass,examplecells);
end
    
for iTWtype = 1:nTWtypes
    
    JCTWtype = JCTWtypes{iTWtype};
    SOTWtype = SOTWtypes{iTWtype};
    
    eval(['ind_nonzero_JCSO = ind_nonzero_JC',num2str(iTWtype),' & ind_nonzero_SO',num2str(iTWtype),';']) 
    ind_good = ind_noneout & ind_nonzero_JCSO;
    % ind_good = ind_nonzero_JCSO;
    % ind_good = ind_noneout;
    cellnum = sum(ind_good);
              
    %
    plotypes = {'FR', 'slope', 'actrange'};
    titletypes = {'firing rate', 'tuning slope \beta', 'activity range'};
%     plotypes = {'intcept', 'slope', 'actrange'};
%     titletypes = {'tuning intercept', 'tuning slope \beta', 'activity range'};
    nplots = length(plotypes);
    for iplot = 1:nplots
        plottype = plotypes{iplot};
        titletype = titletypes{iplot};                
        subplot(nTWtypes,nplots+1,iplot+(iTWtype-1)*(nplots+1));
        if isequal(plottype,'Rsq')
            ind_inf = isinf(Rsq_inJC(:,iTWtype)) | isinf(Rsq_inSO(:,iTWtype));            
        else
            ind_inf = logical(zeros(size(ind_noneout)));
        end  
        try
            ind_consislope = slope_inJC(:,iTWtype).*slope_inSO(:,iTWtype)>0 & ind_rightslope(:,iTWtype);
        catch
            ind_consislope = slope_inJC(:,iTWtype).*slope_inSO(:,iTWtype)>0;
        end
        % ind_consislope = ones(size(slope_inJC(:,iTWtype)));
        eval(['XXX = ',plottype,'_inJC(ind_good & ~ind_inf & ind_consislope, iTWtype);'])
        eval(['YYY = ',plottype,'_inSO(ind_good & ~ind_inf & ind_consislope, iTWtype);'])        
        [bb, bbrange] = regress(YYY,XXX); % assuming 0 intercept
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
            scat_col = scatcols{iTWtype};
            sscat1 = scatter(XXX(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                             YYY(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                             80,'o','markerfacecolor',scat_col,'markeredgecolor','k');
            sscat1.MarkerFaceAlpha = .6;
        end
        end
        % title(titletype);
        xlabel([titletype, ' in Task 1']);
        ylabel([titletype, ' in Task 2']);
        text(minXY+(maxXY-minXY)/10,maxXY-(maxXY-minXY)/10,...
            {['t test: p=',num2str(p_ttest,'%1.1g')],['Wilcoxon: p=',num2str(p_wil,'%1.1g')],...
             ['Y = ', num2str(bb,'%.2f'),'X'],...
             ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
        box off
        axis([minXY maxXY minXY maxXY])
        axis square     
        set(gca,'FontSize',14)
%         if isequal(plottype,'FR')
%             set(gca, 'XTick', [0 10 20 30 40], ...                  
%                      'XTickLabel', {'0','10','20','30','40'}, ...                    
%                      'YTick', [0 10 20 30 40], ...
%                      'YTickLabel', {'0','10','20','30','40'});
%         elseif isequal(plottype,'slope')
%             set(gca, 'XTick', [0 2 4], ...                  
%                      'XTickLabel', {'0','2','4'}, ...                    
%                      'YTick', [0 2 4], ...
%                      'YTickLabel', {'0','2','4'});
%         end
    end   
    %
    subplot(nTWtypes,nplots+1,nplots+1+(iTWtype-1)*(nplots+1));
    XXX = [(steepness_inSO(ind_good & ~ind_inf & ind_consislope))-(steepness_inJC(ind_good & ~ind_inf & ind_consislope))];
%     YYY = [(slope_inSO(ind_good & ~ind_inf & ind_consislope,iTWtype))-(slope_inJC(ind_good & ~ind_inf & ind_consislope,iTWtype))]; 
%     YYY2 = [abs(slope_inSO(ind_good & ~ind_inf & ind_consislope,iTWtype))+abs(slope_inJC(ind_good & ~ind_inf & ind_consislope,iTWtype))]/2; 
%     YYY = YYY1./YYY2;
    YYY = [(actrange_inSO(ind_good & ~ind_inf & ind_consislope,iTWtype))-(actrange_inJC(ind_good & ~ind_inf & ind_consislope,iTWtype))]; 
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
    XX = [floor(min(XXX))-2,ceil(max(XXX))+2];
    YY = [floor(min(YYY)),ceil(max(YYY))];
    Yfit = aa(2)*XX+aa(1);
    [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
    [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
    hold on; 
    plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
    plot(XXX, YYY, 'ko','MarkerSize',8);
    if doexamplecell
        scat_col = scatcols{iTWtype};
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
    % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
    % title('\Delta steepness and \Delta tuning slope');
    xlabel('\Delta steepness \eta (Task 2 - Task 1)');
    % ylabel('\Delta slope \beta (Task 2 - Task 1)');
    ylabel('\Delta activity range (Task 2 - Task 1)');
    text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
        {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
         ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
         ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
    box off
    axis([XX YY])
    axis square    
    set(gca,'FontSize',14);  
    %
    axes('position',[.02 1-0.25*iTWtype .2 .05]);
    text(0,0,{[JCTWtype,' in Task 1'];[SOTWtype,' in Task 2']},'fontsize',11);
    axis off
    hold on
end % for iTWtypes   

% % % 
% % %
% OV CV: postoffer1 and postoffer2 merge
% CJ: postjuice
%
if ~isequal(classname,'CJ')
    iTWtype = [1 2];
    scat_col = {'b','c'};
    JCTWtype = 'postoffer';
    SOTWtype = 'postoffer 1&2 merged';
    ind_nonzero_JCSO = [ind_nonzero_JC1 & ind_nonzero_SO1; ind_nonzero_JC2 & ind_nonzero_SO2];
    ind_noneout_merg = [ind_noneout;ind_noneout];
    ind_G_merg = [ind_G; ind_G];
    ind_J_merg = [ind_J; ind_J];
    ind_B_merg = [ind_B; ind_B];
    ind_good = ind_noneout_merg & ind_nonzero_JCSO;
    ind_good = ind_good & ind_B_merg;
    % ind_good = ind_nonzero_JCSO;
    % ind_good = ind_noneout_merg;
    try
        ind_rightslopemerg = [ind_rightslope(:,1); ind_rightslope(:,2)];
    end
    Rsqmerg_inJC = [Rsq_inJC(:,1); Rsq_inJC(:,2)];
    Rsqmerg_inSO = [Rsq_inSO(:,1); Rsq_inSO(:,2)];
    FRmerg_inJC = [FR_inJC(:,1); FR_inJC(:,2)];
    FRmerg_inSO = [FR_inSO(:,1); FR_inSO(:,2)];
    slopemerg_inJC = [slope_inJC(:,1); slope_inJC(:,2)];
    slopemerg_inSO = [slope_inSO(:,1); slope_inSO(:,2)];
    intceptmerg_inJC = [intcept_inJC(:,1); intcept_inJC(:,2)];
    intceptmerg_inSO = [intcept_inSO(:,1); intcept_inSO(:,2)];
    steepnessmerg_inJC = [steepness_inJC; steepness_inJC];
    steepnessmerg_inSO = [steepness_inSO; steepness_inSO];   
%     rhomerg_inJC = [rho_inJC; rho_inJC];
%     rhomerg_inSO = [rho_inSO; rho_inSO];
    actrangemerg_inJC = [actrange_inJC(:,1); actrange_inJC(:,2)];
    actrangemerg_inSO = [actrange_inSO(:,1); actrange_inSO(:,2)];
    cellnamesmerg_iclass = [cellnames_iclass; cellnames_iclass];
    try
        chhystmerg_inJC = [chhyst_inJC; chhyst_inJC];
        chhystmerg_inSO = [chhyst_inSO; chhyst_inSO];
    end
    %
    ind_off1 = logical([ones(size(Rsq_inSO,1),1); zeros(size(Rsq_inSO,1),1)]);
    ind_off2 = logical([zeros(size(Rsq_inSO,1),1); ones(size(Rsq_inSO,1),1)]);
elseif isequal(classname,'CJ')
    iTWtype = 3;
    scat_col = {'g'};
    JCTWtype = 'postjuice';
    SOTWtype = 'postjuice';
    eval(['ind_nonzero_JCSO = ind_nonzero_JC',num2str(iTWtype),' & ind_nonzero_SO',num2str(iTWtype),';']) 
    ind_noneout_merg = [ind_noneout];
    ind_good = ind_noneout_merg & ind_nonzero_JCSO;
    % ind_good = ind_nonzero_JCSO;
    % ind_good = ind_noneout_merg;
    try
        ind_rightslopemerg = ind_rightslope(:,3);
    end
    Rsqmerg_inJC = Rsq_inJC(:,3);
    Rsqmerg_inSO = Rsq_inSO(:,3);
    FRmerg_inJC = FR_inJC(:,3);
    FRmerg_inSO = FR_inSO(:,3);
    slopemerg_inJC = slope_inJC(:,3);
    slopemerg_inSO = slope_inSO(:,3);
    intceptmerg_inJC = intcept_inJC(:,3);
    intceptmerg_inSO = intcept_inSO(:,3);
    steepnessmerg_inJC = steepness_inJC;
    steepnessmerg_inSO = steepness_inSO;
    try
        chhystmerg_inJC = [chhyst_inJC];
        chhystmerg_inSO = [chhyst_inSO];
    end
    rhomerg_inJC = rho_inJC;
    rhomerg_inSO = rho_inSO;
    actrangemerg_inJC = actrange_inJC(:,3);
    actrangemerg_inSO = actrange_inSO(:,3);
    cellnamesmerg_iclass = [cellnames_iclass];
    ind_off1 = logical([ones(size(Rsq_inSO,1),1)]);
    ind_off2 = logical([zeros(size(Rsq_inSO,1),1)]);
end

if doexamplecell
    ind_example = ismember(cellnamesmerg_iclass,examplecells);
end

%
figure;
set(gcf,'position',[110 65 1550 550], 'PaperPositionMode','auto')
% set(gcf,'position',[110 65 1050 550], 'PaperPositionMode','auto')
%
plotypes = {'FRmerg','actrangemerg'};
titletypes = {'Mean activity', 'Activity range'};
% plotypes = {'FRmerg', 'slopemerg', 'actrangemerg'};
% titletypes = {'firing rate', 'tuning slope \beta', 'activity range'};
% plotypes = {'intceptmerg', 'slopemerg', 'actrangemerg'};
% titletypes = {'tuning intercept', 'tuning slope \beta', 'activity range'};
nplots = length(plotypes);
for iplot = 1:nplots
    plottype = plotypes{iplot};
    titletype = titletypes{iplot};                
    subplot(1,nplots+1,iplot);
%     subplot(1,nplots,iplot);
%     subplot(2,2,iplot);
    if isequal(plottype,'Rsq')
        ind_inf = isinf(Rsqmerg_inJC) | isinf(Rsqmerg_inSO);            
    else
        ind_inf = logical(zeros(size(ind_noneout_merg)));
    end  
    try 
        ind_consislope = slopemerg_inJC.*slopemerg_inSO>0 & ind_rightslopemerg;
    catch
        ind_consislope = slopemerg_inJC.*slopemerg_inSO>0;
    end
    % ind_consislope = ones(size(slopemerg_inJC));
    eval(['XXX = ',plottype,'_inJC(ind_good & ~ind_inf & ind_consislope, :);'])
    eval(['YYY = ',plottype,'_inSO(ind_good & ~ind_inf & ind_consislope, :);'])   
    % [bb, bbrange] = regress(YYY,XXX); % assuming 0 intercept
    [bb,aa]=demingRegression(XXX,YYY,1,1);
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
    %
%     % do not separate off1 and off2
%     plot(XXX, YYY, 'ko','MarkerSize',8);
%     if doexamplecell
%     try
%         if length(iTWtype)==1
%             sscat1 = scatter(XXX(ind_example(ind_good & ~ind_inf & ind_consislope)),...
%                              YYY(ind_example(ind_good & ~ind_inf & ind_consislope)),...
%                              80,'o','markerfacecolor',scat_col,'markeredgecolor','k');
%             sscat1.MarkerFaceAlpha = .6;
%         else
%             ind_exams = find(ind_example(ind_good & ~ind_inf & ind_consislope)==1);
%             sscat1 = scatter(XXX(ind_exams(1)),YYY(ind_exams(1)),...
%                              80,'o','markerfacecolor',scat_col{1},'markeredgecolor','k');
%             sscat1.MarkerFaceAlpha = .6;
%             sscat2 = scatter(XXX(ind_exams(2)),YYY(ind_exams(2)),...
%                              80,'o','markerfacecolor',scat_col{2},'markeredgecolor','k');
%             sscat2.MarkerFaceAlpha = .6;
%         end
%     end
%     end
    %
    % separate off1 and off2
    % off1
    plot(XXX(ind_off1(ind_good & ~ind_inf & ind_consislope)), YYY(ind_off1(ind_good & ~ind_inf & ind_consislope)), 'ko','MarkerSize',8);
    % off2
    plot(XXX(ind_off2(ind_good & ~ind_inf & ind_consislope)), YYY(ind_off2(ind_good & ~ind_inf & ind_consislope)), 'kd','MarkerSize',8);
    if doexamplecell
    try
        if length(iTWtype)==1
            sscat1 = scatter(XXX(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                             YYY(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                             80,'o','markerfacecolor',scat_col,'markeredgecolor','k');
            sscat1.MarkerFaceAlpha = .6;
        else
            ind_exams = find(ind_example(ind_good & ~ind_inf & ind_consislope)==1);
            sscat1 = scatter(XXX(ind_exams(1)),YYY(ind_exams(1)),...
                             80,'o','markerfacecolor',scat_col{1},'markeredgecolor','k');
            sscat1.MarkerFaceAlpha = .6;
            sscat2 = scatter(XXX(ind_exams(2)),YYY(ind_exams(2)),...
                             80,'d','markerfacecolor',scat_col{2},'markeredgecolor','k');
            sscat2.MarkerFaceAlpha = .6;
        end
    end
    end
    % title(titletype);
    xlabel([titletype, ', Task 1 (sp/s)']);
    ylabel([titletype, ', Task 2 (sp/s)']);
    text(minXY+(maxXY-minXY)/10,maxXY-(maxXY-minXY)/10,...
        {['t test: p=',num2str(p_ttest,'%1.1g')];['Wilcoxon: p=',num2str(p_wil,'%1.1g')];...
         ['Y = ', num2str(bb,'%.2f'),'X'];...
        }, 'fontsize', 12);
%          ['N = ',num2str(anacellnum),' responses']}, 'fontsize', 12);
    box off
    axis([minXY maxXY minXY maxXY])
    axis square     
    set(gca,'FontSize',14)
end
% %
subplot(1,nplots+1,nplots+1);
% subplot(1,2,nplots+1);
XXX = [(steepnessmerg_inSO(ind_good & ~ind_inf & ind_consislope))-(steepnessmerg_inJC(ind_good & ~ind_inf & ind_consislope))];
% YYY = [(slopemerg_inSO(ind_good & ~ind_inf & ind_consislope))-(slopemerg_inJC(ind_good & ~ind_inf & ind_consislope))]; 
% YYY2 = [abs(slopemerg_inSO(ind_good & ~ind_inf & ind_consislope))+abs(slopemerg_inJC(ind_good & ~ind_inf & ind_consislope))]./2; 
% YYY = YYY1./YYY2;
YYY = [(actrangemerg_inSO(ind_good & ~ind_inf & ind_consislope))-(actrangemerg_inJC(ind_good & ~ind_inf & ind_consislope))]; 
anacellnum =  sum([ind_good & ~ind_inf & ind_consislope]); 
% aa = deming(XXX,YYY);
% % 
aa_mdl1 = fitlm(XXX,YYY);
aa1 = aa_mdl1.Coefficients.Estimate;
aa_mdl2 = fitlm(YYY,XXX);
aa2 = aa_mdl2.Coefficients.Estimate;
aa2 = [(-aa2(1)/aa2(2));1/aa2(2)];
% aa = (aa1+aa2)./2;
aa = aa1;
% % 
XX = [min(XXX)-((max(XXX)-min(XXX))*0.1), max(XXX)+((max(XXX)-min(XXX))*0.1)];
YY = [min(YYY)-((max(YYY)-min(YYY))*0.1), max(YYY)+((max(YYY)-min(YYY))*0.1)];
Yfit = aa(2)*XX+aa(1);
[RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
[RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
hold on; 
plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
% % do not separate off1 and off2
% plot(XXX, YYY, 'ko','MarkerSize',8);
% separate off1 and off2
% off1
plot(XXX(ind_off1(ind_good & ~ind_inf & ind_consislope)), YYY(ind_off1(ind_good & ~ind_inf & ind_consislope)), 'ko','MarkerSize',8);
% off2
plot(XXX(ind_off2(ind_good & ~ind_inf & ind_consislope)), YYY(ind_off2(ind_good & ~ind_inf & ind_consislope)), 'kd','MarkerSize',8);
%   
plot(XX,[0 0], 'k--','LineWidth',1);
plot([0 0],YY, 'k--','LineWidth',1);
if doexamplecell
try
    if length(iTWtype)==1
        sscat1 = scatter(XXX(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                         YYY(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                         80,'o','markerfacecolor',scat_col,'markeredgecolor','k');
        sscat1.MarkerFaceAlpha = .6;
    else
        ind_exams = find(ind_example(ind_good & ~ind_inf & ind_consislope)==1);
        sscat1 = scatter(XXX(ind_exams(1)),YYY(ind_exams(1)),...
                         80,'o','markerfacecolor',scat_col{1},'markeredgecolor','k');
        sscat1.MarkerFaceAlpha = .6;
        sscat2 = scatter(XXX(ind_exams(2)),YYY(ind_exams(2)),...
                         80,'d','markerfacecolor',scat_col{2},'markeredgecolor','k');
        sscat2.MarkerFaceAlpha = .6;
    end
end
end
Sigma_ell = cov(XXX, YYY);
mu_ell(1) = mean(XXX);
mu_ell(2) = mean(YYY);  
% hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
% title('\Delta steepness and \Delta tuning slope');
% xlabel('\Delta steepness \eta (Task 2 - Task 1)');
xlabel({['Difference in sigmoid steepness'];['(\eta Task 2 - \eta Task 1)']});
% ylabel('\Delta slope \beta (Task 2 - Task 1)');
% ylabel('\Delta activity range (Task 2 - Task 1)');
ylabel({['Difference in activity range (sp/s)'];['(\Deltar Task 2 - \Deltar Task 1)']});
text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
    {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
     ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
     }, 'fontsize', 12);
     % ['N = ',num2str(anacellnum),' responses']}, 'fontsize', 12);
box off
axis([XX YY])
axis square    
set(gca,'FontSize',14)
% %
% subplot(1,nplots+2,nplots+2);
% subplot(1,2,nplots+1);
% XXX = 2*[(rhomerg_inSO(ind_good & ~ind_inf & ind_consislope))-(rhomerg_inJC(ind_good & ~ind_inf & ind_consislope))]./...
%         [(rhomerg_inSO(ind_good & ~ind_inf & ind_consislope))+(rhomerg_inJC(ind_good & ~ind_inf & ind_consislope))];
% YYY = [(actrangemerg_inSO(ind_good & ~ind_inf & ind_consislope))-(actrangemerg_inJC(ind_good & ~ind_inf & ind_consislope))]; 
% anacellnum =  sum([ind_good & ~ind_inf & ind_consislope]); 
% % aa = deming(XXX,YYY);
% % % 
% aa_mdl1 = fitlm(XXX,YYY);
% aa1 = aa_mdl1.Coefficients.Estimate;
% aa_mdl2 = fitlm(YYY,XXX);
% aa2 = aa_mdl2.Coefficients.Estimate;
% aa2 = [(-aa2(1)/aa2(2));1/aa2(2)];
% % aa = (aa1+aa2)./2;
% aa = aa1;
% % % 
% XX = [min(XXX)-((max(XXX)-min(XXX))*0.1), max(XXX)+((max(XXX)-min(XXX))*0.1)];
% YY = [min(YYY)-((max(YYY)-min(YYY))*0.1), max(YYY)+((max(YYY)-min(YYY))*0.1)];
% Yfit = aa(2)*XX+aa(1);
% [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
% [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
% hold on; 
% plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
% % % do not separate off1 and off2
% % plot(XXX, YYY, 'ko','MarkerSize',8);
% % separate off1 and off2
% % off1
% plot(XXX(ind_off1(ind_good & ~ind_inf & ind_consislope)), YYY(ind_off1(ind_good & ~ind_inf & ind_consislope)), 'ko','MarkerSize',8);
% % off2
% plot(XXX(ind_off2(ind_good & ~ind_inf & ind_consislope)), YYY(ind_off2(ind_good & ~ind_inf & ind_consislope)), 'kd','MarkerSize',8);
% %   
% plot(XX,[0 0], 'k--','LineWidth',1);
% plot([0 0],YY, 'k--','LineWidth',1);
% if doexamplecell
% try
%     if length(iTWtype)==1
%         sscat1 = scatter(XXX(ind_example(ind_good & ~ind_inf & ind_consislope)),...
%                          YYY(ind_example(ind_good & ~ind_inf & ind_consislope)),...
%                          80,'o','markerfacecolor',scat_col,'markeredgecolor','k');
%         sscat1.MarkerFaceAlpha = .6;
%     else
%         ind_exams = find(ind_example(ind_good & ~ind_inf & ind_consislope)==1);
%         sscat1 = scatter(XXX(ind_exams(1)),YYY(ind_exams(1)),...
%                          80,'o','markerfacecolor',scat_col{1},'markeredgecolor','k');
%         sscat1.MarkerFaceAlpha = .6;
%         sscat2 = scatter(XXX(ind_exams(2)),YYY(ind_exams(2)),...
%                          80,'d','markerfacecolor',scat_col{2},'markeredgecolor','k');
%         sscat2.MarkerFaceAlpha = .6;
%     end
% end
% end
% Sigma_ell = cov(XXX, YYY);
% mu_ell(1) = mean(XXX);
% mu_ell(2) = mean(YYY);  
% xlabel({['preference bias index']});
% ylabel({['Difference in activity range (sp/s)'];['(\Deltar Task 2 - \Deltar Task 1)']});
% text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
%     {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
%      ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
%      }, 'fontsize', 12);
%      % ['N = ',num2str(anacellnum),' responses']}, 'fontsize', 12);
% box off
% axis([XX YY])
% axis square    
% set(gca,'FontSize',14)
% % %
axes('position',[.06 .45 .2 .05]);
H = text(0,0,{[classname,', ',slopesignname]; [JCTWtype,' in Task 1']; [SOTWtype,' in Task 2']},'fontsize',11);
axis off
set(H,'Rotation',90);


try
% slope difference and ch hyst difference 
figure

XXX = [(chhystmerg_inSO(ind_good & ~ind_inf & ind_consislope))-(chhystmerg_inJC(ind_good & ~ind_inf & ind_consislope))];
% YYY = [(slopemerg_inSO(ind_good & ~ind_inf & ind_consislope))-(slopemerg_inJC(ind_good & ~ind_inf & ind_consislope))]; 
% YYY2 = [abs(slopemerg_inSO(ind_good & ~ind_inf & ind_consislope))+abs(slopemerg_inJC(ind_good & ~ind_inf & ind_consislope))]./2; 
% YYY = YYY1./YYY2;
YYY = [(actrangemerg_inSO(ind_good & ~ind_inf & ind_consislope))-(actrangemerg_inJC(ind_good & ~ind_inf & ind_consislope))]; 
anacellnum =  sum([ind_good & ~ind_inf & ind_consislope]); 
% aa = deming(XXX,YYY);
% % 
aa_mdl1 = fitlm(XXX,YYY);
aa1 = aa_mdl1.Coefficients.Estimate;
aa_mdl2 = fitlm(YYY,XXX);
aa2 = aa_mdl2.Coefficients.Estimate;
aa2 = [(-aa2(1)/aa2(2));1/aa2(2)];
% aa = (aa1+aa2)./2;
aa = aa1;
% % 
XX = [min(XXX)-((max(XXX)-min(XXX))*0.1), max(XXX)+((max(XXX)-min(XXX))*0.1)];
YY = [min(YYY)-((max(YYY)-min(YYY))*0.1), max(YYY)+((max(YYY)-min(YYY))*0.1)];
Yfit = aa(2)*XX+aa(1);
[RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
[RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
hold on; 
plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
% % do not separate off1 and off2
% plot(XXX, YYY, 'ko','MarkerSize',8);
% separate off1 and off2
% off1
plot(XXX(ind_off1(ind_good & ~ind_inf & ind_consislope)), YYY(ind_off1(ind_good & ~ind_inf & ind_consislope)), 'ko','MarkerSize',8);
% off2
plot(XXX(ind_off2(ind_good & ~ind_inf & ind_consislope)), YYY(ind_off2(ind_good & ~ind_inf & ind_consislope)), 'kd','MarkerSize',8);
%   
plot(XX,[0 0], 'k--','LineWidth',1);
plot([0 0],YY, 'k--','LineWidth',1);
if doexamplecell
try
    if length(iTWtype)==1
        sscat1 = scatter(XXX(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                         YYY(ind_example(ind_good & ~ind_inf & ind_consislope)),...
                         80,'o','markerfacecolor',scat_col,'markeredgecolor','k');
        sscat1.MarkerFaceAlpha = .6;
    else
        ind_exams = find(ind_example(ind_good & ~ind_inf & ind_consislope)==1);
        sscat1 = scatter(XXX(ind_exams(1)),YYY(ind_exams(1)),...
                         80,'o','markerfacecolor',scat_col{1},'markeredgecolor','k');
        sscat1.MarkerFaceAlpha = .6;
        sscat2 = scatter(XXX(ind_exams(2)),YYY(ind_exams(2)),...
                         80,'d','markerfacecolor',scat_col{2},'markeredgecolor','k');
        sscat2.MarkerFaceAlpha = .6;
    end
end
end
Sigma_ell = cov(XXX, YYY);
mu_ell(1) = mean(XXX);
mu_ell(2) = mean(YYY);  
% hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
% title('\Delta steepness and \Delta tuning slope');
xlabel({['Difference in choice hysteresis'];['(\xi Task 2 - \xi Task 1)']});
% ylabel('\Delta slope \beta (Task 2 - Task 1)');
% ylabel('\Delta activity range (Task 2 - Task 1)');
ylabel({['Difference in activity range (sp/s)'];['(\Deltar Task 2 - \Deltar Task 1)']});
text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
    {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
     ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
     }, 'fontsize', 12);
     % ['N = ',num2str(anacellnum),' responses']}, 'fontsize', 12);
box off
axis([XX YY])
axis square    
set(gca,'FontSize',14)
    
end

% % % 
% % %
% OV CV: <SO postoffer1 - JC postoffer> VS <SO postoffer2 - JC postoffer>
%
if ~isequal(classname,'CJ')
    %
    ind_nonzero_JCSO = [ind_nonzero_JC1 & ind_nonzero_SO1 & ind_nonzero_SO2];
    ind_noneout_merg = ind_noneout;
    ind_good = ind_noneout & ind_nonzero_JCSO;
    ind_good = ind_good & ind_B;
    % ind_good = ind_nonzero_JCSO;
    % ind_good = ind_noneout_merg;
    try
        ind_rightslopeJCSO = [ind_rightslope(:,1); ind_rightslope(:,2)];
    end
    Rsqmerg_inJC = [Rsq_inJC(:,1), Rsq_inJC(:,2)];
    Rsqmerg_inSO = [Rsq_inSO(:,1), Rsq_inSO(:,2)];
    FRmerg_inJC = [FR_inJC(:,1), FR_inJC(:,2)];
    FRmerg_inSO = [FR_inSO(:,1), FR_inSO(:,2)];
    slopemerg_inJC = [slope_inJC(:,1), slope_inJC(:,2)];
    slopemerg_inSO = [slope_inSO(:,1), slope_inSO(:,2)];
    intceptmerg_inJC = [intcept_inJC(:,1), intcept_inJC(:,2)];
    intceptmerg_inSO = [intcept_inSO(:,1), intcept_inSO(:,2)];
    steepnessmerg_inJC = [steepness_inJC, steepness_inJC];
    steepnessmerg_inSO = [steepness_inSO, steepness_inSO];   
%     rhomerg_inJC = [rho_inJC, rho_inJC];
%     rhomerg_inSO = [rho_inSO, rho_inSO];
    actrangemerg_inJC = [actrange_inJC(:,1), actrange_inJC(:,2)];
    actrangemerg_inSO = [actrange_inSO(:,1), actrange_inSO(:,2)];
    cellnamesmerg_iclass = [cellnames_iclass, cellnames_iclass];
    %
    ind_off1 = logical([ones(size(Rsq_inSO,1),1), zeros(size(Rsq_inSO,1),1)]);
    ind_off2 = logical([zeros(size(Rsq_inSO,1),1), ones(size(Rsq_inSO,1),1)]);
    
    %
    figure;
    set(gcf,'position',[110 65 450 450], 'PaperPositionMode','auto')
    %
%     plotypes = {'FRmerg','actrangemerg'};
%     titletypes = {'Mean activity', 'Activity range'};
    plotypes = {'actrangemerg'};
    titletypes = {'Activity range'};
    nplots = length(plotypes);
    for iplot = 1:nplots
        plottype = plotypes{iplot};
        titletype = titletypes{iplot};                
        subplot(1,nplots,iplot);
        if isequal(plottype,'Rsq')
            ind_inf = logical(sum(isinf(Rsqmerg_inJC) | isinf(Rsqmerg_inSO),2));            
        else
            ind_inf = logical(zeros(size(ind_noneout_merg)));
        end  
        try 
            ind_consislope = slopemerg_inJC.*slopemerg_inSO>0 & ind_rightslopemerg;
        catch
            ind_consislope = slopemerg_inJC.*slopemerg_inSO>0;
        end
        % ind_consislope = ones(size(slopemerg_inJC));
        ind_consislope = sum(ind_consislope,2)>1;
        eval(['dat_JC = ',plottype,'_inJC(ind_good & ~ind_inf & ind_consislope, :);'])
        eval(['dat_SO = ',plottype,'_inSO(ind_good & ~ind_inf & ind_consislope, :);'])
        XXX = dat_SO(:,1) - dat_JC(:,1);
        YYY = dat_SO(:,2) - dat_JC(:,2);        
        % [bb, bbrange] = regress(YYY,XXX); % assuming 0 intercept
        [bb,aa]=demingRegression(XXX,YYY,1,1);
        anacellnum =  sum([ind_good & ~ind_inf & ind_consislope]); 
        maxXY = ceil(max([XXX; YYY]));
        minXY = floor(min([XXX; YYY]));
        if isequal(plottype,'FR')
            maxXY = maxXY + 5;
            minXY = minXY - 5;
        end
        hold on
        plot([minXY maxXY], [minXY maxXY],'--','Color', [0 0 0],'LineWidth',0.5);
        plot([minXY maxXY], [0 0],'--','Color', [0 0 0],'LineWidth',0.5);
        plot([0 0], [minXY maxXY],'--','Color', [0 0 0],'LineWidth',0.5);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);       
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
        [~, p_ttest] = ttest(XXX, YYY);
        p_wil = signrank(XXX, YYY);
        plot(XXX, YYY, 'ko','MarkerSize',9);       
        % title(titletype);
        xlabel({['Difference in ',titletype], ['<Task 2 postoffer1 - Task 1 postoffer>']});
        ylabel({['Difference in ',titletype,] ['<Task 2 postoffer2 - Task 1 postoffer>']});
        text(minXY+(maxXY-minXY)/10,maxXY-(maxXY-minXY)/10,...
            {['t test: p=',num2str(p_ttest,'%1.1g')];['Wilcoxon: p=',num2str(p_wil,'%1.1g')];...
             ['Y = ', num2str(bb,'%.2f'),'X'];...
             ['\Delta(Y-X) = ',num2str(nanmean(YYY-XXX),'%.2f'),''];...
            }, 'fontsize', 12);
    %          ['N = ',num2str(anacellnum),' responses']}, 'fontsize', 12);
        box off
        axis([minXY maxXY minXY maxXY])
        axis square     
        set(gca,'FontSize',14)
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
    
