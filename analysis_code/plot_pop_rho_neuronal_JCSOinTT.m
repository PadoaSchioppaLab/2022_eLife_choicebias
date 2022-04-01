% plot_pop_rho_neuronal_JCSOinTT.m

% script used for plotting in pop_rho_neuronal_JCSOinTT.m

% % %
% for CV cells
% plot rho_neuronal related analysis
%
dononzeroslopes = 1; % change ind_good: none zero slopes of both A and B
doconsislopes   = 1; % change ind_consislope: consistent slopes 

if ismember(classname,{'CV','CV_JCSOoverlap','CV_SOonly','allneurons','taskrelatedneurons'})
    
    % % % 
    % rho neuronal v.s. rho behavioral
    % for Task 1 (JC)
    if 1
    figure;
    set(gcf,'position',[110 65 1550 550], 'PaperPositionMode','auto')
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname];['Task 1']},'fontsize',11);
    axis off
    %
    % three time windows: JC: postoffer latedelay postjuice               
    twins = {'postoffer', 'latedelay', 'postjuice'};
    ntwins = length(twins);
    for iTW = 1:ntwins
        twin = twins{iTW};
        %       
        subplot(1,ntwins,iTW);
        eval(['ind_good = nonzeros_QAB_inJC.',twin,'(:,1).*nonzeros_QAB_inJC.',twin,'(:,2) == 1;']) 
        % eval(['ind_good = nonzeros_QAB_inJC.',twin,'(:,1) + nonzeros_QAB_inJC.',twin,'(:,2) > 0;']) 
        if ~dononzeroslopes, ind_good = logical(ones(size(ind_good))); end
        ind_good = ind_good & ind_noneout;
        eval(['ind_cosislope = slopes_QAB_inJC.',twin,'(:,1).*slopes_QAB_inJC.',twin,'(:,2) > 0;'])
        if ~doconsislopes, ind_cosislope = logical(ones(size(ind_cosislope))); end
        eval(['YYY = (slopes_QAB_inJC.',twin,'(ind_good & ind_cosislope,1)./slopes_QAB_inJC.',twin,'(ind_good & ind_cosislope,2));']) 
        XXX = (rhos_inJC(ind_good & ind_cosislope));   
%         %
        % remove outlier
        kout = 1.5;
        quantiles_XXX = quantile(XXX,[0.25 0.5 0.75]);
        IQR_XXX = quantiles_XXX(3) - quantiles_XXX(1);
        outlier_XXX = [quantiles_XXX(1)-kout*IQR_XXX, quantiles_XXX(3)+kout*IQR_XXX];
        ind_goodXXX = (XXX > outlier_XXX(1)  & XXX < outlier_XXX(2));
        quantiles_YYY = quantile(YYY,[0.25 0.5 0.75]);
        IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
        outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
        ind_goodYYY = (YYY > outlier_YYY(1)  & YYY < outlier_YYY(2));
        ind_goodXY = ind_goodXXX & ind_goodYYY;
        %
        XXX = XXX(ind_goodXY);
        YYY = YYY(ind_goodXY);
        %
        anacellnum =  sum(length(XXX)); 
        %
        maxXY = ceil(max([XXX; YYY]));
        minXY = floor(min([XXX; YYY]));           
        XX = [minXY,maxXY];
        YY = [minXY,maxXY];
        mdl = fitlm(XXX,YYY);
        aa = mdl.Coefficients.Estimate;
        paa1 = coefTest(mdl,[1,0],0);
        paa2 = coefTest(mdl,[0,1],1);
        Yfit = aa(2)*XX+aa(1);
        %
        [~, p_ttest] = ttest(XXX, YYY);
        p_wil = signrank(XXX, YYY);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        %
        hold on
        plot(XX, YY, '--','Color', [0.8 0.8 0.8],'LineWidth',1);
        plot(XX, Yfit, '-','Color', [0.4 0.4 0.4],'LineWidth',3);
        plot(XXX, YYY, 'ko','MarkerSize',8);
        ylabel([' \rho(neuronal) in ',twin]);
        xlabel([' \rho(behavioral) of Task 1']);
        text(minXY+(maxXY-minXY)/20,maxXY-(maxXY-minXY)/9,...
            {['t test: p=',num2str(p_ttest,'%1.1g')];['Wilcoxon: p=',num2str(p_wil,'%1.1g')];...
             ['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 9);
        text(minXY+(maxXY-minXY)/2,minXY+(maxXY-minXY)/9,...
            {['y = ',num2str(aa(1),'%.2f'),'+',num2str(aa(2),'%.2f'),'x'];...
             ['slope: p (null=1) =',num2str(paa2,'%1.1g')]; ...
             ['intercept: p (null=0) =',num2str(paa1,'%1.1g')]; ...
            }, 'fontsize', 9); 
        box off
        axis([minXY maxXY minXY maxXY])
        axis square     
        set(gca,'FontSize',12)
    end   
        % 
    % plot rho neuronal vs rho task1 and task2 together
    h =  findobj('type','figure');
    %
    figure(length(h)+1);
    set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname]; ['Task 1']},'fontsize',11);
    axis off
    %
    twins = {'postoffer', 'latedelay', 'postjuice'};
    ntwins = length(twins);
    for iTW = 1:ntwins
        twin = twins{iTW};
        %       
        subplot(1,ntwins,iTW);
        eval(['ind_good = nonzeros_QAB_inJC.',twin,'(:,1).*nonzeros_QAB_inJC.',twin,'(:,2) == 1;']) 
        % eval(['ind_good = nonzeros_QAB_inJC.',twin,'(:,1) + nonzeros_QAB_inJC.',twin,'(:,2) > 0;']) 
        if ~dononzeroslopes, ind_good = logical(ones(size(ind_good))); end
        ind_good = ind_good & ind_noneout;
        eval(['ind_cosislope = slopes_QAB_inJC.',twin,'(:,1).*slopes_QAB_inJC.',twin,'(:,2) > 0;'])
        if ~doconsislopes, ind_cosislope = logical(ones(size(ind_cosislope))); end
        eval(['YYY = (slopes_QAB_inJC.',twin,'(ind_good & ind_cosislope,1)./slopes_QAB_inJC.',twin,'(ind_good & ind_cosislope,2));']) 
        XXX1 = (rhos_inJC(ind_good & ind_cosislope)); 
        XXX2 = (rhos_inSO(ind_good & ind_cosislope));  
%         %
        % remove outlier
        kout = 1.5;
        quantiles_XXX1 = quantile(XXX1,[0.25 0.5 0.75]);
        IQR_XXX1 = quantiles_XXX1(3) - quantiles_XXX1(1);
        outlier_XXX1 = [quantiles_XXX1(1)-kout*IQR_XXX1, quantiles_XXX1(3)+kout*IQR_XXX1];
        ind_goodXXX1 = (XXX1 > outlier_XXX1(1)  & XXX1 < outlier_XXX1(2));
        quantiles_XXX2 = quantile(XXX2,[0.25 0.5 0.75]);
        IQR_XXX2 = quantiles_XXX2(3) - quantiles_XXX2(1);
        outlier_XXX2 = [quantiles_XXX2(1)-kout*IQR_XXX2, quantiles_XXX2(3)+kout*IQR_XXX2];
        ind_goodXXX2 = (XXX2 > outlier_XXX2(1)  & XXX2 < outlier_XXX2(2));
        quantiles_YYY = quantile(YYY,[0.25 0.5 0.75]);
        IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
        outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
        ind_goodYYY = (YYY > outlier_YYY(1)  & YYY < outlier_YYY(2));
        ind_goodXY = ind_goodXXX1 & ind_goodXXX2 & ind_goodYYY;
%         ind_goodXY = logical(ones(size(XXX1)));
        %
        XXX1 = XXX1(ind_goodXY);
        XXX2 = XXX2(ind_goodXY);
        YYY = YYY(ind_goodXY);
        anacellnum =  sum(length(XXX1)); 
        %
        maxXY = ceil(max([XXX1; XXX2; YYY]));
        minXY = floor(min([XXX1; XXX2; YYY]));
        %
        mdl1 = fitlm(XXX1,YYY);
        aa1 = mdl1.Coefficients.Estimate;
        paa1_1 = coefTest(mdl1,[1,0],0);
        paa2_1 = coefTest(mdl1,[0,1],1);
        mdl2 = fitlm(XXX2,YYY);
        aa2 = mdl2.Coefficients.Estimate;
        paa1_2 = coefTest(mdl2,[1,0],0);
        paa2_2 = coefTest(mdl2,[0,1],1);
        XX = [minXY maxXY];
        Yfit1 = aa1(2)*XX+aa1(1);
        Yfit2 = aa2(2)*XX+aa2(1);
        %
        Sigma_ell1 = cov(XXX1, YYY);
        mu_ell1(1) = mean(XXX1);
        mu_ell1(2) = mean(YYY);
        Sigma_ell2 = cov(XXX2, YYY);
        mu_ell2(1) = mean(XXX2);
        mu_ell2(2) = mean(YYY);
        %
        [RR_Spe1,pp_Spe1] = corr(XXX1,YYY,'Type','Spearman');
        [RR_Pea1,pp_Pea1] = corr(XXX1,YYY,'Type','Pearson');
        [RR_Spe2,pp_Spe2] = corr(XXX2,YYY,'Type','Spearman');
        [RR_Pea2,pp_Pea2] = corr(XXX2,YYY,'Type','Pearson');
        %
        [~, p_tt1] = ttest(XXX1,YYY);
        [p_wil1,~] = signrank(XXX1,YYY);
        [~, p_tt2] = ttest(XXX2,YYY);
        [p_wil2,~] = signrank(XXX2,YYY);
        %
        % ancova
        YYY_all = [YYY; YYY]; % rho_neuronal
        XXX_all = [XXX1; XXX2]; % rho_JC; rho_SO
        taskgroup = [repelem({'Task1'},length(XXX1)),repelem({'Task2'},length(XXX2))]';
        [h,atab,ctab,stats] = aoctool(XXX_all,YYY_all,taskgroup,0.05,'','','','off','separate lines');
        p_slp = atab{4,6};
        %
        hold on        
        plot(XXX1, YYY, 'bo','MarkerSize',8);
        plot(XXX2, YYY, 'ro','MarkerSize',8);
        plot(XX,Yfit1,'b-','LineWidth',3);
        plot(XX,Yfit2,'r-','LineWidth',3);
        plotErrorEllipse(mu_ell1, Sigma_ell1, 0.90, 'b')
        plotErrorEllipse(mu_ell2, Sigma_ell2, 0.90, 'r')
        plot([minXY maxXY], [minXY maxXY],'--','Color', [0.5 0.5 0.5],'LineWidth',1);   
        ylabel([' \rho(neuronal) in ',twin]);
        xlabel(['\rho(behavioral)']);
        legend({'Task 1','Task 2'});
        text(minXY+(maxXY-minXY)/2,minXY+(maxXY-minXY)/7,...
            {['Task 1: y = ',num2str(aa1(1),'%.2f'),'+',num2str(aa1(2),'%.2f'),'x'];...
             ['Task 1 slope: p (null=1) =',num2str(paa2_1,'%1.1g')]; ...
             ['Task 1 intercept: p (null=0) =',num2str(paa1_1,'%1.1g')]; ...
             ['Task 1 Spearman: r=',num2str(RR_Spe1,'%1.1g'), ', p=',num2str(pp_Spe1,'%1.1g')]; ...
             ['Task 1 Pearson: r=',num2str(RR_Pea1,'%1.1g'),', p=',num2str(pp_Pea1,'%1.1g')];...
             ['Task 1 t test: p=',num2str(p_tt1,'%1.1g')];...
             ['Task 1 Wilcoxon: p=',num2str(p_wil1,'%1.1g')];...
             ['Task 2: y = ',num2str(aa2(1),'%.2f'),'+',num2str(aa2(2),'%.2f'),'x'];...
             ['Task 2 slope: p (null=1) =',num2str(paa2_2,'%1.1g')]; ...
             ['Task 2 intercept: p (null=0) =',num2str(paa1_2,'%1.1g')]; ...
             ['Task 2 Spearman: r=',num2str(RR_Spe2,'%1.1g'), ', p=',num2str(pp_Spe2,'%1.1g')]; ...
             ['Task 2 Pearson: r=',num2str(RR_Pea2,'%1.1g'),', p=',num2str(pp_Pea2,'%1.1g')];...
             ['Task 2 t test: p=',num2str(p_tt2,'%1.1g')];...
             ['Task 2 Wilcoxon: p=',num2str(p_wil2,'%1.1g')];...
             ['slope difference: p (ANCOVA) =',num2str(p_slp,'%1.1g')]; ...
             ['N = ', num2str(length(YYY)),'cells']; ...
            }, 'fontsize', 9); 
        box off
        axis([minXY maxXY minXY maxXY])
        axis square     
        set(gca,'FontSize',14)
    end % for itype   
    end
    
    
    % % % 
    % rho neuronal v.s. rho behavioral
    % for Task 2 (SO)
    % compare with rho task 1 and rho task 2
    rhos_bhv = {'inJC','inSO'};
    bhvnames = {'Task 1', 'Task 2'};
    for ibhv = 1:length(rhos_bhv)
        rho_bhv = rhos_bhv{ibhv};
        bhvname = bhvnames{ibhv};
        %
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        axes('position',[.02 .97 .2 .05]);
        text(0,0,{[classname,', ',slopesignname]; ['Task 2']},'fontsize',11);
        axis off
        %
        % six types of rho neuronal               
%         plottypes = {'inSOoff1', 'inSOoff2', 'inAB12', 'inBA12', 'inABpj', 'inBApj'};
%         labeltypes = {'post-offer1', 'post-offer2','AB trials', 'BA trials', 'AB (post-juice)', 'BA (post-juice)'};
        plottypes = {'inSOoff1', 'inSOoff2', 'inAB12', 'inBA12'};
        labeltypes = {'post-offer1', 'post-offer2','AB trials', 'BA trials'};
        ntypes = length(plottypes);
        for itype = 1:ntypes       
            plottype = plottypes{itype};
            labeltype = labeltypes{itype};
            %       
            subplot(2,ceil(ntypes/2),itype);
            eval(['ind_good = nonzeros_QAB_inSO.',plottype,'(:,1).*nonzeros_QAB_inSO.',plottype,'(:,2) == 1;']) 
            % eval(['ind_good = nonzeros_QAB_inSO.',plottype,'(:,1) + nonzeros_QAB_inSO.',plottype,'(:,2) > 0;'])
            if ~dononzeroslopes, ind_good = logical(ones(size(ind_good))); end
            ind_good = ind_good & ind_noneout;
            eval(['ind_cosislope = slopes_QAB_inSO.',plottype,'(:,1).*slopes_QAB_inSO.',plottype,'(:,2) > 0;'])
            if ~doconsislopes, ind_cosislope = logical(ones(size(ind_cosislope))); end
            eval(['YYY = (slopes_QAB_inSO.',plottype,'(ind_good & ind_cosislope,1)./slopes_QAB_inSO.',plottype,'(ind_good & ind_cosislope,2));']) 
            eval(['XXX = (rhos_',rho_bhv,'(ind_good & ind_cosislope));'])            
%             %
            % remove outlier
            kout = 1.5;
            quantiles_XXX = quantile(XXX,[0.25 0.5 0.75]);
            IQR_XXX = quantiles_XXX(3) - quantiles_XXX(1);
            outlier_XXX = [quantiles_XXX(1)-kout*IQR_XXX, quantiles_XXX(3)+kout*IQR_XXX];
            ind_goodXXX = (XXX > outlier_XXX(1)  & XXX < outlier_XXX(2));
            quantiles_YYY = quantile(YYY,[0.25 0.5 0.75]);
            IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
            outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
            ind_goodYYY = (YYY > outlier_YYY(1)  & YYY < outlier_YYY(2));
            ind_goodXY = ind_goodXXX & ind_goodYYY;
%             ind_goodXY = logical(ones(size(XXX1)));
            %
            XXX = XXX(ind_goodXY);
            YYY = YYY(ind_goodXY);
%             %
            anacellnum =  sum(length(XXX)); 
            %
            maxXY = ceil(max([XXX; YYY]));
            minXY = floor(min([XXX; YYY]));           
            XX = [minXY,maxXY];
            YY = [minXY,maxXY];
            mdl = fitlm(XXX,YYY);
            aa = mdl.Coefficients.Estimate;
            paa1 = coefTest(mdl,[1,0],0);
            paa2 = coefTest(mdl,[0,1],1);
            Yfit = aa(2)*XX+aa(1);
            %
            [~, p_ttest] = ttest(XXX, YYY);
            p_wil = signrank(XXX, YYY);
            [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
            [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
            %
            hold on
            plot(XX, YY, '--','Color', [0.8 0.8 0.8],'LineWidth',1);
            plot(XX, Yfit, '-','Color', [0.4 0.4 0.4],'LineWidth',3);
            plot(XXX, YYY, 'ko','MarkerSize',8);
            ylabel(['\rho(neuronal) ',labeltype]);
            xlabel(['\rho(behavioral) of ',bhvname]);
            text(minXY+(maxXY-minXY)/20,maxXY-(maxXY-minXY)/9,...
                {['t test: p=',num2str(p_ttest,'%1.1g')];['Wilcoxon: p=',num2str(p_wil,'%1.1g')];...
                 ['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
                 ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
                 ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 9);
            text(minXY+(maxXY-minXY)/2,minXY+(maxXY-minXY)/9,...
                {['y = ',num2str(aa(1),'%.2f'),'+',num2str(aa(2),'%.2f'),'x'];...
                 ['slope: p (null=1) =',num2str(paa2,'%1.1g')]; ...
                 ['intercept: p (null=0) =',num2str(paa1,'%1.1g')]; ...
                }, 'fontsize', 9); 
            box off
            axis([minXY maxXY minXY maxXY])
            axis square     
            set(gca,'FontSize',12)
        end % for itype
    end % for ibhv
    
    % 
    % plot rho neuronal vs rho task1 and task2 together
    h =  findobj('type','figure');
    %
    figure(length(h)+1);
    set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname]; ['Task 2']},'fontsize',11);
    axis off
    %
    % six types of rho neuronal               
    plottypes = {'inSOoff1', 'inSOoff2', 'inAB12', 'inBA12', 'inABpj', 'inBApj'};
    labeltypes = {'post-offer1', 'post-offer2','AB trials', 'BA trials', 'AB (post-juice)', 'BA (post-juice)'};
    ntypes = length(plottypes);
    for itype = 1:ntypes       
        plottype = plottypes{itype};
        labeltype = labeltypes{itype};
        %       
        subplot(2,ceil(ntypes/2),itype);
        eval(['ind_good = nonzeros_QAB_inSO.',plottype,'(:,1).*nonzeros_QAB_inSO.',plottype,'(:,2) == 1;'])
        % eval(['ind_good = nonzeros_QAB_inSO.',plottype,'(:,1) + nonzeros_QAB_inSO.',plottype,'(:,2) > 0;'])
        if ~dononzeroslopes, ind_good = logical(ones(size(ind_good))); end
        ind_good = ind_good & ind_noneout;
        eval(['ind_cosislope = slopes_QAB_inSO.',plottype,'(:,1).*slopes_QAB_inSO.',plottype,'(:,2) > 0;'])
        if ~doconsislopes,  ind_cosislope = logical(ones(size(ind_cosislope))); end
        eval(['YYY = ((slopes_QAB_inSO.',plottype,'(ind_good & ind_cosislope,1)./slopes_QAB_inSO.',plottype,'(ind_good & ind_cosislope,2)));']) 
        XXX1 = (rhos_inJC(ind_good & ind_cosislope)); 
        XXX2 = (rhos_inSO(ind_good & ind_cosislope)); 
        % 
        % remove outlier
        kout = 1.5;
        quantiles_XXX1 = quantile(XXX1,[0.25 0.5 0.75]);
        IQR_XXX1 = quantiles_XXX1(3) - quantiles_XXX1(1);
        outlier_XXX1 = [quantiles_XXX1(1)-kout*IQR_XXX1, quantiles_XXX1(3)+kout*IQR_XXX1];
        ind_goodXXX1 = (XXX1 > outlier_XXX1(1)  & XXX1 < outlier_XXX1(2));
        quantiles_XXX2 = quantile(XXX2,[0.25 0.5 0.75]);
        IQR_XXX2 = quantiles_XXX2(3) - quantiles_XXX2(1);
        outlier_XXX2 = [quantiles_XXX2(1)-kout*IQR_XXX2, quantiles_XXX2(3)+kout*IQR_XXX2];
        ind_goodXXX2 = (XXX2 > outlier_XXX2(1)  & XXX2 < outlier_XXX2(2));
        quantiles_YYY = quantile(YYY,[0.25 0.5 0.75]);
        IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
        outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
        ind_goodYYY = (YYY > outlier_YYY(1)  & YYY < outlier_YYY(2));
        ind_goodXY = ind_goodXXX1 & ind_goodXXX2 & ind_goodYYY;
%         ind_goodXY = logical(ones(size(XXX1)));
        %
        XXX1 = XXX1(ind_goodXY);
        XXX2 = XXX2(ind_goodXY);
        YYY = YYY(ind_goodXY);
        %
        anacellnum =  sum(length(XXX1)); 
        %
        maxXY = ceil(max([XXX1; XXX2; YYY]));
        minXY = floor(min([XXX1; XXX2; YYY]));
        %
        mdl1 = fitlm(XXX1,YYY);
        aa1 = mdl1.Coefficients.Estimate;
        paa1_1 = coefTest(mdl1,[1,0],0);
        paa2_1 = coefTest(mdl1,[0,1],1);
        mdl2 = fitlm(XXX2,YYY);
        aa2 = mdl2.Coefficients.Estimate;
        paa1_2 = coefTest(mdl2,[1,0],0);
        paa2_2 = coefTest(mdl2,[0,1],1);
        XX = [minXY maxXY];
        Yfit1 = aa1(2)*XX+aa1(1);
        Yfit2 = aa2(2)*XX+aa2(1);
        %
        Sigma_ell1 = cov(XXX1, YYY);
        mu_ell1(1) = mean(XXX1);
        mu_ell1(2) = mean(YYY);
        Sigma_ell2 = cov(XXX2, YYY);
        mu_ell2(1) = mean(XXX2);
        mu_ell2(2) = mean(YYY);
        %
        [RR_Spe1,pp_Spe1] = corr(XXX1,YYY,'Type','Spearman');
        [RR_Pea1,pp_Pea1] = corr(XXX1,YYY,'Type','Pearson');
        [RR_Spe2,pp_Spe2] = corr(XXX2,YYY,'Type','Spearman');
        [RR_Pea2,pp_Pea2] = corr(XXX2,YYY,'Type','Pearson');
        %
        [~, p_tt1] = ttest(XXX1,YYY);
        [p_wil1,~] = signrank(XXX1,YYY);
        [~, p_tt2] = ttest(XXX2,YYY);
        [p_wil2,~] = signrank(XXX2,YYY);
        %
        % ancova
        YYY_all = [YYY; YYY]; % rho_neuronal
        XXX_all = [XXX1; XXX2]; % rho_JC; rho_SO
        taskgroup = [repelem({'Task1'},length(XXX1)),repelem({'Task2'},length(XXX2))]';
        [h,atab,ctab,stats] = aoctool(XXX_all,YYY_all,taskgroup,0.05,'','','','off','separate lines');
        p_slp = atab{4,6};
        %
        hold on        
        plot(XXX1, YYY, 'bo','MarkerSize',8);
        plot(XXX2, YYY, 'ro','MarkerSize',8);
        plot(XX,Yfit1,'b-','LineWidth',3);
        plot(XX,Yfit2,'r-','LineWidth',3);
        plotErrorEllipse(mu_ell1, Sigma_ell1, 0.90, 'b')
        plotErrorEllipse(mu_ell2, Sigma_ell2, 0.90, 'r')
        plot([minXY maxXY], [minXY maxXY],'--','Color', [0.5 0.5 0.5],'LineWidth',1);   
        ylabel(['\rho(neuronal) ',labeltype]);
        xlabel(['\rho(behavioral)']);
        legend({'Task 1','Task 2'});
        text(minXY+(maxXY-minXY)/2,minXY+(maxXY-minXY)/7,...
            {['Task 1: y = ',num2str(aa1(1),'%.2f'),'+',num2str(aa1(2),'%.2f'),'x'];...
             ['Task 1 slope: p (null=1) =',num2str(paa2_1,'%1.1g')]; ...
             ['Task 1 intercept: p (null=0) =',num2str(paa1_1,'%1.1g')]; ...
             ['Task 1 Spearman: r=',num2str(RR_Spe1,'%1.1g'), ', p=',num2str(pp_Spe1,'%1.1g')]; ...
             ['Task 1 Pearson: r=',num2str(RR_Pea1,'%1.1g'),', p=',num2str(pp_Pea1,'%1.1g')];...
             ['Task 1 t test: p=',num2str(p_tt1,'%1.1g')];...
             ['Task 1 Wilcoxon: p=',num2str(p_wil1,'%1.1g')];...
             ['Task 2: y = ',num2str(aa2(1),'%.2f'),'+',num2str(aa2(2),'%.2f'),'x'];...
             ['Task 2 slope: p (null=1) =',num2str(paa2_2,'%1.1g')]; ...
             ['Task 2 intercept: p (null=0) =',num2str(paa1_2,'%1.1g')]; ...
             ['Task 2 Spearman: r=',num2str(RR_Spe2,'%1.1g'), ', p=',num2str(pp_Spe2,'%1.1g')]; ...
             ['Task 2 Pearson: r=',num2str(RR_Pea2,'%1.1g'),', p=',num2str(pp_Pea2,'%1.1g')];...
             ['Task 2 t test: p=',num2str(p_tt2,'%1.1g')];...
             ['Task 2 Wilcoxon: p=',num2str(p_wil2,'%1.1g')];...
             ['slope difference: p (ANCOVA) =',num2str(p_slp,'%1.1g')]; ...
             ['N = ', num2str(length(YYY)),'cells']; ...
            }, 'fontsize', 9); 
        box off
        axis([minXY maxXY minXY maxXY])
        axis square     
        set(gca,'FontSize',14)
    end % for itype
    
    % 
    if 0
    % plot rho neuronal vs rho task1 and task2 together - merged postoffer1 and postoffer2
    h =  findobj('type','figure');
    %
    figure(length(h)+1);
    set(gcf,'position',[110 65 850 850], 'PaperPositionMode','auto')
    %          
    ind_good_SOoff1 =  nonzeros_QAB_inSO.inSOoff1(:,1).*nonzeros_QAB_inSO.inSOoff1(:,2) == 1;
    ind_good_SOoff2 =  nonzeros_QAB_inSO.inSOoff2(:,1).*nonzeros_QAB_inSO.inSOoff2(:,2) == 1;
    ind_good = [ind_good_SOoff1; ind_good_SOoff2];
    ind_good = ind_good & [ind_noneout;ind_noneout];
    ind_cosislope_SOoff1 = slopes_QAB_inSO.inSOoff1(:,1).*slopes_QAB_inSO.inSOoff1(:,2) > 0;
    ind_cosislope_SOoff2 = slopes_QAB_inSO.inSOoff2(:,1).*slopes_QAB_inSO.inSOoff2(:,2) > 0;
    ind_cosislope = [ind_cosislope_SOoff1; ind_cosislope_SOoff2];
    XXX1 = ([rhos_inJC; rhos_inJC]); 
    XXX1 = XXX1(ind_good & ind_cosislope);
    XXX2 = ([rhos_inSO; rhos_inSO]); 
    XXX2 = XXX2(ind_good & ind_cosislope);
    sloperatios_SOoff1 = (slopes_QAB_inSO.inSOoff1(:,1)./slopes_QAB_inSO.inSOoff1(:,2));
    sloperatios_SOoff2 = (slopes_QAB_inSO.inSOoff2(:,1)./slopes_QAB_inSO.inSOoff2(:,2));
    YYY = ([sloperatios_SOoff1; sloperatios_SOoff2]);
    YYY = YYY(ind_good & ind_cosislope);
    %
%     % remove outlier
    kout = 1.5;
    quantiles_XXX1 = quantile(XXX1,[0.25 0.5 0.75]);
    IQR_XXX1 = quantiles_XXX1(3) - quantiles_XXX1(1);
    outlier_XXX1 = [quantiles_XXX1(1)-kout*IQR_XXX1, quantiles_XXX1(3)+kout*IQR_XXX1];
    ind_goodXXX1 = (XXX1 > outlier_XXX1(1)  & XXX1 < outlier_XXX1(2));
    quantiles_XXX2 = quantile(XXX2,[0.25 0.5 0.75]);
    IQR_XXX2 = quantiles_XXX2(3) - quantiles_XXX2(1);
    outlier_XXX2 = [quantiles_XXX2(1)-kout*IQR_XXX2, quantiles_XXX2(3)+kout*IQR_XXX2];
    ind_goodXXX2 = (XXX2 > outlier_XXX2(1)  & XXX2 < outlier_XXX2(2));
    quantiles_YYY = quantile(YYY,[0.25 0.5 0.75]);
    IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
    outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
    ind_goodYYY = (YYY > outlier_YYY(1)  & YYY < outlier_YYY(2));
    ind_goodXY = ind_goodXXX1 & ind_goodXXX2 & ind_goodYYY;
%     ind_goodXY = logical(ones(size(XXX1)));
    %
    XXX1 = XXX1(ind_goodXY);
    XXX2 = XXX2(ind_goodXY);
    YYY = YYY(ind_goodXY);
    %
    anacellnum =  sum(length(XXX1)); 
    %
    maxXY = ceil(max([XXX1; XXX2; YYY]));
    minXY = floor(min([XXX1; XXX2; YYY]));
    %
    mdl1 = fitlm(XXX1,YYY);
    aa1 = mdl1.Coefficients.Estimate;
    paa1_1 = coefTest(mdl1,[1,0],0);
    paa2_1 = coefTest(mdl1,[0,1],1);
    mdl2 = fitlm(XXX2,YYY);
    aa2 = mdl2.Coefficients.Estimate;
    paa1_2 = coefTest(mdl2,[1,0],0);
    paa2_2 = coefTest(mdl2,[0,1],1);
    XX = [minXY maxXY];
    Yfit1 = aa1(2)*XX+aa1(1);
    Yfit2 = aa2(2)*XX+aa2(1);
    %
    Sigma_ell1 = cov(XXX1, YYY);
    mu_ell1(1) = mean(XXX1);
    mu_ell1(2) = mean(YYY);
    Sigma_ell2 = cov(XXX2, YYY);
    mu_ell2(1) = mean(XXX2);
    mu_ell2(2) = mean(YYY);
    %
    [RR_Spe1,pp_Spe1] = corr(XXX1,YYY,'Type','Spearman');
    [RR_Pea1,pp_Pea1] = corr(XXX1,YYY,'Type','Pearson');
    [RR_Spe2,pp_Spe2] = corr(XXX2,YYY,'Type','Spearman');
    [RR_Pea2,pp_Pea2] = corr(XXX2,YYY,'Type','Pearson');
    %
    % ancova
    YYY_all = [YYY; YYY]; % rho_neuronal
    XXX_all = [XXX1; XXX2]; % rho_JC; rho_SO
    taskgroup = [repelem({'Task1'},length(XXX1)),repelem({'Task2'},length(XXX2))]';
    [h,atab,ctab,stats] = aoctool(XXX_all,YYY_all,taskgroup,0.05,'','','','off','separate lines');
    p_slp = atab{4,6};
    %
    hold on        
    plot(XXX1, YYY, 'bo','MarkerSize',8);
    plot(XXX2, YYY, 'ro','MarkerSize',8);
    plot(XX,Yfit1,'b-','LineWidth',3);
    plot(XX,Yfit2,'r-','LineWidth',3);
    plotErrorEllipse(mu_ell1, Sigma_ell1, 0.90, 'b')
    plotErrorEllipse(mu_ell2, Sigma_ell2, 0.90, 'r')
    plot([minXY maxXY], [minXY maxXY],'--','Color', [0.5 0.5 0.5],'LineWidth',1);   
    ylabel([' \rho(neuronal)']);
    xlabel([' \rho(behavioral)']);
    legend({'Task 1','Task 2'});
    text(minXY+(maxXY-minXY)/2,minXY+(maxXY-minXY)/7,...
        {['Task 1: y = ',num2str(aa1(1),'%.2f'),'+',num2str(aa1(2),'%.2f'),'x'];...
         ['Task 1 slope: p (null=1) =',num2str(paa2_1,'%1.1g')]; ...
         ['Task 1 intercept: p (null=0) =',num2str(paa1_1,'%1.1g')]; ...
         ['Task 1 Spearman: r=',num2str(RR_Spe1,'%1.1g'), ', p=',num2str(pp_Spe1,'%1.1g')]; ...
         ['Task 1 Pearson: r=',num2str(RR_Pea1,'%1.1g'),', p=',num2str(pp_Pea1,'%1.1g')];...
         ['Task 2: y = ',num2str(aa2(1),'%.2f'),'+',num2str(aa2(2),'%.2f'),'x'];...
         ['Task 2 slope: p (null=1) =',num2str(paa2_2,'%1.1g')]; ...
         ['Task 2 intercept: p (null=0) =',num2str(paa1_2,'%1.1g')]; ...
         ['Task 2 Spearman: r=',num2str(RR_Spe2,'%1.1g'), ', p=',num2str(pp_Spe2,'%1.1g')]; ...
         ['Task 2 Pearson: r=',num2str(RR_Pea2,'%1.1g'),', p=',num2str(pp_Pea2,'%1.1g')];...
         ['slope difference: p (ANCOVA) =',num2str(p_slp,'%1.1g')]; ...
        }, 'fontsize', 9); 
    box off
    axis([minXY maxXY minXY maxXY])
    axis square     
    set(gca,'FontSize',14)
    %
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname]; ['Task 2']; ['postoffer1 and postoffer2 merged']},'fontsize',11);
    axis off
    end
    
    
    % % % 
    % rho neuronal v.s. order bias
    if 0
    %
    figure;
    set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname]; ['\rho(neuronal) and order bias \epsilon']},'fontsize',11);
    axis off
    %
    % three comparison; three correlation               
    xplottypes = {'inSOoff1', 'inAB12', 'inABpj'};
    yplottypes = {'inSOoff2', 'inBA12', 'inBApj'};
    xlabeltypes = {'postoffer1', 'AB postoffer', 'AB postjuice'};
    ylabeltypes = {'postoffer2', 'BA postoffer', 'BA postjuice'};
    ntypes = length(xplottypes);
    for itype = 1:ntypes       
        xplottype = xplottypes{itype};
        yplottype = yplottypes{itype};
        xlabeltype = xlabeltypes{itype};
        ylabeltype = ylabeltypes{itype};
        %       
        subplot(2,ntypes,itype);
        eval(['ind_goodx = nonzeros_QAB_inSO.',xplottype,'(:,1).*nonzeros_QAB_inSO.',xplottype,'(:,2) == 1;']) 
        eval(['ind_goody = nonzeros_QAB_inSO.',yplottype,'(:,1).*nonzeros_QAB_inSO.',yplottype,'(:,2) == 1;']) 
        eval(['ind_cosislopex = slopes_QAB_inSO.',xplottype,'(:,1).*slopes_QAB_inSO.',xplottype,'(:,2) > 0;'])
        eval(['ind_cosislopey = slopes_QAB_inSO.',yplottype,'(:,1).*slopes_QAB_inSO.',yplottype,'(:,2) > 0;'])       
        eval(['rho_XXX = log(slopes_QAB_inSO.',xplottype,'(:,1)./slopes_QAB_inSO.',xplottype,'(:,2));']) 
        eval(['rho_YYY = log(slopes_QAB_inSO.',yplottype,'(:,1)./slopes_QAB_inSO.',yplottype,'(:,2));']) 
        ind_good = ind_goodx & ind_goody;
        if ~dononzeroslopes, ind_good = logical(ones(size(ind_good))); end
        ind_good = ind_good & ind_noneout;
        ind_cosislope = ind_cosislopex & ind_cosislopey;
        if ~donconsislopes, ind_cosislope = logical(ones(size(ind_cosislope))); end
        %
        XXX = (rho_XXX(ind_good & ind_cosislope,:));
        YYY = (rho_YYY(ind_good & ind_cosislope,:));
        % 
        % remove outlier
        kout = 1.5;
        quantiles_XXX = quantile(XXX,[0.25 0.5 0.75]);
        IQR_XXX = quantiles_XXX(3) - quantiles_XXX(1);
        outlier_XXX = [quantiles_XXX(1)-kout*IQR_XXX, quantiles_XXX(3)+kout*IQR_XXX];
        ind_goodXXX = (XXX > outlier_XXX(1)  & XXX < outlier_XXX(2));
        quantiles_YYY = quantile(YYY,[0.25 0.5 0.75]);
        IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
        outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
        ind_goodYYY = (YYY > outlier_YYY(1)  & YYY < outlier_YYY(2));
        ind_goodXY = ind_goodXXX & ind_goodYYY;
        %
        XXX = XXX(ind_goodXY);
        YYY = YYY(ind_goodXY);
        %
        anacellnum =  sum(length(XXX)); 
        %
        maxXY = ceil(max([XXX; YYY]));
        minXY = floor(min([XXX; YYY]));
        hold on
        plot([minXY maxXY], [minXY maxXY],'--','Color', [0.8 0.8 0.8],'LineWidth',1);   
        aa = deming(XXX, YYY);
        XX = [minXY maxXY];
        Yfit = aa(2)*XX+aa(1);
        [~, p_ttest] = ttest(XXX, YYY);
        p_wil = signrank(XXX, YYY);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        plot(XXX, YYY, 'ko','MarkerSize',8);
        plot(XX, Yfit, '-','Color', [0.4 0.4 0.4],'LineWidth',3);
        xlabel(['log \rho(neuronal) of ',xlabeltype]);
        ylabel(['log \rho(neuronal) of ',ylabeltype]);
        text(minXY+(maxXY-minXY)/10,maxXY-(maxXY-minXY)/10,...
            {['t test: p=',num2str(p_ttest,'%1.1g')];['Wilcoxon: p=',num2str(p_wil,'%1.1g')];...
             ['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
         
        box off
        axis([minXY maxXY minXY maxXY])
        axis square     
        set(gca,'FontSize',14)
        
        % % %
        %       
        subplot(2,ntypes,itype+ntypes);
        eval(['ind_goodx = nonzeros_QAB_inSO.',xplottype,'(:,1).*nonzeros_QAB_inSO.',xplottype,'(:,2) == 1;']) 
        eval(['ind_goody = nonzeros_QAB_inSO.',yplottype,'(:,1).*nonzeros_QAB_inSO.',yplottype,'(:,2) == 1;']) 
        eval(['ind_cosislopex = slopes_QAB_inSO.',xplottype,'(:,1).*slopes_QAB_inSO.',xplottype,'(:,2) > 0;'])
        eval(['ind_cosislopey = slopes_QAB_inSO.',yplottype,'(:,1).*slopes_QAB_inSO.',yplottype,'(:,2) > 0;'])       
        eval(['rho_XXX = log(slopes_QAB_inSO.',xplottype,'(:,1)./slopes_QAB_inSO.',xplottype,'(:,2));']) 
        eval(['rho_YYY = log(slopes_QAB_inSO.',yplottype,'(:,1)./slopes_QAB_inSO.',yplottype,'(:,2));']) 
        ind_good = ind_goodx & ind_goody;
        if ~dononzeroslopes, ind_good = logical(ones(size(ind_good))); end
        ind_good = ind_good & ind_noneout;
        ind_cosislope = ind_cosislopex & ind_cosislopey;
        if ~doconsislopes, ind_cosislope = logical(ones(size(ind_cosislope))); end
        %
        XXX = orderbias_inSO(ind_good & ind_cosislope,:);
        YYY = rho_YYY(ind_good & ind_cosislope,:) - rho_XXX(ind_good & ind_cosislope,:);
        anacellnum =  sum(length(XXX));    
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
        plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
        hold on; plot(XXX, YYY, 'ko','MarkerSize',8);
        xlabel(['order bias \epsilon in Task 2']);
        ylabel(['\Delta log \rho(neuronal) (',ylabeltype,'-',xlabeltype,')']);
        text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
            {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
        box off
        axis([XX YY])
        axis square     
        set(gca,'FontSize',14)
    end % for itype
    end
    
    % % % 
    % rho neuronal v.s. preference bias
    %
    figure;
    set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname]; ['\rho(neuronal) and preference bias']},'fontsize',11);
    axis off
    %
    % three comparison; three correlation               
    xplottypes = {'inSOoff1', 'inAB12', 'inABpj'};
    yplottypes = {'inSOoff2', 'inBA12', 'inBApj'};
    xlabeltypes = {'postoffer1', 'AB postoffer', 'AB postjuice'};
    ylabeltypes = {'postoffer2', 'BA postoffer', 'BA postjuice'};
    ntypes = length(xplottypes);
    for itype = 1:ntypes       
        xplottype = xplottypes{itype};
        yplottype = yplottypes{itype};
        xlabeltype = xlabeltypes{itype};
        ylabeltype = ylabeltypes{itype};
        %       
        subplot(2,ntypes,itype);
        eval(['ind_goodx = nonzeros_QAB_inSO.',xplottype,'(:,1).*nonzeros_QAB_inSO.',xplottype,'(:,2) == 1;']) 
        eval(['ind_goody = nonzeros_QAB_inSO.',yplottype,'(:,1).*nonzeros_QAB_inSO.',yplottype,'(:,2) == 1;']) 
        eval(['ind_cosislopex = slopes_QAB_inSO.',xplottype,'(:,1).*slopes_QAB_inSO.',xplottype,'(:,2) > 0;'])
        eval(['ind_cosislopey = slopes_QAB_inSO.',yplottype,'(:,1).*slopes_QAB_inSO.',yplottype,'(:,2) > 0;'])       
        eval(['rho_XXX = (slopes_QAB_inSO.',xplottype,'(:,1)./slopes_QAB_inSO.',xplottype,'(:,2));']) 
        eval(['rho_YYY = (slopes_QAB_inSO.',yplottype,'(:,1)./slopes_QAB_inSO.',yplottype,'(:,2));']) 
        ind_good = ind_goodx & ind_goody;
        if ~dononzeroslopes, ind_good = logical(ones(size(ind_good))); end
        ind_good = ind_good & ind_noneout;
        ind_cosislope = ind_cosislopex & ind_cosislopey;
        if ~doconsislopes, ind_cosislope = logical(ones(size(ind_cosislope))); end
        %
        XXX = (rho_XXX(ind_good & ind_cosislope,:));
        YYY = (rho_YYY(ind_good & ind_cosislope,:));
        %
        % remove outlier
        kout = 1.5;
        quantiles_XXX = quantile(XXX,[0.25 0.5 0.75]);
        IQR_XXX = quantiles_XXX(3) - quantiles_XXX(1);
        outlier_XXX = [quantiles_XXX(1)-kout*IQR_XXX, quantiles_XXX(3)+kout*IQR_XXX];
        ind_goodXXX = (XXX > outlier_XXX(1)  & XXX < outlier_XXX(2));
        quantiles_YYY = quantile(YYY,[0.25 0.5 0.75]);
        IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
        outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
        ind_goodYYY = (YYY > outlier_YYY(1)  & YYY < outlier_YYY(2));
        ind_goodXY = ind_goodXXX & ind_goodYYY;
        %
        XXX = XXX(ind_goodXY);
        YYY = YYY(ind_goodXY);
        %
        anacellnum =  sum(length(XXX)); 
        %
        maxXY = ceil(max([XXX; YYY]));
        minXY = floor(min([XXX; YYY]));
        hold on
        plot([minXY maxXY], [minXY maxXY],'--','Color', [0.8 0.8 0.8],'LineWidth',1);   
        aa = deming(XXX, YYY);
        XX = [minXY maxXY];
        Yfit = aa(2)*XX+aa(1);
        [~, p_ttest] = ttest(XXX, YYY);
        p_wil = signrank(XXX, YYY);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        plot(XXX, YYY, 'ko','MarkerSize',8);
        plot(XX, Yfit, '-','Color', [0.4 0.4 0.4],'LineWidth',3);
        xlabel(['\rho(neuronal) of ',xlabeltype]);
        ylabel(['\rho(neuronal) of ',ylabeltype]);
        text(minXY+(maxXY-minXY)/10,maxXY-(maxXY-minXY)/10,...
            {['t test: p=',num2str(p_ttest,'%1.1g')];['Wilcoxon: p=',num2str(p_wil,'%1.1g')];...
             ['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
         
        box off
        % axis([minXY maxXY minXY maxXY])
        axis([minXY maxXY minXY maxXY])
        axis square     
        set(gca,'FontSize',14)
        
        % % %
        %       
        subplot(2,ntypes,itype+ntypes);
        eval(['ind_goodx = nonzeros_QAB_inSO.',xplottype,'(:,1).*nonzeros_QAB_inSO.',xplottype,'(:,2) == 1;']) 
        eval(['ind_goody = nonzeros_QAB_inSO.',yplottype,'(:,1).*nonzeros_QAB_inSO.',yplottype,'(:,2) == 1;']) 
        eval(['ind_cosislopex = slopes_QAB_inSO.',xplottype,'(:,1).*slopes_QAB_inSO.',xplottype,'(:,2) > 0;'])
        eval(['ind_cosislopey = slopes_QAB_inSO.',yplottype,'(:,1).*slopes_QAB_inSO.',yplottype,'(:,2) > 0;'])       
        eval(['rho_XXX = (slopes_QAB_inSO.',xplottype,'(:,1)./slopes_QAB_inSO.',xplottype,'(:,2));']) 
        eval(['rho_YYY = (slopes_QAB_inSO.',yplottype,'(:,1)./slopes_QAB_inSO.',yplottype,'(:,2));']) 
        ind_good = ind_goodx & ind_goody;
        if ~dononzeroslopes, ind_good = logical(ones(size(ind_good))); end
        ind_good = ind_good & ind_noneout;
        ind_cosislope = ind_cosislopex & ind_cosislopey;
        if ~doconsislopes, ind_cosislope = logical(ones(size(ind_cosislope))); end
        %
        % % remove outlier
        XXX_out = (rho_XXX(ind_good & ind_cosislope,:));
        YYY_out = (rho_YYY(ind_good & ind_cosislope,:));        
        %     
        kout = 1.5;
        quantiles_XXX = quantile(XXX_out,[0.25 0.5 0.75]);
        IQR_XXX = quantiles_XXX(3) - quantiles_XXX(1);
        outlier_XXX = [quantiles_XXX(1)-kout*IQR_XXX, quantiles_XXX(3)+kout*IQR_XXX];
        ind_goodXXX = (rho_XXX > outlier_XXX(1)  & rho_XXX < outlier_XXX(2));
        quantiles_YYY = quantile(YYY_out,[0.25 0.5 0.75]);
        IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
        outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
        ind_goodYYY = (rho_YYY > outlier_YYY(1)  & rho_YYY < outlier_YYY(2));
        ind_goodXY = ind_goodXXX & ind_goodYYY;
        % %
        % ind_goodXY = logical(ones(size(rho_XXX)));
        
        %
        XXX = (rhos_inJC(ind_good & ind_cosislope & ind_goodXY,:));
%         XXX = ((rhos_inSO(ind_good & ind_cosislope & ind_goodXY,:)) - (rhos_inJC(ind_good & ind_cosislope & ind_goodXY,:)))./...
%               ((rhos_inSO(ind_good & ind_cosislope & ind_goodXY,:)) + (rhos_inJC(ind_good & ind_cosislope & ind_goodXY,:)));
        YYY = (rho_YYY(ind_good & ind_cosislope & ind_goodXY,:))-(rho_XXX(ind_good & ind_cosislope & ind_goodXY,:));
        anacellnum =  sum(length(XXX));    
        % % 
        aa_mdl1 = fitlm(XXX,YYY);
        aa1 = aa_mdl1.Coefficients.Estimate;
        aa_mdl2 = fitlm(YYY,XXX);
        aa2 = aa_mdl2.Coefficients.Estimate;
        aa2 = [(-aa2(1)/aa2(2));1/aa2(2)];
        aa = aa1;
        % % 
        XX = [floor(min(XXX)),ceil(max(XXX))];
%         XX = [-0.2,0.2];
        YY = [floor(min(YYY)),ceil(max(YYY))]; 
        Yfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
        hold on; plot(XXX, YYY, 'ko','MarkerSize',8);
        xlabel(['\rho in Task 1']);
%         xlabel(['normalized \Delta\rho']);
        ylabel(['\Delta\rho(neuronal) (',ylabeltype,' - ',xlabeltype,')']);
        text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
            {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
        box off
        axis([XX YY])
        axis square     
        set(gca,'FontSize',14)
    end % for itype
    
    
    
    
    % % % 
    % rho neuronal: Task 1 and Task 2
    %
    figure;
    set(gcf,'position',[110 65 1650 850], 'PaperPositionMode','auto')
    axes('position',[.02 .97 .2 .05]);
    text(0,0,{[classname,', ',slopesignname]; ['\rho(neuronal) Task 1 and Task 2']},'fontsize',11);
    axis off
    %
    % two comparison       
    xplottypes = {'postoffer', 'postoffer', 'postoffer', 'postoffer'};
    % xplottypes = {'latedelay', 'latedelay', 'latedelay', 'latedelay'};
    % xplottypes = {'postjuice', 'postjuice', 'postjuice', 'postjuice'};
    yplottypes = {'inSOoff1',  'inSOoff2',  'inAB12',    'inBA12'};
    xlabeltypes = {'JC postoffer',  'JC postoffer',  'JC postoffer', 'JC postoffer'};
    % xlabeltypes = {'JC latedelay',  'JC latedelay',  'JC latedelay', 'JC latedelay'};
    % xlabeltypes = {'JC postjuice',  'JC postjuice',  'JC postjuice', 'JC postjuice'};
    ylabeltypes = {'SO postoffer1', 'SO postoffer2', 'AB postoffer', 'BA postoffer'};
    ntypes = length(xplottypes);
    for itype = 1:ntypes       
        xplottype = xplottypes{itype};
        yplottype = yplottypes{itype};
        xlabeltype = xlabeltypes{itype};
        ylabeltype = ylabeltypes{itype};
        %       
        subplot(2,ntypes,itype);
        eval(['ind_goodx = nonzeros_QAB_inJC.',xplottype,'(:,1).*nonzeros_QAB_inJC.',xplottype,'(:,2) == 1;']) 
        eval(['ind_goody = nonzeros_QAB_inSO.',yplottype,'(:,1).*nonzeros_QAB_inSO.',yplottype,'(:,2) == 1;']) 
        eval(['ind_cosislopex = slopes_QAB_inJC.',xplottype,'(:,1).*slopes_QAB_inJC.',xplottype,'(:,2) > 0;'])
        eval(['ind_cosislopey = slopes_QAB_inSO.',yplottype,'(:,1).*slopes_QAB_inSO.',yplottype,'(:,2) > 0;'])       
        eval(['rho_XXX = (slopes_QAB_inJC.',xplottype,'(:,1)./slopes_QAB_inJC.',xplottype,'(:,2));']) 
        eval(['rho_YYY = (slopes_QAB_inSO.',yplottype,'(:,1)./slopes_QAB_inSO.',yplottype,'(:,2));']) 
        ind_good = ind_goodx & ind_goody;
        if ~dononzeroslopes, ind_good = logical(ones(size(ind_good))); end
        ind_good = ind_good & ind_noneout;
        ind_cosislope = ind_cosislopex & ind_cosislopey;
        if ~doconsislopes, ind_cosislope = logical(ones(size(ind_cosislope))); end
        %
        XXX = (rho_XXX(ind_good & ind_cosislope,:));
        YYY = (rho_YYY(ind_good & ind_cosislope,:));        
        %    
        % % remove outlier
        kout = 1.5;
        quantiles_XXX = quantile(XXX,[0.25 0.5 0.75]);
        IQR_XXX = quantiles_XXX(3) - quantiles_XXX(1);
        outlier_XXX = [quantiles_XXX(1)-kout*IQR_XXX, quantiles_XXX(3)+kout*IQR_XXX];
        ind_goodXXX = (XXX > outlier_XXX(1)  & XXX < outlier_XXX(2));
        quantiles_YYY = quantile(YYY,[0.25 0.5 0.75]);
        IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
        outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
        ind_goodYYY = (YYY > outlier_YYY(1)  & YYY < outlier_YYY(2));
        ind_goodXY = ind_goodXXX & ind_goodYYY;
        %
        XXX = XXX(ind_goodXY);
        YYY = YYY(ind_goodXY);
        %
        anacellnum =  length(XXX);    
        maxXY = ceil(max([XXX; YYY]));
        minXY = floor(min([XXX; YYY]));
        hold on
        plot([minXY maxXY], [minXY maxXY],'--','Color', [0.8 0.8 0.8],'LineWidth',1);   
        %
        [aa,~,~,~,stat] = deming(XXX, YYY);
        b_ci = stat.b_ci;
%         aa_mdl1 = fitlm(XXX,YYY);
%         aa1 = aa_mdl1.Coefficients.Estimate;
%         aa_mdl2 = fitlm(YYY,XXX);
%         aa2 = aa_mdl2.Coefficients.Estimate;
%         aa2 = [(-aa2(1)/aa2(2));1/aa2(2)];
%         aa = [aa1+aa2]/2;
        %
        XX = [minXY maxXY];
        Yfit = aa(2)*XX+aa(1);
        [~, p_ttest] = ttest(XXX, YYY);
        p_wil = signrank(XXX, YYY);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        plot(XXX, YYY, 'ko','MarkerSize',8);
        plot(XX, Yfit, '-','Color', [0.4 0.4 0.4],'LineWidth',3);
        xlabel(['\rho(neuronal) of ',xlabeltype]);
        ylabel(['\rho(neuronal) of ',ylabeltype]);
        text(minXY+(maxXY-minXY)/10,maxXY-(maxXY-minXY)/10,...
            {['t test: p=',num2str(p_ttest,'%1.1g')];['Wilcoxon: p=',num2str(p_wil,'%1.1g')];...
             ['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
             ['slope (95% interval) = ',num2str(b_ci(2,1),'%.2f'),'~',num2str(b_ci(2,2),'%.2f')];...
             ['intercept (95% interval) = ',num2str(b_ci(1,1),'%.2f'),'~',num2str(b_ci(1,2),'%.2f')];...
             ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 9);
         
        box off
        % axis([minXY maxXY minXY maxXY])
        axis([minXY maxXY minXY maxXY])
        axis square     
        set(gca,'FontSize',14)
        
        % % %
        %       
        subplot(2,ntypes,itype+ntypes);
        eval(['ind_goodx = nonzeros_QAB_inJC.',xplottype,'(:,1).*nonzeros_QAB_inJC.',xplottype,'(:,2) == 1;']) 
        eval(['ind_goody = nonzeros_QAB_inSO.',yplottype,'(:,1).*nonzeros_QAB_inSO.',yplottype,'(:,2) == 1;']) 
        eval(['ind_cosislopex = slopes_QAB_inJC.',xplottype,'(:,1).*slopes_QAB_inJC.',xplottype,'(:,2) > 0;'])
        eval(['ind_cosislopey = slopes_QAB_inSO.',yplottype,'(:,1).*slopes_QAB_inSO.',yplottype,'(:,2) > 0;'])       
        eval(['rho_XXX = (slopes_QAB_inJC.',xplottype,'(:,1)./slopes_QAB_inJC.',xplottype,'(:,2));']) 
        eval(['rho_YYY = (slopes_QAB_inSO.',yplottype,'(:,1)./slopes_QAB_inSO.',yplottype,'(:,2));']) 
        ind_good = ind_goodx & ind_goody;
        if ~dononzeroslopes, ind_good = logical(ones(size(ind_good))); end
        ind_good = ind_good & ind_noneout;
        ind_cosislope = ind_cosislopex & ind_cosislopey;
        if ~doconsislopes, ind_cosislope = logical(ones(size(ind_cosislope))); end
        %
        % % remove outlier
        XXX_out = (rho_XXX(ind_good & ind_cosislope,:));
        YYY_out = (rho_YYY(ind_good & ind_cosislope,:));        
        %     
        kout = 1.5;
        quantiles_XXX = quantile(XXX_out,[0.25 0.5 0.75]);
        IQR_XXX = quantiles_XXX(3) - quantiles_XXX(1);
        outlier_XXX = [quantiles_XXX(1)-kout*IQR_XXX, quantiles_XXX(3)+kout*IQR_XXX];
        ind_goodXXX = (rho_XXX > outlier_XXX(1)  & rho_XXX < outlier_XXX(2));
        quantiles_YYY = quantile(YYY_out,[0.25 0.5 0.75]);
        IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
        outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
        ind_goodYYY = (rho_YYY > outlier_YYY(1)  & rho_YYY < outlier_YYY(2));
        ind_goodXY = ind_goodXXX & ind_goodYYY;
        % %
        % ind_goodXY = logical(ones(size(rho_XXX)));
        
        XXX = (rhos_inJC(ind_good & ind_cosislope & ind_goodXY,:));
%         XXX = ((rhos_inSO(ind_good & ind_cosislope & ind_goodXY,:)) - (rhos_inJC(ind_good & ind_cosislope & ind_goodXY,:)))./...
%               ((rhos_inSO(ind_good & ind_cosislope & ind_goodXY,:)) + (rhos_inJC(ind_good & ind_cosislope & ind_goodXY,:)));
        YYY = (rho_YYY(ind_good & ind_cosislope & ind_goodXY,:)) - (rho_XXX(ind_good & ind_cosislope & ind_goodXY,:));
        anacellnum = length(XXX);    
        % % 
        aa_mdl1 = fitlm(XXX,YYY);
        aa1 = aa_mdl1.Coefficients.Estimate;
        aa_mdl2 = fitlm(YYY,XXX);
        aa2 = aa_mdl2.Coefficients.Estimate;
        aa2 = [(-aa2(1)/aa2(2));1/aa2(2)];
        aa = aa1;
        % % 
        XX = [floor(min(XXX)),ceil(max(XXX))];
%         XX = [-0.2,0.2];
        YY = [floor(min(YYY)),ceil(max(YYY))]; 
        Yfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
        hold on; plot(XXX, YYY, 'ko','MarkerSize',8);
        xlabel(['\rho in Task 1']);
%         xlabel(['normalized \Delta\rho']);
        ylabel(['\rho(neuronal) (',ylabeltype,' - ',xlabeltype,')']);
        text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
            {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
        box off
        axis([XX YY])
        axis square     
        set(gca,'FontSize',14)
        
    end
   
end





% % % 
function plotErrorEllipse(mu_ell, Sigma_ell, p_ell, clr_codes)
    s = -2 * log(1 - p_ell);
    [V, D] = eig(Sigma_ell * s);
    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    plot(a(1, :) + mu_ell(1), a(2, :) + mu_ell(2), 'LineWidth',1.5, 'Color', clr_codes);
end
    
