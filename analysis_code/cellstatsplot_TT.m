function nplots = cellstatsplot_TT(cellname, cellstats)
%cellstatsplot
%
% plots all the cell stats (3-way anova + tuning fit)
%

% author: camillo
% revisions:	november 2009: now also displays anova bytrialtype
%				december 2009: adapted for DT
%				January 2017: adapted for SO
%               October 2018: adapted for TT -WS

%keyboard


if 0
	clear all
	% 	cellname = 'L091216c11';
	cellname = 'L100310a11';
	session = cellname(1:8); readsession
	filename = [dirroot, cellname, '_cellstats'];
	eval(['load ',filename])
end



disp(['   ... plotting stats of cell ', cellname])
fittoplot = 'tuningfit';
Rsqtouse = 'Rsq';
nonzero_confint_touse = 'p90';
fsize=9; limp=0.01;

JCSO = {'JC','SO'};

for iJCorSO = 1:length(JCSO)
    JCorSO = JCSO{iJCorSO};
    
    if isequal(JCorSO,'JC')
        pairnames = {'AB'};
        npairs = 1;
    elseif isequal(JCorSO,'SO')
        pairnames = {'AB','BA','ABA'};
        npairs = 3;
    end
    
    session = cellname(1:8); readsession_TT
    %#ok<NASGU>
    [timewindows,~] = get_timewindows_TT(JCorSO);
    filename = [dirroot, cellname, '_bhvParams']; eval(['load ',filename])
    
    ntwins = length(timewindows);
    
    anovastats = cellstats.anovastats;
    anovastats_both = cellstats.anovastats_both;
    tuningfit = cellstats.tuningfit;  %plot all models
    % tuningfit = cellstats.cellclass;    %plot only selected models

    for ipair = npairs:-1:1
    %pairname = pairs{ipair};
        pairname = pairnames{ipair};
	
		rnames = [];
		for iwin = 1:ntwins
			rnames = [rnames; timewindows{iwin}(1)];
        end
	
        if		isequal(pairname,pairnames{1}), pairnum = 1;
        elseif	isequal(pairname,pairnames{2}), pairnum = 2;
        elseif	isequal(pairname,pairnames{3}), pairnum = 3;
        end
	
        
        
        if isequal(JCorSO,'JC')
        
            %plot anovastat % two task types together
            hf(1)=figure(1);
            pairnum==1; 
            set(gca,'position',[.020 .71 .28 .25]);

            xx = [0,.24:.2:.65 .9];
            yy = [.855:-.08:0];
	
            text(xx(2),	yy(1),	'offtype',		'fontsize', fsize)
            text(xx(3),	yy(1),	'order',		'fontsize', fsize)
            text(xx(4),	yy(1),	'tasktype',			'fontsize', fsize)
            text(xx(5),	yy(1),	'trialtype',		'fontsize', fsize)
            timewindows_both = {'postoffer','postjuice'};
            for iwin = 1:2
                text(xx(1),	yy(iwin+1),	timewindows_both{iwin},	'fontsize', fsize)
            end
	
            %bytaskNtrialtype
            for ifactor = 1:3
                for iwin = 1:2
                    eval(['pp = anovastats_both.JCSO.bytaskNtrialtype.pval;']); pp = pp';
                    ppp = pp(ifactor,iwin);
                    ht = text(xx(ifactor+1),yy(iwin+1),sprintf('%0.4f',ppp),'fontsize', fsize-1);
                    if (ppp<.01),	set(ht,'color',[.7 0 0]); end
                    if (ppp<.001),	set(ht,'color',[1 0 0]); end
                end
            end
            %bytrialtype
            for iwin = 1:2
                eval(['pp = anovastats_both.JCSO.bytrialtype.pval;']); pp = pp';
                ppp = pp(1,iwin);
                ht = text(xx(5),yy(iwin+1),sprintf('%0.4f',ppp),'fontsize', fsize-1);
                if (ppp<.01),	set(ht,'color',[.7 0 0]); end
                if (ppp<.001),	set(ht,'color',[1 0 0]); end
            end
            axis off
            
            %add cellname
            axes('position',[.05 .925 .5 .017])
            text(0,3,['cell: ',cellname,' JC and SO together'],'fontweight','bold')
            text(0,1.5, 'pvalues of 3-way anova ', ...
                'fontweight','bold','fontsize', fsize)
            axis off
            
            
            
            
            %plot anovastat
            hf(2)=figure(2);
            pairnum==1; 
            set(gca,'position',[.020 .71 .28 .25]);

            xx = [0,.24:.2:.65 .9];
            yy = [.855:-.08:0];
	
            text(xx(2),	yy(1),	'offtype',		'fontsize', fsize)
            text(xx(3),	yy(1),	'pos',		'fontsize', fsize)
            text(xx(4),	yy(1),	'movdir',			'fontsize', fsize)
            text(xx(5),	yy(1),	'trtype',		'fontsize', fsize)
            for iwin = 1:ntwins
                text(xx(1),	yy(iwin+1),	timewindows{iwin}(1),	'fontsize', fsize)
            end
	
            %byofftypeNside
            for ifactor = 1:3
                for iwin = 1:ntwins
                    eval(['pp = anovastats.',JCorSO,'.',pairname,'.byoffertypeNside.pval;']); pp = pp';
                    ppp = pp(ifactor,iwin);
                    ht = text(xx(ifactor+1),yy(iwin+1),sprintf('%0.4f',ppp),'fontsize', fsize-1);
                    if (ppp<.01),	set(ht,'color',[.7 0 0]); end
                    if (ppp<.001),	set(ht,'color',[1 0 0]); end
                end
            end
            %bytrialtype
            eval(['pp = anovastats.',JCorSO,'.',pairname,'.bytrialtype.pval;']); pp = pp';
            pp_new = pp;
            eval(['pp_both = anovastats_both.JCSO.bytrialtype.pval;']); pp_both = pp_both';
            for iwin = 1:ntwins
               
                ppp = pp(1,iwin);
                twin = timewindows{iwin}{1};
                
                if isequal(twin,'postoffer')
                    pp_new(1,iwin) = min(ppp,pp_both(1,1));
                elseif isequal(twin,'postoffer1')
                    pp_new(1,iwin) = min(ppp,pp_both(1,1));
                elseif isequal(twin,'postoffer2')
                    pp_new(1,iwin) = min(ppp,pp_both(1,1));
                elseif isequal(twin,'postjuice')
                    pp_new(1,iwin) = min(ppp,pp_both(1,2));
                end
                
                
                ht = text(xx(5),yy(iwin+1),sprintf('%0.4f',ppp),'fontsize', fsize-1);
                if (ppp<.01),	set(ht,'color',[.7 0 0]); end
                if (ppp<.001),	set(ht,'color',[1 0 0]); end
            end
            axis off
	
            
            
            %
            %retrieve Rsq, slope, nonzero and modellabs
            eval(['Rsq = ',fittoplot,'.',JCorSO,'.',pairname,'.',Rsqtouse,';'])	
            eval(['slopes = ',fittoplot,'.',JCorSO,'.',pairname,'.slope;'])	
            eval(['nonzero = ',fittoplot,'.',JCorSO,'.',pairname,'.nonzero.',nonzero_confint_touse,';'])
            eval(['modellabs = ',fittoplot,'.',JCorSO,'.',pairname,'.modellabs;'])
            Rsq(isnan(Rsq)) = 0;
            slopes(isnan(slopes)) = 0;
            slopes=sign(slopes);
            nonzero(isnan(nonzero)) = 0;
        %
%             %select and reorder models
%             nmodels = size(Rsq,1);
%             if nmodels>4
%                 ind = reorder_models;
%                 Rsq = Rsq(ind,:);
%                 nonzero = nonzero(ind,:);
%                 modellabs = modellabs(ind);
%             end
            nmodels = size(Rsq,1);
            %
            %plot tuningfit
            if nmodels>4
                axesyy = [.4 .1 .0];
                axes('position',[.03 axesyy(pairnum) .9 .22])
                xx = [.02, .10:.9/nmodels:1];
                yy = [.855:-.08:0];
            else
                axesyy = [.4 .1 .0];
                axes('position',[.03 axesyy(pairnum) .9 .22])
                xx = [.02, .10:.9/15:1];
                yy = [.855:-.08:0];
            end
    
            ht = [];
            for iwin = 1:ntwins
                ht(iwin) = text(xx(1),	yy(iwin+1),	timewindows{iwin}(1),	'fontsize', fsize-1);
            end
            for iwin = 1:ntwins
                if (pp_new(1,iwin)<.01)		set(ht(iwin),'color',[.7 0 0]);	end
                if (pp_new(1,iwin)<.001)	set(ht(iwin),'color',[1 0 0]);	end
            end
            iwin = iwin+1;
            ht(iwin) = text(xx(1),	yy(iwin+1),	['sum'],	'fontsize', fsize-1);
             
            for imodel = 1:nmodels
                text(xx(imodel+1),yy(1),modellabs{imodel},'fontsize', fsize-2,'rotation',30);
            end
	
            for iwin = 1:ntwins+1
                ht = [];
                for imodel = 1:nmodels
                    if iwin < ntwins+1
                    %show only variables that explain the response
                    % if pp(1,iwin)<limp && nonzero(imodel,iwin)
                    if pp_new(1,iwin)<limp %&& nonzero(imodel,iwin)
                        txt = sprintf('%0.2f',Rsq(imodel,iwin));
                    else
                        if Rsq(imodel,iwin)>0
                            txt = sprintf('%0.2f',Rsq(imodel,iwin));
                        else
                            txt = '-';
                        end
                        % txt = sprintf('%0.2f',Rsq(imodel,iwin));
                    end
                    
                    sumRsq = sum(Rsq,2);
                    % sumRsq = sum(Rsq(:,[2,3,7]),2);
                    elseif iwin == ntwins+1
                        txt = sprintf('%0.2f',sumRsq(imodel,:));
                    end
                    
                    ht(imodel) = text(xx(imodel+1),yy(iwin+1),txt,'fontsize', fsize-1);
                end
                %add color codes
                if iwin < ntwins+1
                [Rmax,bestmod] = max(Rsq(:,iwin));
                if nonzero(bestmod,iwin)
                    if		(pp_new(1,iwin)<.001)	set(ht(bestmod),'color',[1 0 0]); 
                    elseif	(pp_new(1,iwin)<.01)	set(ht(bestmod),'color',[.7 0 0]); 
                    else						set(ht(bestmod),'color',[0 0 .7]); 
                    end
                end
                elseif iwin == ntwins+1
                    [Rmax,bestmod] = max(sumRsq(:,1));
                    set(ht(bestmod),'color',[1 0 0]); 
                end
            end
            axis off
            
            
            % for cell classes, only when all models are shown for tuningfit above
            if nmodels > 4
                axesyy = [.1 .1 .0];
                axes('position',[.03 axesyy(pairnum) .9 .22])
                xx = [.02, .10:.9/15:1];
                yy = [.855:-.08:0];
                
                % consider sign of slopes in Rsq
                Rsq = Rsq.*slopes;
                
                ht = [];
                twins_select = [2,3,6,7]; % postoffer latedelay prejuice postjuice
                ntwins_select = length(twins_select);
                
                for iwin = 1:ntwins_select
                ht(iwin) = text(xx(1),	yy(iwin+1),	timewindows{twins_select(iwin)}(1),	'fontsize', fsize-1);
                end
                iwin = iwin+1;
                ht(iwin) = text(xx(1),	yy(iwin+1),	['sum'],	'fontsize', fsize-1);
                
                models_select = [12 13 6 14]; % OVA OVB CV CJ
                nmodels_select = length(models_select);
                for imodel = 1:nmodels_select
                    text(xx(imodel+1),yy(1),modellabs{models_select(imodel)},'fontsize', fsize-2,'rotation',30);
                end
	
                for iwin = 1:ntwins_select+1
                    ht = [];
                    for imodel = 1:nmodels_select
                        if iwin < ntwins_select+1
                        if pp_new(1,twins_select(iwin))<limp %&& nonzero(imodel,iwin)
                            txt = sprintf('%0.2f',Rsq(models_select(imodel),twins_select(iwin)));
                        else
                            if abs(Rsq(models_select(imodel),twins_select(iwin)))>0
                                txt = sprintf('%0.2f',Rsq(models_select(imodel),twins_select(iwin)));
                            else
                                txt = '-';
                            end
                        end
                    
                        sumRsq = sum(Rsq(models_select,twins_select),2);
                        elseif iwin == ntwins_select+1
                            txt = sprintf('%0.2f',sumRsq(imodel,:));
                        end
                    
                        ht(imodel) = text(xx(imodel+1),yy(iwin+1),txt,'fontsize', fsize-1);
                    end
                end
                axis off
            end
            
            
            
            % for cell classes, only when all models are shown for
            % tuningfit above and consider slopes and nonzero or not 
            if nmodels > 4
                axesyy = [.1 .1 .0];
                axes('position',[.03 axesyy(pairnum) .9 .22])
                xx = [.02, .10:.9/15:1]+0.5;
                yy = [.855:-.08:0];
                
                % consider sign of slopes in Rsq
                Rsq = Rsq.*nonzero; % slopes has considered before
                
                ht = [];
                twins_select = [2,3,6,7]; % postoffer latedelay prejuice postjuice
                ntwins_select = length(twins_select);
                
                for iwin = 1:ntwins_select
                ht(iwin) = text(xx(1),	yy(iwin+1),	timewindows{twins_select(iwin)}(1),	'fontsize', fsize-1);
                end
                iwin = iwin+1;
                ht(iwin) = text(xx(1),	yy(iwin+1),	['sum'],	'fontsize', fsize-1);
                
                models_select = [12 13 6 14]; % OVA OVB CV CJ
                nmodels_select = length(models_select);
                for imodel = 1:nmodels_select
                    text(xx(imodel+1),yy(1),modellabs{models_select(imodel)},'fontsize', fsize-2,'rotation',30);
                end
	
                for iwin = 1:ntwins_select+1
                    ht = [];
                    for imodel = 1:nmodels_select
                        if iwin < ntwins_select+1
                        if pp_new(1,twins_select(iwin))<limp && abs(Rsq(models_select(imodel),twins_select(iwin)))>0
                            txt = sprintf('%0.2f',Rsq(models_select(imodel),twins_select(iwin)));
                        else
                            if abs(Rsq(models_select(imodel),twins_select(iwin)))>0
                                txt = sprintf('%0.2f',Rsq(models_select(imodel),twins_select(iwin)));
                            else
                                txt = '-';
                            end
                        end
                    
                        sumRsq = sum(Rsq(models_select,twins_select),2);
                        elseif iwin == ntwins_select+1
                            txt = sprintf('%0.2f',sumRsq(imodel,:));
                        end
                    
                        ht(imodel) = text(xx(imodel+1),yy(iwin+1),txt,'fontsize', fsize-1);
                    end
                end
                axis off
            end
            
            
            
            
            %add cellname
            axes('position',[.05 .925 .5 .017])
            text(0,3,['cell: ',cellname,' JC trials'],'fontweight','bold')
            text(0,1.5, 'pvalues of 3-way anova ', ...
                'fontweight','bold','fontsize', fsize)
            axis off
		
            if 		(pairnum==1) yytxt = .655;
            elseif 	(pairnum==2) yytxt = .355;
            elseif 	(pairnum==3) yytxt = .055;
            end
            starty=0.6;
            axes('position',[.05 starty+0.1 .5 .017])
            text(0,1.025, ['R^2 of linear models (without nonzero)'],'fontweight','bold','fontsize', fsize)
            axis off
            
            starty=0.6;
            axes('position',[.05 starty+0.1 .5 .017])
            text(1,-20.025, ['R^2 of linear models(with nonzero)'],'fontweight','bold','fontsize', fsize)
            axis off
		
            set(gcf,'Units','normalized','position',[0.5 0 0.5 1])
            set(gcf,'PaperPositionMode','auto')
        
        
        
        elseif ipair==3 & isequal(JCorSO,'SO')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ANOVA I
            %%%%%%%%%%%%%% ANOVA byOffertypeOrderChoosenJuiceOrder
            hf(3)=figure(3);
%           eval(['cc1 = anovastats.',pairname,'.byOffertypeOrderSide.tab.preoffer1(2:end-2,1);']);
%           cnames = [cc1; 'trialtypes'];
%           eval(['pp1 = anovastats.',pairname,'.byOffertypeOrderSide.pval;']);
%           eval(['pp2 = anovastats.',pairname,'.bytrialtype.pval;']);
%           pp=[pp1 pp2];

            eval(['cc1 = anovastats.',JCorSO,'.',pairname,'.byOffertypeOrderSideChoosen.tab.preoffers(2:end-2,1);']);
		
            try
                cnames = {'Offertypes','ChosenOrder','ChosenJuice','ChosenSide','Offtype*Order','Trialtypes'};
            catch
                for n=1:10
                    disp(['ERROR:' cellname])
                end
                    keyboard
                break
            end

          eval(['pp1 = anovastats.',JCorSO,'.',pairname,'.byOffertypeOrderSideChoosen.pval;']);
          eval(['pp2 = anovastats.',JCorSO,'.',pairname,'.bytrialtype.pval;']);
          eval(['pp3 = anovastats.',JCorSO,'.',pairname,'.byOfOrChinteractions.pval;']);
          pp=[pp1 pp3(:,1) pp2];
         

%             eval('pp = anovastats.SO.ABA.bytrialtype.pval;');

            starty=0.95; limy=0.75;
            stp=1.2/(size(cnames,2)+1);
            xx = [-.1:stp:1.1];
            stp=(starty-limy)/(size(rnames,1)+1);
            yy = [starty:-stp:limy];
            for ifactor = 1:size(pp,2)
                for iwin = 1:ntwins
                    ppp = pp(iwin,ifactor);
                    % if (ppp<.01)
                        ht = text(xx(ifactor+1),yy(iwin+1),sprintf('%0.4f',ppp),'fontsize', fsize);
                        if (ppp<.01),	set(ht,'color',[.7 0 0]); end
                        if (ppp<.001),	set(ht,'color',[1 0 0]); end
                    % else
                        % try
                        %     ht = text(xx(ifactor+1),yy(iwin+1),'--','fontsize', fsize-2);
                        % catch
                        %     % keyboard
                        % end
                   %  end
                end
            end
		
% 				keyboard 
				
            for ifactor = 1:size(pp,2)
                text(xx(ifactor+1),	yy(1),	cnames{ifactor},		'fontsize', fsize-2)
            end
            for iwin = 1:ntwins
                text(xx(1),	yy(iwin+1),	rnames{iwin},		'fontsize', fsize-1)
            end
            axis off

		
% 		
% 		%%%%%%%%%%%%%ANOVA bySideChoosenSide trialtype
% 		starty=limy-0.1; limy=limy-0.30;
% 		stp=1.2/(size(cnames,1)+1);
% 		xx = [-.1:stp:1.1];
% 		stp=(starty-limy)/(size(rnames,1)+1);
% 		yy = [starty:-stp:limy];
% 		eval(['cc1 = anovastats.',pairname,'.byOfOrChinteractions.tab.preoffers(2:end-2,1);']);
% 		cnames = [cc1];
% 		eval(['pp = anovastats.',pairname,'.byOfOrChinteractions.pval;']);
% 		
% 		for ifactor = 1:size(pp,2)
% 			for iwin = 1:ntwins
% 				ppp = pp(iwin,ifactor);
% 				if (ppp<.01)
% 					ht = text(xx(ifactor+1),yy(iwin+1),sprintf('%0.4f',ppp),'fontsize', fsize);
% 					if (ppp<.01),	set(ht,'color',[.7 0 0]); end
% 					if (ppp<.001),	set(ht,'color',[1 0 0]); end
% 				else
% 					ht = text(xx(ifactor+1),yy(iwin+1),'--','fontsize', fsize-2);
% 				end
% 			end
% 		end
% 		
		
		
            for ifactor = 1:size(pp,2)
			
                try
                    text(xx(ifactor+1),	yy(1),	cnames{ifactor},		'fontsize', fsize-2)
                catch 
                    keyboard 
                end
			
            end
            for iwin = 1:ntwins
                text(xx(1),	yy(iwin+1),	rnames{iwin},		'fontsize', fsize-1)
            end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %retrieve Rsq, nonzero and modellabs
		
%           eval(['pp1 = anovastats.',pairname,'.byOffertypeOrderSideChoosen.pval;']); 
%           eval(['pp2 = anovastats.',pairname,'.bytrialtype.pval;']);
%           eval(['pp3 = anovastats.',pairname,'.byOfOrChinteractions.pval;']);
%           pp=horzcat(pp1,pp2,pp3); pp=nanmin(pp'); pp=pp';

%           eval('pp1 = anovastats.ABA.bytrialtype.pval;');
%           eval('pp2 = anovastats.AB.bytrialtype.pval;');
%           eval('pp3 = anovastats.BA.bytrialtype.pval;');
%           pp=horzcat(pp1,pp2,pp3); pp=nanmin(pp'); pp=pp';

            eval('pp = anovastats.SO.ABA.bytrialtype.pval;');
            pp_new = pp;
            eval('pp_both = anovastats_both.JCSO.bytrialtype.pval;');
            for iwin = 1:ntwins
                twin = timewindows{iwin}{1};
                
                if isequal(twin,'postoffer')
                    pp_new(iwin) = min(pp(iwin),pp_both(1));
                elseif isequal(twin,'postoffer1')
                    pp_new(iwin) = min(pp(iwin),pp_both(1));
                elseif isequal(twin,'postoffer2')
                    pp_new(iwin) = min(pp(iwin),pp_both(1));
                elseif isequal(twin,'postjuice')
                    pp_new(iwin) = min(pp(iwin),pp_both(2));
                end
            end
            
            
            
            eval(['Rsq = ',fittoplot,'.',JCorSO,'.',pairname,'.',Rsqtouse,';'])
            eval(['slopes = ',fittoplot,'.',JCorSO,'.',pairname,'.slope;'])
            eval(['nonzero = ',fittoplot,'.',JCorSO,'.',pairname,'.nonzero.',nonzero_confint_touse,';'])
            eval(['modellabs = ',fittoplot,'.',JCorSO,'.',pairname,'.modellabs;'])
            Rsq(isnan(Rsq)) = 0;
            slopes(isnan(slopes)) = 0;
            slopes=sign(slopes);
            nonzero(isnan(nonzero)) = 0;
            %
            %select and reorder models
            nmodels = size(Rsq,1);

%           if nmodels>4
%               ind = reorder_models;
%               Rsq = Rsq(ind,:);
%               nonzero = nonzero(ind,:);
%               modellabs = modellabs(ind);
%               nmodels = size(Rsq,1);
%           end

            %plot tuningfit
            starty=limy-0.15; limy=limy-0.30;
            stp=1.2/(nmodels+2);
            xx = [-.15:stp:1.15];
            stp=(starty-limy)/(size(rnames,1)+1);
            yy = [starty:-stp:limy];
            ht = [];
            for iwin = 1:ntwins
                ht(iwin) = text(xx(1),	yy(iwin+1),	timewindows{iwin}(1),	'fontsize', fsize-1);
            end
            for iwin = 1:ntwins
                if (pp_new(iwin)<.01)	set(ht(iwin),'color',[.7 0 0]);	end
                if (pp_new(iwin)<.001)	set(ht(iwin),'color',[1 0 0]);	end
            end
            for imodel = 1:nmodels
                text(xx(imodel+2),yy(1),modellabs{imodel},'fontsize', fsize-2,'rotation',45)
            end
            for iwin = 1:ntwins
                ht = [];
                for imodel = 1:nmodels
                    %show only variables that explain the response
                    % if pp(iwin)<limp && nonzero(imodel,iwin)
                    if pp_new(iwin)<limp % && nonzero(imodel,iwin)
                        txt = sprintf('%0.2f',Rsq(imodel,iwin));
                    else
                        if Rsq(imodel,iwin)>0
                            txt = sprintf('%0.2f',Rsq(imodel,iwin));
                        else
                            txt = '-';
                        end
                        % txt = sprintf('%0.2f',Rsq(imodel,iwin));
                    end
                    ht(imodel) = text(xx(imodel+2),yy(iwin+1),txt,'fontsize', fsize-1);
                end
                %add color codes
                [Rmax,bestmod] = max(Rsq(:,iwin));
                if nonzero(bestmod,iwin)
                    if		(pp_new(iwin)<.001)	set(ht(bestmod),'color',[1 0 0]);
                    elseif	(pp_new(iwin)<.01)	set(ht(bestmod),'color',[.7 0 0]);
                    else					set(ht(bestmod),'color',[0 0 .7]);
                    end
                end
            end
            axis off
		
            % for cell classes (sequences), only when all models are shown for tuningfit above
            if nmodels > 4
                
                starty=limy-0.15; limy=limy-0.30;
                stp=1.2/(nmodels+2);
                xx = [-.15:stp:1.15];
                stp=(starty-limy)/(size(rnames,1)+1);
                yy = [starty:-stp:limy];
                ht = [];
                
                Rsq(20,2)=-Rsq(20,2); % reverse the sign of AB|BA in postoffer1 twins (compared with postoffer2 twins)
                % consider sign of slopes in Rsq
                Rsq = Rsq.*slopes;
                
                sequence_twins = [2 4 8];
                sequence_models = [ 2  3  16; %OVA
                                    6  5  17; %OVB
                                    7  9  15; %CV
                                    20 20 19; %CJ
                                    ];      
                sequence_labels = {'offer value A', 'offer value B', 'chosen value', 'chosen juice'};
                
                Rsq_sequence = [];
                ntwins_sequence = size(sequence_twins,2);
                nmodels_sequence = size(sequence_models,1);
                
                for iwin = 1:ntwins_sequence
                    ht(iwin) = text(xx(1),	yy(iwin+1),	timewindows{sequence_twins(iwin)}(1),	'fontsize', fsize-1);
                end
                    
                for imodel = 1:nmodels_sequence
                    text(xx(imodel+2),yy(1),sequence_labels{imodel},'fontsize', fsize-2,'rotation',45)
                end
            
                for iwin = 1:ntwins_sequence
                    ht=[];
                    for imodel = 1:nmodels_sequence
                        Rsq_sequence(imodel,iwin) = Rsq(sequence_models(imodel,iwin),sequence_twins(iwin));
                         if pp_new(sequence_twins(iwin))<limp
                             txt = sprintf('%0.2f',Rsq(sequence_models(imodel,iwin),sequence_twins(iwin)));
                         elseif abs(Rsq(sequence_models(imodel,iwin),sequence_twins(iwin)))>0
                             txt = sprintf('%0.2f',Rsq(sequence_models(imodel,iwin),sequence_twins(iwin)));
                         else
                             txt = '-';
                         end
                         ht(imodel) = text(xx(imodel+2),yy(iwin+1),txt,'fontsize', fsize-1);
                        
                    end
                end
                        
                iwin = iwin+1;
                ht=[];
                Rsq_sequence_sum = sum(Rsq_sequence,2);
                for imodel = 1:nmodels_sequence
                    txt = sprintf('%0.2f',Rsq_sequence_sum(imodel));
                    ht(imodel) = text(xx(imodel+2),yy(iwin+1),txt,'fontsize', fsize-1);
                end
            end
            axis off
            
            
            
             % for cell classes (sequences), only when all models are shown
             % for tuningfit above and consider slope and nonzero
            if nmodels > 4
                
                starty=limy-0.15; limy=limy-0.30;
                stp=1.2/(nmodels+2);
                xx = [-.15:stp:1.15];
                stp=(starty-limy)/(size(rnames,1)+1);
                yy = [starty:-stp:limy];
                ht = [];
                
                % Rsq(20,2)=-Rsq(20,2); % reverse the sign of AB|BA in postoffer1 twins (compared with postoffer2 twins) % has reversed already
                % consider sign of slopes in Rsq
                Rsq = Rsq.*nonzero; % slopes has considered before
                
                sequence_twins = [2 4 8];
                sequence_models = [ 2  3  16; %OVA
                                    6  5  17; %OVB
                                    7  9  15; %CV
                                    20 20 19; %CJ
                                    ];      
                sequence_labels = {'offer value A', 'offer value B', 'chosen value', 'chosen juice'};
                
                Rsq_sequence = [];
                ntwins_sequence = size(sequence_twins,2);
                nmodels_sequence = size(sequence_models,1);
                
                for iwin = 1:ntwins_sequence
                    ht(iwin) = text(xx(1),	yy(iwin+1),	timewindows{sequence_twins(iwin)}(1),	'fontsize', fsize-1);
                end
                    
                for imodel = 1:nmodels_sequence
                    text(xx(imodel+2),yy(1),sequence_labels{imodel},'fontsize', fsize-2,'rotation',45)
                end
            
                for iwin = 1:ntwins_sequence
                    ht=[];
                    for imodel = 1:nmodels_sequence
                        Rsq_sequence(imodel,iwin) = Rsq(sequence_models(imodel,iwin),sequence_twins(iwin));
                         if pp_new(sequence_twins(iwin))<limp & abs(Rsq(sequence_models(imodel,iwin),sequence_twins(iwin)))>0
                             txt = sprintf('%0.2f',Rsq(sequence_models(imodel,iwin),sequence_twins(iwin)));
                         elseif abs(Rsq(sequence_models(imodel,iwin),sequence_twins(iwin)))>0
                             txt = sprintf('%0.2f',Rsq(sequence_models(imodel,iwin),sequence_twins(iwin)));
                         else
                             txt = '-';
                         end
                         ht(imodel) = text(xx(imodel+2),yy(iwin+1),txt,'fontsize', fsize-1);
                        
                    end
                end
                        
                iwin = iwin+1;
                ht=[];
                Rsq_sequence_sum = sum(Rsq_sequence,2);
                for imodel = 1:nmodels_sequence
                    txt = sprintf('%0.2f',Rsq_sequence_sum(imodel));
                    ht(imodel) = text(xx(imodel+2),yy(iwin+1),txt,'fontsize', fsize-1);
                end
            end
            axis off

            
            
            %add cellname
            axes('position',[.05 .925 .5 .017])
            text(0,3,['cell: ',cellname,' SO trials ABA pairs'],'fontweight','bold')
            text(0,1.5, 'pvalues of 3-way anova ', ...
                'fontweight','bold','fontsize', fsize)
            axis off
		
            if 		(pairnum==1) yytxt = .655;
            elseif 	(pairnum==2) yytxt = .355;
            elseif 	(pairnum==3) yytxt = .055;
            end
            axes('position',[.05 limy+0.55 .5 .017])
            text(0,15.025, ['R^2 of linear models (without nonzero)'],'fontweight','bold','fontsize', fsize)
            axis off
            
            axes('position',[.05 limy+0.55 .5 .017])
            text(0,-11.025, ['R^2 of linear models(with nonzero)'],'fontweight','bold','fontsize', fsize)
            axis off
            
            set(gcf,'Units','normalized','position',[0.5 0 0.5 1])
            set(gcf,'PaperPositionMode','auto')
		
		
		
		
		
		
		
		
		
		
		
		
		
        elseif 0 % else pairs in SO

            hf(4)=figure(4); %hold on
            if ipair==1
                starty=0.95; limy=0.6;
            elseif ipair==2
                starty=0.45; limy=0.1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %retrieve Rsq, nonzero and modellabs
%           eval(['pp1 = anovastats.',pairname,'.byOffertypeOrderSideChoosen.pval;']); 
%           eval(['pp2 = anovastats.',pairname,'.bytrialtype.pval;']);
%           eval(['pp3 = anovastats.',pairname,'.byOfOrChinteractions.pval;']);
%           pp=horzcat(pp1,pp2,pp3); pp=nanmin(pp'); pp=pp';

            eval(['pp = anovastats.',JCorSO,'.',pairname,'.bytrialtype.pval;']);

		
            eval(['Rsq = ',fittoplot,'.',JCorSO,'.',pairname,'.',Rsqtouse,';'])
            eval(['nonzero = ',fittoplot,'.',JCorSO,'.',pairname,'.nonzero.',nonzero_confint_touse,';'])
            eval(['modellabs = ',fittoplot,'.',JCorSO,'.',pairname,'.modellabs;'])
            Rsq(isnan(Rsq)) = 0;
            nonzero(isnan(nonzero)) = 0;
		
%           %
%           %select and reorder models
%           nmodels = size(Rsq,1);
%           if nmodels>4
%               ind = reorder_models;
%               Rsq = Rsq(ind,:);
%               nonzero = nonzero(ind,:);
%               modellabs = modellabs(ind);
%               nmodels = size(Rsq,1);
%           end
		
            %
            %plot tuningfit
            stp=1.2/(nmodels+2);
            xx = [-.1:stp:1.1];
            stp=(starty-limy)/(size(rnames,1)+1);
            yy = [starty:-stp:limy];
            ht = [];
            for iwin = 1:ntwins
                ht(iwin) = text(xx(1),	yy(iwin+1),	timewindows{iwin}(1),	'fontsize', fsize-1);
            end
            for iwin = 1:ntwins
                if (pp(iwin)<.01)		set(ht(iwin),'color',[.7 0 0]);	end
                if (pp(iwin)<.001)	set(ht(iwin),'color',[1 0 0]);	end
            end
		
            for imodel = 1:nmodels
                text(xx(imodel+2),yy(1),modellabs{imodel},'fontsize', fsize-2,'rotation',45)
            end
		
            for iwin = 1:ntwins
                ht = [];
                for imodel = 1:nmodels
				%show only variables that explain the response
                    % if pp(iwin)<limp && nonzero(imodel,iwin)
                    if pp(iwin)<limp % && nonzero(imodel,iwin)
                        txt = sprintf('%0.2f',Rsq(imodel,iwin));
                    else
                        if Rsq(imodel,iwin)>0
                            txt = sprintf('%0.2f',Rsq(imodel,iwin));
                        else
                            txt = '-';
                        end
                        % txt = sprintf('%0.2f',Rsq(imodel,iwin));
                    end
                    ht(imodel) = text(xx(imodel+2),yy(iwin+1),txt,'fontsize', fsize-1);
                end
                %add color codes
                [Rmax,bestmod] = max(Rsq(:,iwin));
                if nonzero(bestmod,iwin)
                    if		(pp(iwin)<.001)	set(ht(bestmod),'color',[1 0 0]);
                    elseif	(pp(iwin)<.01)	set(ht(bestmod),'color',[.7 0 0]);
                    else						set(ht(bestmod),'color',[0 0 .7]);
                    end
                end
            end
		

            if ipair == 1
                %axis off
                %add cellname
                %axes('position',[.05 1 .5 .017])
                text(-.15,starty+0.1,[ 'cell: ',cellname, ' SO trials ',pairname,' pairs'],'fontweight','bold')
                text(-.15,starty+0.075, 'R^2 of linear models','fontweight','bold','fontsize', fsize)
                axis off
		
		
                set(gcf,'Units','normalized','position',[0.5 0 0.5 1])
                set(gcf,'PaperPositionMode','auto')
                axis off
            elseif ipair == 2 
                text(-.15,starty+0.1,[ 'cell: ',cellname, ' SO trials ',pairname,' pairs'],'fontweight','bold')
                text(-.15,starty+0.075, 'R^2 of linear models','fontweight','bold','fontsize', fsize)
                axis off
            end
				
        end
		
    end
end

nplots = length(hf);
% nplots = 3;