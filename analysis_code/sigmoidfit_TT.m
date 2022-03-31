function [pairlist, psyphycell, mdl] = sigmoidfit_TT(cellname,fitt, varargin)
%
% Fits the choice patterns with a sigmoid normcdf(x,x0,w), and returns
% psyphycell, a matrix that contains one row for each juice pair:
%		psyphycell(pair,:) = [fittedmodel.x0, fittedmodel.w, Rsq]
% The fit is done in the space of number ratios nB:nA. It can be done in
% either linear space or log space (better).
% For verbose = 1 the function plots the data and fitted sigmoid, as well
% as the probability distribution of relative value (dotted line)



disp('   ... fitting sigmoid (by cell)');

if isempty(varargin) verbose = 0;
else verbose = varargin{1};
end


% fl = 0; leg{1}='AB';leg{2}='BA'; leg{3}='JC'; %; leg{2}='AB'; leg{3}='BA';
% fl=1; %fl+1;

% leg{1}='AB';leg{2}='BA'; leg{3}='JC'; %; leg{2}='AB'; leg{3}='BA';
leg{1}='Task 2 (AB)';leg{2}='Task 2 (BA)'; leg{3}='Task 1'; %; leg{2}='AB'; leg{3}='BA';
htt =  findobj('type','figure');
nplots = length(htt);
fl = nplots + 1;

space = 'log';      % 'log' or 'linear'
dx = .005;			% for Riemann's integral
epsilon = 0.0001;	% for integration domain (how small the distr. before we ignore)

%load data
session = cellname(1:8); readsession_TT
filename = [dirroot, cellname, '_data']; eval(['load ',filename])
filename = [dirroot, cellname, '_bhvParams']; eval(['load ',filename])
npairs_JC = size(pairs_JC,1);
npairs_SO = size(pairs_SO,1);
if npairs_JC~=npairs_SO
    disp('error: pair numbers are different between JC and SO');
else
    npairs = npairs_JC;
end
if pairlist_JC~=pairlist_SO
    disp('error: pair lists are different between JC and SO');
else
    pairlist = pairlist_JC;
end

%probability ratios
prob_ratios = zeros(npairs,1);
for ipair = 1:npairs
	prob_ratios(ipair) = sessionParams.goods(pairlist(ipair,1)).probability / sessionParams.goods(pairlist(ipair,2)).probability;
end

%alloffers_num
%alloffers_num_JC
alloffers_num_JC = [];
for ipair = 1:npairs
	table01_JC = table01_all_JC.pooldir{ipair};
	alloffers_num_JC = [alloffers_num_JC; table01_JC(:,1:2)];
end
alloffers_num_JC = unique(abs(alloffers_num_JC),'rows');
eps = 0.001;
aux = alloffers_num_JC + eps;
%keyboard
[~, jnd] = sort(aux(:,2)./aux(:,1));
alloffers_num_JC = alloffers_num_JC(jnd,:);
%
%alloffers_num_SO
alloffers_num_SO = [];
for ipair = 1:npairs
	table01_SO = table01_all_SO.pooldir{ipair};
	alloffers_num_SO = [alloffers_num_SO; table01_SO(:,1:2)];
end
alloffers_num_SO = unique(abs(alloffers_num_SO),'rows');
eps = 0.001;
aux = alloffers_num_SO + eps;
%keyboard
[~, jnd] = sort(aux(:,2)./aux(:,1));
alloffers_num_SO = alloffers_num_SO(jnd,:);
%
if alloffers_num_JC~=alloffers_num_SO
    disp('error: all offers are different between JC and SO');
else
    alloffers_num = alloffers_num_JC;
end



for ipair = 1:npairs
	for dirswitch=1
		
        %
        % SO
        %
		table01_SO = [ abs((table01_all_SO.bydir{ipair,1})); abs((table01_all_SO.bydir{ipair,2}));];
		colrs = [0 0 0]; linew=3; dsize=120;
		
		%separate non-forced choice trials and forced choices
		table1mod_SO = table01_SO(logical(table01_SO(:,1) & table01_SO(:,2)),:);
		forcedAtab_SO = table01_SO(logical(table01_SO(:,1) & ~table01_SO(:,2)),:);
		forcedBtab_SO = table01_SO(logical(~table01_SO(:,1) & table01_SO(:,2)),:);
		nfA = size(forcedAtab_SO,1);
		nfB = size(forcedBtab_SO,1);
		%
		xx_SO = table1mod_SO(:,2)./table1mod_SO(:,1);
		if		isequal(space,'linear')
		elseif	isequal(space,'log')
			xx_SO = log(xx_SO);
		end
		yy_SO = table1mod_SO(:,3);
		yy2_SO = table1mod_SO(:,3).*table1mod_SO(:,4);
		
		%ORDER
		for n=1:numel(yy_SO)
			if n<=(numel(yy_SO))/2; ordere='AB'; else ordere='BA'; end
			xx3_SO{n}=ordere;
		end
		xx3a_SO=strcmp(xx3_SO,'AB');
		xx3b_SO=strcmp(xx3_SO,'BA');
		xx3_SO=xx3a_SO-xx3b_SO;
        %xx3=['AB';'AB';'AB';'AB';'AB';'AB';'AB';'BA';'BA';'BA';'BA';'BA';'BA';'BA'];
		
        %
        % JC
        %
		table01_JC = [ abs((table01_all_JC.pooldir{ipair,1})) ];
		colrs = [0 0 0]; linew=3; dsize=140;
		
		%separate non-forced choice trials and forced choices
		table1mod_JC = table01_JC(logical(table01_JC(:,1) & table01_JC(:,2)),:);
		forcedAtab_JC = table01_JC(logical(table01_JC(:,1) & ~table01_JC(:,2)),:);
		forcedBtab_JC = table01_JC(logical(~table01_JC(:,1) & table01_JC(:,2)),:);
		nfA = size(forcedAtab_JC,1);
		nfB = size(forcedBtab_JC,1);
		%
		xx_JC = table1mod_JC(:,2)./table1mod_JC(:,1);
		if		isequal(space,'linear')
		elseif	isequal(space,'log')
			xx_JC = log(xx_JC);
		end
		yy_JC = table1mod_JC(:,3);
		yy2_JC = table1mod_JC(:,3).*table1mod_JC(:,4);
		
        
        %keyboard
		% CHOOSING THE FITTING FUNCTION
		if strcmp(fitt,'normcdf')
            % NOT USING ANYMORE!! -WS 09/2019
            disp('fitt is normcdf which is not used anymore!');
            keyboard
        else
			%keyboard
			% colrs = colrs+0.15;
            %
            % SO
            %
			Binosize_SO=table1mod_SO(:,4);
			%keyboard
			try				
				tbl=table(xx_SO,xx3_SO',yy2_SO,'VariableNames',{'offer','order','choice'});
 				[tbl idx]=sortrows(tbl,'order','descend');
				Binosize_SO=Binosize_SO(idx);
				mdl_SO=fitglm(tbl, 'choice ~ offer + order ','Distribution','binomial','link',fitt,'BinomialSize',Binosize_SO); % tell glmfit to use the binomial response
			catch
				keyboard
            end
			coefit_SO = [mdl_SO.Coefficients.Estimate(1), mdl_SO.Coefficients.Estimate(2) mdl_SO.Coefficients.Estimate(3)];
			
            %
            % JC
            %
            Binosize_JC=table1mod_JC(:,4);
			%keyboard
			try				
				tbl=table(xx_JC,yy2_JC,'VariableNames',{'offer','choice'});
				mdl_JC=fitglm(tbl, 'choice ~ offer ','Distribution','binomial','link',fitt,'BinomialSize',Binosize_JC); % tell glmfit to use the binomial response
			catch
				keyboard
            end
			coefit_JC = [mdl_JC.Coefficients.Estimate(1), mdl_JC.Coefficients.Estimate(2)];
			
			%
            % BOTH (assume the same width)
            %
            Binosize_both = [Binosize_SO;Binosize_JC];
            try
                xx_both = [xx_SO;xx_JC];
                xx2_both = [-ones(size(xx_SO));ones(size(xx_JC))];
                xx3_both = [xx3_SO';zeros(size(xx_JC))];
                yy2_both = [yy2_SO; yy2_JC];
                tbl=table(xx_both,xx2_both,xx3_both,yy2_both,'VariableNames',{'offer','task','order','choice'});
                [tbl idx]=sortrows(tbl,'task','descend');
                Binosize_both=Binosize_both(idx);
				mdl_both=fitglm(tbl, 'choice ~ offer + task + order ','Distribution','binomial','link',fitt,'BinomialSize',Binosize_both); % tell glmfit to use the binomial response
			catch
				keyboard
            end
            coefit_both = [mdl_both.Coefficients.Estimate(1), mdl_both.Coefficients.Estimate(2) mdl_both.Coefficients.Estimate(3),mdl_both.Coefficients.Estimate(4)];
			           
            
			%find the end of integration domain, imposing value distr. <
			%epsilon %%%% ????? to change with logit????
			leftlim = min(xx_SO);
			while normpdf(leftlim,coefit_SO(1),coefit_SO(2))>epsilon
				leftlim = leftlim-dx;
			end
			rightlim = max(xx_SO);
			while normpdf(rightlim,coefit_SO(1),coefit_SO(2))>epsilon
				rightlim = rightlim+dx;
			end

			%compute relvalue (NB: in log space, it is not just the center of the value distrib.)
			if	isequal(space,'linear')
                %
                % SO
				relvalue_SO = coefit_SO(1);
				width_SO = coefit_SO(2);
                %
                % JC
				relvalue_JC = coefit_JC(1);
				width_JC = coefit_JC(2);
                %
                % both
                relvalue_both = coefit_both(1);
				width_both = coefit_both(2);
                
			elseif	isequal(space,'log')
				%keyboard
                %
                % SO
				relvalue_SO = exp(-1*coefit_SO(1)/coefit_SO(2));
				width_SO = exp(1/coefit_SO(2));
                orderbias_SO = 2*relvalue_SO*coefit_SO(3)/coefit_SO(2);
                
				delta1=1+(coefit_SO(3)/coefit_SO(2));
				delta2=1-(coefit_SO(3)/coefit_SO(2));
				
				relvalueord(1) = relvalue_SO*delta2;
				relvalueord(2) = relvalue_SO*delta1;
                
                %
                % JC
                relvalue_JC = exp(-1*coefit_JC(1)/coefit_JC(2));
				width_JC = exp(1/coefit_JC(2));
                
                %
                % both
                width_both = exp(1/coefit_both(2));
                
                relvalue_both_JC = exp(-1*(coefit_both(1)+coefit_both(3))/coefit_both(2));
                relvalue_both_SO = exp(-1*(coefit_both(1)-coefit_both(3))/coefit_both(2));
                
                delta1_both = 1+(coefit_both(4)/coefit_both(2));
				delta2_both = 1-(coefit_both(4)/coefit_both(2));
                relvalueord_both(1) = relvalue_both_SO*delta2_both;
				relvalueord_both(2) = relvalue_both_SO*delta1_both;

            end
            %
            % SO
			Rsq_SO = mdl_SO.Rsquared.Adjusted;
			% JC
			Rsq_JC = mdl_JC.Rsquared.Adjusted;
			% both
            Rsq_both = mdl_both.Rsquared.Adjusted;
            
			%Distribution	Link Function Name	Link Function	       Mean (Inverse) Function
			%'binomial'	    'logit'				f(?) = log(?/(1–?))	   ? = exp(Xb) / (1 + exp(Xb))
			%x = [min(xx)-log(2)/2 : .05 : max(xx)+log(2)/2];
% 			ord2={'AB';'BA'}; 
            %
            % SO
			ord2=[1;-1]; x=[-2:0.025:3];
			for ord=1:2
				for n=1:numel(x)
					order2(n)=ord2(ord);
				end
				tbl_pred{ord}=table(x',order2','VariableNames',{'offer','order'});

				[y_choicefit{ord}, ystd{ord}] = predict(mdl_SO,tbl_pred{ord});	
            end
			%
            % JC
            tbl_pred{3}=table(x','VariableNames',{'offer'});
			[y_choicefit{3}, ystd{3}] = predict(mdl_JC,tbl_pred{3});
                
            
			psyphycell.SO{1} = mdl_SO;
			psyphycell.SO{2}(dirswitch,:) = coefit_SO;
			psyphycell.SO{3}(dirswitch,:) = [relvalue_SO, width_SO, Rsq_SO, relvalueord,delta1,delta2];	
            
            psyphycell.JC{1} = mdl_JC;
			psyphycell.JC{2}(dirswitch,:) = coefit_JC;
			psyphycell.JC{3}(dirswitch,:) = [relvalue_JC, width_JC, Rsq_JC];	
            
            mdl={mdl_SO,mdl_JC};
            
            psyphycell.both{1} = mdl_both;
            psyphycell.both{2}(dirswitch,:) = coefit_both;
            psyphycell.both{3}(dirswitch,:) = [relvalue_both_JC,relvalue_both_SO,relvalueord_both,width_both];	
            
		end
		
		
		%take care of outliers
		if (isequal('space','log') && confinterval(2)>5) || (isequal('space','linear') && confinterval(2)>100)
			confinterval(2) = NaN;
			disp(['manually corrected confinterval']);
		end
		
				
		%this should be made part of the function plot_sigmoid
		if verbose
			%initialize figures
			if ipair==1 && dirswitch==1
				hf = figure(fl); set(gcf, 'position',[650 620 425 350], 'PaperPositionMode','auto');
				if dospikecheck
				axes('position',[.1 .57 .8 .4]); 
				else
				axes('position',[.1 .37 .8 .6]); 
				end
				hold on; box on; set(gca,'fontsize',8);
            end
			
			%
			%plot fitted function
			%keyboard
			
			ptsord={'v';'^';'o'};
 			colord={[1. .2 .2];[ .2 .2 1.];[.5 .5 .5]}; lord=(numel(yy_SO)/2); % or numel(yy_JC)
			
			% for ord=1:3 % AB BA JC
%             for ord=[3 1 2] % JC AB BA
            for ord=[1 2] % AB BA
				if ord < 3 % AB BA
                    intval=0;
                    p{ord,fl}=plot(x,y_choicefit{ord},'-','color',colord{ord},'linewidth',linew,'DisplayName',leg{ord});
                    hht(ord) = p{ord,fl};
                    % plot(xx_SO(lord*(ord-1)+1:lord*(ord)),yy_SO(lord*(ord-1)+1:lord*(ord)),ptsord{ord},'color',colord{ord},'markersize',dsize,'markerfacecolor',colord{ord},'HandleVisibility','off')
                    sscat1 = scatter(xx_SO(lord*(ord-1)+1:lord*(ord)),yy_SO(lord*(ord-1)+1:lord*(ord)),dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord},'DisplayName','off');
                    sscat1.MarkerFaceAlpha = .4;
                    %plot choice pattern + forced choices
                    % xx_forcedA = min(xx_SO)-log(2)*[1:nfA];	xx_forcedA = sort(xx_forcedA)';
                    % xx_forcedB = max(xx_SO)+log(2)*[1:nfB];	xx_forcedB = sort(xx_forcedB)';
                    xx_forcedA = min(xx_SO)-log(1.3)*[1:nfA];	xx_forcedA = sort(xx_forcedA)';
                    xx_forcedB = max(xx_SO)+log(1.3)*[1:nfB];	xx_forcedB = sort(xx_forcedB)';
                    % plot(xx_forcedA,zeros(1,nfA),ptsord{ord},'color',colord{ord},'markersize',dsize,'markerfacecolor',colord{ord},'HandleVisibility','off')
                    % plot(xx_forcedB, ones(1,nfB),ptsord{ord},'color',colord{ord},'markersize',dsize,'markerfacecolor',colord{ord},'HandleVisibility','off')
                    sscat2 = scatter(xx_forcedA,zeros(1,nfA),dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord});
                    sscat2.MarkerFaceAlpha = .4;
                    sscat3 = scatter(xx_forcedB, ones(1,nfB),dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord},'DisplayName','off');
                    sscat3.MarkerFaceAlpha = .4;    
                else
                    intval=0;
                    p{ord,fl}=plot(x,y_choicefit{ord},'-','color',colord{ord},'linewidth',linew,'DisplayName',leg{ord});
                    hht(ord) = p{ord,fl};
                    % plot(xx_JC,yy_JC,ptsord{ord},'color',colord{ord},'markersize',dsize,'markerfacecolor',colord{ord},'HandleVisibility','off');
                    sscat1 = scatter(xx_JC,yy_JC,dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord});
                    sscat1.MarkerFaceAlpha = .4;
                    %plot choice pattern + forced choices
                    % xx_forcedA = min(xx_JC)-log(2)*[1:nfA];	xx_forcedA = sort(xx_forcedA)';
                    % xx_forcedB = max(xx_JC)+log(2)*[1:nfB];	xx_forcedB = sort(xx_forcedB)';
                    xx_forcedA = min(xx_JC)-log(1.3)*[1:nfA];	xx_forcedA = sort(xx_forcedA)';
                    xx_forcedB = max(xx_JC)+log(1.3)*[1:nfB];	xx_forcedB = sort(xx_forcedB)';
                    % plot(xx_forcedA,zeros(1,nfA),ptsord{ord},'color',colord{ord},'markersize',dsize,'markerfacecolor',colord{ord},'HandleVisibility','off')
                    % plot(xx_forcedB, ones(1,nfB),ptsord{ord},'color',colord{ord},'markersize',dsize,'markerfacecolor',colord{ord},'HandleVisibility','off')
                    sscat2 = scatter(xx_forcedA,zeros(1,nfA),dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord});
                    sscat2.MarkerFaceAlpha = .4;
                    sscat3 = scatter(xx_forcedB, ones(1,nfB),dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord});
                    sscat3.MarkerFaceAlpha = .4;                   
                end
            end
			legend 'show';
            % legend(hht([3,1,2]));
            legend(hht([1,2]));
			%
			%cosmetics
			set(gca,'xlim',[min(xx_SO)-1.1*log(1.3) max(xx_SO)+1.1*log(1.3)])
            box off; axis square
			[~,ind,~] = unique(xx_SO);	%remove doubles in xx
			xxx = xx_SO(ind);
			
% 			text((min(x)+dirswitch*.66),.925,['\rho = ',	sprintf('%1.3g',psyphycell{3}(1,1)+intval)],'fontsize',9,'Color',colrs(ipair,:))
% 			text((min(x)+dirswitch*.66),.855,['\mu = ',sprintf('%1.3g',psyphycell{3}(1,2))],	'fontsize',9,'Color',colrs(ipair,:))
% 			text((min(x)+dirswitch*.66),.725,['Rsq = ',sprintf('%1.3g',psyphycell{3}(1,3))],	'fontsize',9,'Color',colrs(ipair,:))
% 			text((min(x)+dirswitch*.66),.655,['\delta AB = ',sprintf('%1.3g',relvalueord(1)), '  \delta BA = ',sprintf('%1.3g',relvalueord(2))],	'fontsize',9,'Color',colrs(ipair,:))
% 			text((min(x)+dirswitch*.66),.585,['p.' '\delta = ',sprintf('%1.3g',mdl.Coefficients.pValue(3))],	'fontsize',9,'Color',colrs(ipair,:))
			

			xpos=min(xx_SO);

% 			text(xpos,.925,['\rho(JC) = ',	sprintf('%1.3g',psyphycell.JC{3}(1,1)+intval)],'fontsize',9,'Color',colrs(ipair,:))
% 			text(xpos,.855,['\mu(JC) = ',sprintf('%1.3g',psyphycell.JC{3}(1,2))],	'fontsize',9,'Color',colrs(ipair,:))
% 			% text(xpos,.785,['Rsq(JC) = ',sprintf('%1.3g',psyphycell.JC{3}(1,3))],	'fontsize',9,'Color',colrs(ipair,:))
			
			if mdl_SO.Coefficients.pValue(3)<=0.05
			text(xpos,.655,['\rho(SO) = ',sprintf('%1.3g',relvalue_SO), ' \rho BA = ',sprintf('%1.3g',relvalueord(2)), '  \rho AB = ',sprintf('%1.3g',relvalueord(1))],	'fontsize',9,'Color',colrs(ipair,:))
            text(xpos,.585,['\mu(SO) = ',sprintf('%1.3g',psyphycell.SO{3}(1,2))],	'fontsize',9,'Color',colrs(ipair,:))
			text(xpos,.515,['\epsilon(SO) = ',sprintf('%1.3g',orderbias_SO)],	'fontsize',9,'Color',colrs(ipair,:))
            % text(xpos,.585,['p.' '\rho = ',sprintf('%1.3g',mdl_SO.Coefficients.pValue(3))],	'fontsize',9,'Color','b')
			else
			text(xpos,.655,['\rho(SO) = ',sprintf('%1.3g',relvalue_SO), ' \rho BA = ',sprintf('%1.3g',relvalueord(2)), '  \rho AB = ',sprintf('%1.3g',relvalueord(1))],	'fontsize',9,'Color',colrs(ipair,:))
            text(xpos,.585,['\mu(SO) = ',sprintf('%1.3g',psyphycell.SO{3}(1,2))],	'fontsize',9,'Color',colrs(ipair,:))
            text(xpos,.515,['\epsilon(SO) = ',sprintf('%1.3g',orderbias_SO)],	'fontsize',9,'Color',colrs(ipair,:))
            % text(xpos,.585,['p.' '\rho = ',sprintf('%1.3g',mdl_SO.Coefficients.pValue(3))],	'fontsize',9,'Color',colrs(ipair,:))
			end
% 			text(xpos,.485,['n trials(SO) = ',sprintf('%1.3g',sum(table01_SO(:,4)))],	'fontsize',9,'Color',colrs(ipair,:))
% 			text(xpos,.415,['n trials(JC) = ',sprintf('%1.3g',sum(table01_JC(:,4)))],	'fontsize',9,'Color',colrs(ipair,:))
			
			if dirswitch==1;
				%
				xlab = [];		%xlabels
				for ifA = 1:nfA
					xlab{ifA} = [num2str(forcedAtab_SO(nfA-ifA+1,2)),':',num2str(forcedAtab_SO(nfA-ifA+1,1))];
				end
				for i = 1:size(xxx,1)
					xlab{nfA+i} = [num2str(table1mod_SO(ind(i),2)),':',num2str(table1mod_SO(ind(i),1))];
				end
				for ifB = 1:nfB
					xlab{nfA+i+ifB} = [num2str(forcedBtab_SO(nfB-ifB+1,2)),':',num2str(forcedBtab_SO(nfB-ifB+1,1))];
				end
				%
				%add forced choices
				xxx = [xx_forcedA;xxx;xx_forcedB];
				set(gca,'xtick',xxx,'xticklabel',xlab)
				set(gca,'ylim',[0,1]);
				set(gca,'ytick',[0:.25:1],'yticklabel',[0:25:100]);
                xtickangle(45);
                if isequal(space,'log')
                    xlabel('log(qB:qA)');
                elseif isequal(space,'linear')
                    xlabel('qB:qA');
                end
                ylabel('B choice %');
				%
				%mark cellname
				%if ipair == npairs
				
				
                if dospikecheck
                    axes('position',[.1 .57 .8 .4]); 
                else
                    axes('position',[.1 .37 .8 .6]); 
				end
				
				title(strcat(session, ' , Fit:',  fitt, ', (Space: ',space,')'), ...
					'fontname','timesnewroman','fontsize', 8)
				axis off
				set(gcf,'PaperPositionMode','auto')
				
			end
			warning off;
			% 			legend('off');
			% 			legend(gca,'show')
		end %if verbose
	end % for dirswitch
end % for ipairs



if verbose
% keyboard
% CARTON MATRIX
	
if dospikecheck
	axes('position',[.1 .31 .8 .22]);
else
	axes('position',[.1 .02 .8 .28]);
end
				
% axes('position',[.1 .02 .8 .28]);
hold on; box on; set(gca,'fontsize',8)
% axis equal
axis([0 10 0 5])
axis square
set(gca,'xtick',0:10,'ytick',0:5)
grid on
colord={[1. .2 .2];[ .2 .2 1.]};
for ord = 1:2
	hold on
	plot([0 10*relvalueord(ord)],[0 10],'--','linewidth',2,'Color',colord{ord})
% 	title([strcat(session, ' , Choice Matrix with \rho BA = ',  num2str(relvalueord(2)), ' and \rho AB = ',  num2str(relvalueord(1))  )], 'fontname','timesnewroman','fontsize', 8)
end
plot([0 10*relvalue_SO],[0 10],'--','linewidth',2,'Color','k')
offtypes=table01_SO(:,1:2);
for n=1:size(offtypes,1)
	hold on
	plot(offtypes(n,2),offtypes(n,1),'ko','markerfacecolor','k','markersize',table01_SO(n,4)/2.5)
end
xlabel('nb of B')
ylabel('nb of A')
end


end
