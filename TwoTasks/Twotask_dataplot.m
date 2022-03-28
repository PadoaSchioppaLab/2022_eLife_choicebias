function [] = Twotask_dataplot(TrialRecord)
%
%

try
    
if ~isfield(TrialRecord.User, 'choicePattern_JC') ||  ~isfield(TrialRecord.User, 'choicePattern_SO')
    return
end

table01_JC = TrialRecord.User.choicePattern_JC;
table01_SO = TrialRecord.User.choicePattern_SO;

if isempty(table01_JC) && isempty(table01_SO)
    return
%     else
%         keyboard
end

cla;
hold on

clrs = reshape([TrialRecord.User.sessionParams.goods.color],3,2)';
clrs(3,:)=mean(clrs);
symbs={'x','+','o'};

temptable=[];
temptable{1}= table01_SO.byorder{:,1};
temptable{2}= table01_SO.byorder{:,2};
temptable{3}= table01_JC.pooldir{:,1};

table01_all=[];
for n=1:3
    for z=1:size(temptable{n},1)
        table01_all=[table01_all; [temptable{n}(z,:) n]];
    end
end

alloffers = (abs(table01_all(:,1:2)));
eps = 0.001;
aux = alloffers + eps;
[~, jnd2] = sort(aux(:,2)./aux(:,1));
for isign=1:3
    xx=[]; yy=[];
    alloffers = abs(table01_all(jnd2,:));
    mask=alloffers(:,5)~=isign;
    alloffers(mask,3)=NaN;
    abs_off = abs(alloffers(:,1:2)); abs_off=unique(abs_off,'rows');
    noffers=size(abs_off,1);
    aux = abs_off + eps;
    [~, jnd3] = sort(aux(:,2)./aux(:,1));
    abs_off=abs_off(jnd3,:);
    for npair=1:noffers
%   [~, ind, jnd] = intersect(abs_off(npair,:), alloffers(:,1:2), 'rows');
  jnd4=find(alloffers(:,1)==abs_off(npair,1) & alloffers(:,2)==abs_off(npair,2));
   xx(npair) = npair;
    yy(npair) = nanmean(alloffers(jnd4,3));
    end
    % colors
    clr = clrs(isign,:);
    symb=symbs{isign};
    % plot
    mask=~isnan(yy);
    plot(xx(mask), yy(mask), '-', 'linewidth',2, 'color',clr); hold on
    plot(xx(mask), yy(mask), symb,'markersize',10,'linewidth',2,'color',clr);
end
    plot([-.25 noffers+.25], [.5 .5], 'k--', 'linewidth',1);

%
%cosmetics
set(gca,'color',[1 1 1])
for ioff = 1:noffers
    xticlab{ioff} = [num2str(abs_off(ioff,2)),':',num2str(abs_off(ioff,1))];
end
set(gca, 'tickdir','out','xtick',[1:noffers])
    axis([.5, noffers+.5, 0, 1]);
    set(gca, 'xticklabel',xticlab, 'fontsize',8,'xticklabelrotation',90)
  set(gca, 'xgrid','on')
set(gca, 'tag', 'rtgraph');
return


catch
    disp('Error in plot function')
%     keyboard
end

end
    
    
    
    
    