function [hitdata] = get_hitdata_TT(goodTrials, pair, ngoods, JCorSO)
%keyboard
%
% from goodTrials makes the matrix hitdata, defined for a single pair of juices
% 
% for juice choice - columns use to be: trialnum, #A, #B, posA, chosenpos
%
% WRONG!!! for sequential offer - columns are: trialnum, #A, #B, posA(-1=Left),Choosen ID(-1=A), OrderA(-1:First) %% wrong!!!
% correct one: for sequential offer - columns are: trialnum, #A, #B, posA(-1=right), Choosen ID(-1=B), OrderA(-1:Second)
%
% author: camillo, november 2008
% revisions: 12/2009, 12/2010, 12/2016 SB
% revisions: 10/2018  WS for TT

% %make hitdata
if isequal(JCorSO, 'JC')
    goodOffers = goodTrials(:,2:ngoods+1);
    jA = pair(1);
    jB = pair(2);
    indAB = goodOffers(:,jA) & goodOffers(:,jB);
    indA0 = goodOffers(:,jA) & ~sum(goodOffers(:,setdiff(1:ngoods,jA)),2);
    indB0 = goodOffers(:,jB) & ~sum(goodOffers(:,setdiff(1:ngoods,jB)),2);
    ind = indA0 | indB0 | indAB;
    tnoff = [goodTrials(ind,1), abs(goodOffers(ind,[jA jB]))];
    pos_A = sign(goodOffers(ind, jA));
    pos_B = sign(goodOffers(ind, jB));
    pos_A(~pos_A) = -pos_B(~pos_A);

    jnd = goodTrials(ind,end-1)==jA;
    chpos = pos_A .* jnd;
    chpos(~jnd) = -pos_A(~jnd);
    hitdata = [tnoff, pos_A, chpos];

elseif isequal(JCorSO, 'SO')

    goodOffers = abs(goodTrials);
    goodOffers(:,6)=(goodTrials(:,2)>0 | goodTrials(:,3)<0)+1;
    for n=4:6
        goodOffers(goodOffers(:,n)==1,n)= 1; %SECOND
        goodOffers(goodOffers(:,n)==2,n)= -1; %FIRST
    end
    %goodOffers(:,6)=goodOffers(:,6).*goodOffers(:,5);

    hitdata=goodOffers;
    
end

 %FILTER FOR PROFILE & TUNING
 %A: Choice by ChoosenID: 5= -1:choosen A OR 1 Choosen B
 %B: Choice by Order: 5*6= 1 choosen First OR -1 choosen Second
 %C: Choice by Side: 5*4= 1 choosen Left OR -1 choosen Right
 
 
%keyboard
