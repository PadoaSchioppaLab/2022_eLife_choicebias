%SequentialOffers (timing script)

%This task requires that either an "eye" input or joystick (attached to the
%eye input channels) is available to perform the necessary responses.
%
%During a real experiment, a task such as this should make use of the
%"eventmarker" command to keep track of key actions and state changes (for
%instance, displaying or extinguishing an object, initiating a movement,
%etc).

%init
% hotkey('j','juiceDeliver([], 2, 1, 65);');%juiceline, n of drops
hotkey('j','reward_function(100,''juiceline'',1,''numReward'',1,''pauseTime'',30,''eventmarker'',50);') % FIXSB
putvalue(DAQ.TTL{1}, 0);
% putvalue(DaqInfo.TTL2, 0);

%filename to save TrialRecord
% filename = TrialRecord.User.DataFile;
pathstr=MLConfig.MLPath.ExperimentDirectory; % FIXSB
filename = MLConfig.MLPath.DataFile;  % FIXSB
% [pathstr,name,ext,versn] = fileparts(filename);
% filesave = [pathstr,filesep,name,'.mat'];
filesave = [pathstr,filename,'.mat']; % FIXSB

%write down and pass start_trial
eventmarker(21);			%start trial
putvalue(DAQ.TTL{1}, 1);	%pass start trial to level variable for trace alignement
currTrial	= TrialRecord.User.currTrial;
trialNumber	= currTrial.trialNumber;

%
%load the parameters used in this session
sessionParams	= Twotasks_parameters;
fakemonkey      = sessionParams.fakemonkey;
p_nofixbreak	= sessionParams.p_nofixbreak;
manualchoice	= sessionParams.manualchoice;
probs       	= [sessionParams.goods.probability]';

%give names to the TaskObjects defined in the taskloop:
prefix_point    = 1;
fix_point_SO	= 3; % switch task cue
fix_point_JC	= 2; % switch task cue
fix_point_dim   = 4;
fixpoint_offer1	= 5;
fixpoint_offer2	= 6;
target1         = 7;
target2     	= 8;
offer1      	= 9;
offer2      	= 10;

%define time intervals (in ms):
% SO
wait_for_fix	= sessionParams.trial.wait_for_fix;
initial_prefix	= sessionParams.trial.initial_prefix;
initial_fix_JC	= sessionParams.trial.initial_fix_JC;
initial_fix_SO	= sessionParams.trial.initial_fix_SO;
initial_idle	= sessionParams.trial.initial_idle;
offer1on_time	= sessionParams.trial.offer1on_time;
interoffer_time	= sessionParams.trial.interoffer_time;
offer2on_time	= sessionParams.trial.offer2on_time;
offeroff_time	= sessionParams.trial.offeroff_time;
tgton_range		= sessionParams.trial.tgton_range;
max_reactime	= sessionParams.trial.max_reactime;
max_sacc_time	= sessionParams.trial.max_sacc_time;
sacctgt_idle	= sessionParams.trial.sacctgt_idle;
tgthold_time	= sessionParams.trial.tgthold_time;
preend_time		= sessionParams.trial.preend_time;
% JC
delay_range		= sessionParams.trial.delay_range;
% tgton_range_JC	= sessionParams.trial.tgton_range;

%compute tgton_time for this trial
tgton_time = round((tgton_range(2)-tgton_range(1))*rand(1) + tgton_range(1));
%compute delay_time for this trial
delay_time = round((delay_range(2)-delay_range(1))*rand(1) + delay_range(1));

%fixation window (in degrees):
prefix_radius = sessionParams.prefix_radius;
fix_radius = sessionParams.fix_radius;
% tgt_radius = sessionParams.tgt_radius;
tgt_radiusL = sessionParams.tgt_radiusL;
tgt_radiusR = sessionParams.tgt_radiusR;

% TASK:
chosenOffer = [];
gotjuice = 0;

%%%%% CHECK IF IT IS A SO OR JC TRIALS
TwoT=currTrial.currentOffers(1).Twotasks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JC
if TwoT==1 %%% JC

    %acquire initial fixation:
toggleobject(prefix_point, 'eventmarker', 24);		%prefix point o
if ~fakemonkey
	ontarget = eyejoytrack('acquirefix', prefix_point, prefix_radius, wait_for_fix);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob of no fixation
	if ontarget
		pause(wait_for_fix*.9*rand(1)/1000)
	else
		pause(wait_for_fix/1000)
	end
end

if ~ontarget,
	toggleobject(prefix_point, 'eventmarker', 28);	%fix point off
	rt = NaN;
	trialerror(4);		%no fixation
	eventmarker(60);	%end of trial
	%putvalue(DaqInfo.TTL1, 0);
    putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

%
%hold prefixation 
eyejoytrack('idle', initial_idle);
if ~fakemonkey
	ontarget = eyejoytrack('holdfix', prefix_point, prefix_radius, initial_prefix-initial_idle);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob fix break
	if ontarget
		pause(initial_fix_SO/1000)
	else
		pause(initial_fix_SO*rand(1)/1000)
	end
end
if ~ontarget,
	toggleobject(prefix_point, 'eventmarker', 28);
	rt = NaN;
	trialerror(3);		%broke fixation
	eventmarker(60);	%end of trial
	%putvalue(DaqInfo.TTL1, 0);
    putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

%
%hold fixation and change prefix to usual fixation point
toggleobject([prefix_point fix_point_JC], 'eventmarker', 25);		%prefix off; fix point on; 
eyejoytrack('idle', initial_idle);
if ~fakemonkey
	ontarget = eyejoytrack('holdfix', fix_point_JC, fix_radius, initial_fix_JC-initial_idle);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob fix break
	if ontarget
		pause(initial_fix_SO/1000)
	else
		pause(initial_fix_SO*rand(1)/1000)
	end
end
if ~ontarget,
	toggleobject(fix_point_JC, 'eventmarker', 27);
	rt = NaN;
	trialerror(3);		%broke fixation
	eventmarker(60);	%end of trial
	%putvalue(DaqInfo.TTL1, 0);
    putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

%
%present offers and hold for delay_time
toggleobject([offer1 offer2], 'eventmarker', 30);	%simultaenously displays 2 offers
% putvalue(DaqInfo.TTL2, 1);	%pass offer on to level variable
if ~fakemonkey
	ontarget = eyejoytrack('holdfix', fix_point_JC, fix_radius, delay_time);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob fix break
	if ontarget
		pause(delay_time/1000)
	else
		pause(delay_time*rand(1)/1000)
	end
end
if ~ontarget,
	toggleobject([fix_point_JC offer1 offer2], 'eventmarker', 27, 'eventmarker', 31);
% 	putvalue(DaqInfo.TTL2, 0);
	rt = NaN;
	trialerror(5);		%broke fixation
	eventmarker(60);	%end of trial
	%putvalue(DaqInfo.TTL1, 0);
    putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end
% 
% %deliver juice
% trialerror(0); %correct
% eyejoytrack('idle', tgthold_time);
% toggleobject([fix_point offer1 offer2], 'eventmarker', 27, 'eventmarker', 36);
% for n=1:2
%     if TrialRecord.User.currTrial.currentOffers(n).quantity>0
%         chosenOffer=TrialRecord.User.currTrial.currentOffers(n);
%     end
% end
% 
%  chosenOffer = currTrial.currentOffers(1);
 aux = [TrialRecord.User.currTrial.currentOffers.goodId];
 chosenOffer = TrialRecord.User.currTrial.currentOffers(aux==1);

%
%present saccade targets and check that a saccade is initiated
toggleobject([fix_point_JC fix_point_dim target1 target2], 'eventmarker', 35);	%2 targets; to fix: dim fix point
%toggleobject([fix_point offer1 offer2 fix_point_dim target1 target2], 'eventmarker', 35);	%2 targets; to fix: dim fix point
if fakemonkey && manualchoice
	[scancode rt] = getkeypress(max_reactime);
	if		scancode==71, ontarget = 1;
	elseif	scancode==73, ontarget = 2;
	else		ontarget = 0;
	end
else
	if ~fakemonkey
		[ontarget rt] = eyejoytrack('holdfix', fix_point_JC, fix_radius, max_reactime);
	else
		tnow = trialtime;
		ontarget = rand(1) < .1;		%prob of no response
		if ontarget
			pause(max_reactime/1000)
		else
			pause(max_reactime*rand(1)/1000)
		end
		rt = trialtime - tnow;
	end
	if ontarget,
		toggleobject([fix_point_dim target1 target2 offer1 offer2], 'eventmarker', 27, 'eventmarker', 36);
% 		toggleobject([fix_point_dim target1 target2], 'eventmarker', 27, 'eventmarker', 36);
% 		putvalue(DaqInfo.TTL2, 0);
		rt = NaN;
		trialerror(1);		%no response
		eventmarker(60);	%end of trial
		%putvalue(DaqInfo.TTL1, 0);
        putvalue(DAQ.TTL{1}, 0);
		eval(['save ',filesave,' TrialRecord'])
		return
	end

	%
	%look for choice
	if ~fakemonkey
		ontarget = eyejoytrack('acquirefix', [target1 target2], [tgt_radiusL tgt_radiusR], max_sacc_time);
	else
		ontarget = rand(1) < p_nofixbreak;		%1 - prob late response
		ontarget = ceil(ontarget*2*rand(1));	%randomly select targtet
		if ontarget
			pause(max_sacc_time*rand(1)/1000)
		else
			pause(max_sacc_time/1000)
		end
	end
end

%
if rt < sessionParams.trial.toofast
	toggleobject([fix_point_dim target1 target2 offer1 offer2], 'eventmarker', 27, 'eventmarker', 36);
% 	toggleobject([fix_point_dim target1 target2], 'eventmarker', 27, 'eventmarker', 36);
% 	putvalue(DaqInfo.TTL2, 0);
	trialerror(3);		%too fast
	eventmarker(60);	%end of trial
	%putvalue(DaqInfo.TTL1, 0);
    putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end
%
if ~ontarget,
	toggleobject([fix_point_dim target1 target2 offer1 offer2], 'eventmarker', 27, 'eventmarker', 36);
% 	toggleobject([fix_point_dim target1 target2], 'eventmarker', 27, 'eventmarker', 36);
% 	putvalue(DaqInfo.TTL2, 0);
	trialerror(2);		%too slow
	eventmarker(60);	%end of trial
	%putvalue(DaqInfo.TTL1, 0);
    putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end
%
eventmarker(39);	%choice made
chosenOffer = currTrial.currentOffers(ontarget);

%
%if chosenOffer == null, abort
if ~chosenOffer.quantity
	toggleobject([fix_point_dim target1 target2 offer1 offer2], 'eventmarker', 27, 'eventmarker', 36);
% 	toggleobject([fix_point_dim target1 target2], 'eventmarker', 27, 'eventmarker', 36);
% 	putvalue(DaqInfo.TTL2, 0);
	trialerror(6);		%incorrect response
	eval(['save ',filesave,' TrialRecord'])
	return
end

%
%hold target then reward
if ~fakemonkey
	if		ontarget==1, chtarget = target1; tgt_radius = tgt_radiusL;
	elseif	ontarget==2, chtarget = target2; tgt_radius = tgt_radiusR;
	end
	ontarget = eyejoytrack('holdfix', chtarget, tgt_radius, tgthold_time);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob target fix break
	if ontarget
		pause(tgthold_time/1000)
	else
		pause(tgthold_time*rand(1)/1000)
	end
end
if ~ontarget,
	toggleobject([fix_point_dim target1 target2 offer1 offer2], 'eventmarker', 27, 'eventmarker', 36);
% 	toggleobject([fix_point_dim target1 target2], 'eventmarker', 27, 'eventmarker', 36);
% 	putvalue(DaqInfo.TTL2, 0);
	trialerror(3);		%broke fixation
	eventmarker(60);	%end of trial
	%putvalue(DaqInfo.TTL1, 0);
    putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end
trialerror(0); %correct

%
%turn off all visual stimuli
toggleobject([fix_point_dim target1 target2 offer1 offer2], 'eventmarker', 27, 'eventmarker', 36);
%toggleobject([fix_point_dim target1 target2], 'eventmarker', 27, 'eventmarker', 36);
% putvalue(DaqInfo.TTL2, 0);

%
%deliver juice?
try
jprob = chosenOffer.good.probability;
deliverj = rand(1)<jprob;
psudorand = sessionParams.rand_outcome.psudorand;
catch
    keyboard
end

if jprob<1 && psudorand
	phat0 = sessionParams.rand_outcome.phat0;
	phat1 = sessionParams.rand_outcome.phat1;
	%
	gotjhist_oneg = TrialRecord.User.summaryStats.gotjuice_hist(1:trialNumber-1,chosenOffer.goodId);
	gotjhist_oneg = gotjhist_oneg(~isnan(gotjhist_oneg));
	gotjhist_oneg = [chosenOffer.good.probability*ones(size(gotjhist_oneg)), gotjhist_oneg];
	%
	aux = TrialRecord.User.summaryStats.gotjuice_hist(1:trialNumber-1,:);
	auxp = ~isnan(aux)*probs;	%probabilities
	ind = auxp>0 & auxp<1;		%successful trials with p<1
	auxp = auxp(ind);
	aux = aux(ind,:);
	aux(isnan(aux)) = 0;
	gotjhist_allg = [auxp, sum(aux,2)];
	%
	gotjhist = {gotjhist_oneg, gotjhist_allg};
	for icrit = 1:2		%all-goods criterion applies last
		gjh = [gotjhist{icrit}; jprob, deliverj];	%add current trial
		if size(gjh)
			n0 = size(gjh,1)-find(gjh(:,2)==1,1,'last');
			n1 = size(gjh,1)-find(gjh(:,2)==0,1,'last');
			if n1
				P1 = prod(gjh(end-n1+1:end,1));
				if P1<phat1
					deliverj = 0;
					TrialRecord.User.summaryStats.psudorand(trialNumber) =-1;
					if icrit==1
						disp('NB: applied pseudorand --  one good, forced p=0')
					elseif icrit==2
						disp('NB: applied pseudorand --  across goods, forced p=0')
					end
				end
			elseif n0		%poor luck criterion applies last
				P0 = prod(1-gjh(end-n0+1:end,1));
				if P0<phat0
					deliverj = 1;
					TrialRecord.User.summaryStats.psudorand(trialNumber) = 1;
					if icrit==1
						disp('NB: applied pseudorand --  one good, forced p=1')
					elseif icrit==2
						disp('NB: applied pseudorand --  across goods, forced p=1')
					end
				end
			end %if n0, n1
		end %if size(gjh)
	end %for icrit
end %if pseudorand

%
tic;		%start preend_time before juice delivery
if ~deliverj
	toggleobject(poorluck_sound, 'eventmarker', 49);
else
	juiceline	= chosenOffer.good.juiceline;
	ndrops	= chosenOffer.quantity;
	quantum	= chosenOffer.good.quantum;
    
 juiceDeliver(DAQ, juiceline, ndrops, quantum); % FIXSB
       
	gotjuice	= 1;
end
%
%finish preend_time
jdt = round(1000*toc);	%juice delivery time in ms
if	jdt > preend_time
	disp(['exceeded preend_time: juice delivery took ',num2str(jdt),' ms!'])
else
% 	idle(preend_time - jdt);
	ontarget = eyejoytrack('idle', preend_time - jdt);
end	
eventmarker(60);	%end of trial
%putvalue(DaqInfo.TTL1, 0);
putvalue(DAQ.TTL{1}, 0);
%update records
trialIndex = currTrial.trialIndex;
nTrialsHalfBlock = size(TrialRecord.User.trialsDoneThisHalfBlock,1);
TrialRecord.User.trialsDoneThisHalfBlock(trialIndex) = 1;
if sum(TrialRecord.User.trialsDoneThisHalfBlock) == nTrialsHalfBlock
	TrialRecord.User.halfBlocksDone = TrialRecord.User.halfBlocksDone + 1;
end
if TrialRecord.User.halfBlocksDone == 2
	TrialRecord.User.blocksDone = TrialRecord.User.blocksDone + 1;
end
%
currTrial.chosenOffer	= chosenOffer;
currTrial.gotjuice		= gotjuice;
TrialRecord.User.currTrial	= currTrial;
TrialRecord.User.trialRecord(trialNumber) = currTrial;
%
%update goodTrials and compute choicePattern
ngoods = TrialRecord.User.ngoods;
currentOffers = currTrial.currentOffers;
thistrial = zeros(1, ngoods+3);
thistrial(1) = trialNumber;
for ioff = 1:2
	goodId = currentOffers(ioff).goodId;
	if goodId
		thistrial(1+goodId) = currentOffers(ioff).quantity * sign(currentOffers(ioff).xyposition(1));
	end
end
thistrial(ngoods+2) = chosenOffer.goodId;
thistrial(ngoods+3) = gotjuice;
TrialRecord.User.goodTrials_JC = [TrialRecord.User.goodTrials_JC; thistrial];
%
%summary stats
ind = (sign(chosenOffer.xyposition(1)) + 3)/2;
TrialRecord.User.summaryStats.LRchoice(ind) = TrialRecord.User.summaryStats.LRchoice(ind) + 1;
TrialRecord.User.summaryStats.gotjuice(:, chosenOffer.goodId) = TrialRecord.User.summaryStats.gotjuice(:, chosenOffer.goodId) + [1; gotjuice];
TrialRecord.User.summaryStats.gotjuice_hist(trialNumber, chosenOffer.goodId) = gotjuice;

if gotjuice,
   xx = TrialRecord.User.currTrial.chosenOffer.goodId;
   TrialRecord.User.summaryStats.quantity(1,xx) = TrialRecord.User.summaryStats.quantity(1,xx) + TrialRecord.User.currTrial.chosenOffer.quantity;
end	
%
%choice pattern
[pairlist, table01_all_JC] = make_choicepattern(TrialRecord.User.goodTrials_JC);
% probs = [sessionParams.goods.probability]';
% TrialRecord.User.problist = probs(pairlist); % FIXSB
TrialRecord.User.choicePattern_JC = table01_all_JC;

%save TrialRecord to file
eval(['save ',filesave,' TrialRecord'])
    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SO
elseif TwoT==2 %%%% SO
 
    
%acquire initial fixation:
toggleobject(prefix_point, 'eventmarker', 24);		%prefix point o
if ~fakemonkey
	ontarget = eyejoytrack('acquirefix', prefix_point, prefix_radius, wait_for_fix);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob of no fixation
	if ontarget
		pause(wait_for_fix*.9*rand(1)/1000)
	else
		pause(wait_for_fix/1000)
	end
end

if ~ontarget,
	toggleobject(prefix_point, 'eventmarker', 28);	%fix point off
	rt = NaN;
	trialerror(4);		%no fixation
	eventmarker(60);	%end of trial
	%putvalue(DaqInfo.TTL1, 0);
    putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

%
%hold prefixation 
eyejoytrack('idle', initial_idle);
if ~fakemonkey
	ontarget = eyejoytrack('holdfix', prefix_point, prefix_radius, initial_prefix-initial_idle);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob fix break
	if ontarget
		pause(initial_fix_SO/1000)
	else
		pause(initial_fix_SO*rand(1)/1000)
	end
end
if ~ontarget,
	toggleobject(prefix_point, 'eventmarker', 28);
	rt = NaN;
	trialerror(3);		%broke fixation
	eventmarker(60);	%end of trial
	%putvalue(DaqInfo.TTL1, 0);
    putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

% %
% %acquire initial fixation:
% toggleobject(fix_point, 'eventmarker', 25);		%fix point on
% if ~fakemonkey
% 	ontarget = eyejoytrack('acquirefix', fix_point, fix_radius, wait_for_fix);
% else
% 	ontarget = rand(1) < p_nofixbreak;	%1 - prob of no fixation
% 	if ontarget
% 		pause(wait_for_fix*rand(1)/1000)
% 	else
% 		pause(wait_for_fix/1000)
% 	end
% end
% if ~ontarget,
% 	toggleobject(fix_point, 'eventmarker', 27);	%fix point off
% 	rt = NaN;
% 	trialerror(4);		%no fixation
% 	eventmarker(60);	%end of trial
% 	putvalue(DAQ.TTL{1}, 0);
% 	eval(['save ',filesave,' TrialRecord'])
% 	return
% end

%if rand(1) > 0.9,
%juiceDeliver(DaqInfo, 2, 1, 60);
%end

%hold fixation and change prefix to usual fixation point
toggleobject([prefix_point fix_point_SO], 'eventmarker', 25);		%prefix off; fix point on; 
eyejoytrack('idle', initial_idle);
if ~fakemonkey
	ontarget = eyejoytrack('holdfix', fix_point_SO, fix_radius, initial_fix_SO-initial_idle);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob fix break
	if ontarget
		pause(initial_fix_SO/1000)
	else
		pause(initial_fix_SO*rand(1)/1000)
	end
end
if ~ontarget,
	toggleobject(fix_point_SO, 'eventmarker', 27);
	rt = NaN;
	trialerror(3);		%broke fixation
	eventmarker(60);	%end of trial
	putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

%
%present offer1 and hold for offer1on_time
toggleobject([fixpoint_offer1 offer1] , 'eventmarker', 30);	%displays first offer and change fp color

%STIM ON
% putvalue(DAQ.TTL{1}, 1);


% putvalue(DaqInfo.TTL2, 1);	%pass offer on to level variable
if ~fakemonkey
	ontarget = eyejoytrack('holdfix', fix_point_SO, fix_radius, offer1on_time);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob fix break
	if ontarget
		pause(offer1on_time/1000)
	else
		pause(offer1on_time*rand(1)/1000)
	end
end
if ~ontarget,
	toggleobject([fix_point_SO offer1 fixpoint_offer1], 'eventmarker', 27, 'eventmarker', 31);
	% 	putvalue(DaqInfo.TTL2, 0);
	rt = NaN;
	trialerror(3);		%broke fixation
	eventmarker(60);	%end of trial
	putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	
%STIM OFF
% putvalue(DAQ.TTL{1}, 0);
	
	return
end


%STIM OFF
% putvalue(DAQ.TTL{1}, 0);

%
%remove offer1 and hold for interoffer_time
toggleobject([fixpoint_offer1 offer1] , 'eventmarker', 31); %And start stimulation	

% putvalue(DaqInfo.TTL2, 0);	%pass offer on to level variable
t1=clock;
if ~fakemonkey	
	ontarget = eyejoytrack('holdfix', fix_point_SO, fix_radius, interoffer_time);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob fix break
	if ontarget
		pause(interoffer_time/1000)
	else
		pause(interoffer_time*rand(1)/1000)
	end
end
% if ~ontarget
% disp('blink during interoffer')
% 	toggleobject(fix_point_SO , 'eventmarker', 27);
% 	rt = NaN;
% 	trialerror(3);		%broke fixation
% 	eventmarker(60);	%end of trial
% 	eventmarker(32);	%stop stimulation
% 	putvalue(DAQ.TTL{1}, 0);
% 	eval(['save ',filesave,' TrialRecord'])
% 	return
% % t2=clock;
% % t3=etime(t2,t1)*1000;
% % if t3>5
% % eyejoytrack('idle', interoffer_time-t3);
% % end
% % eventmarker(55);	%break fixation allowed
% end

%eyejoytrack('idle', interoffer_time); % Instead of fixation

%
%present offer2 and hold for offer2on_time
toggleobject([fixpoint_offer2 offer2] , 'eventmarker', 32);	%displays second offer change fp color

% putvalue(DaqInfo.TTL2, 1);	%pass offer on to level variable
if ~fakemonkey
	ontarget = eyejoytrack('holdfix', fix_point_SO, fix_radius, offer2on_time);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob fix break
	if ontarget
		pause(offer2on_time/1000)
	else
		pause(offer2on_time*rand(1)/1000)
	end
end

if ~ontarget,
	toggleobject([fix_point_SO offer2 fixpoint_offer2], 'eventmarker', 27, 'eventmarker', 33 );
	% 	putvalue(DaqInfo.TTL2, 0);
	rt = NaN;
	trialerror(3);		%broke fixation
	eventmarker(60);	%end of trial
	putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

%
%remove offer2 and hold for offeroff_time
toggleobject([fixpoint_offer2 offer2], 'eventmarker', 33);	%remove offer2

% putvalue(DaqInfo.TTL2, 0);	%pass offer on to level variable
if ~fakemonkey
	ontarget = eyejoytrack('holdfix', fix_point_SO, fix_radius, offeroff_time);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob fix break
	if ontarget
		pause(offeroff_time/1000)
	else
		pause(offeroff_time*rand(1)/1000)
	end
end
if ~ontarget,
	toggleobject(fix_point_SO, 'eventmarker', 27);
	rt = NaN;
	trialerror(3);		%broke fixation
	eventmarker(60);	%end of trial
	putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

%
%present saccade targets and wait for go signal
toggleobject([target1 target2], 'eventmarker', 35);	%display 2 targets
if ~fakemonkey
	ontarget = eyejoytrack('holdfix', fix_point_SO, fix_radius, tgton_time);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob fix break
	if ontarget
		pause(tgton_time/1000)
	else
		pause(tgton_time*rand(1)/1000)
	end
end
if ~ontarget,
	toggleobject([fix_point_SO target1 target2], 'eventmarker', 27, 'eventmarker', 36);
	rt = NaN;
	trialerror(5);		%broke fixation
	eventmarker(60);	%end of trial
	putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

%
%give go signal and check that a saccade is initiated
toggleobject(fix_point_SO, 'eventmarker', 27);
if fakemonkey && manualchoice
	[scancode rt] = getkeypress(max_reactime);
	if		(scancode==71), ontarget = 1;
	elseif	(scancode==73), ontarget = 2;
	else	ontarget = 0;
	end

	%change to LR coordint2s
	tgt1x = currTrial.currentOffers(1).xytargetposition(1);
	if		tgt1x<0
	elseif	tgt1x>0, ontarget = setdiff([1,2],ontarget);
	end

else
	if ~fakemonkey
		[ontarget rt] = eyejoytrack('holdfix', fix_point_SO, fix_radius, max_reactime);
	else
		tnow = trialtime;
		ontarget = rand(1) < 0;		%prob of no response
		if ontarget
			pause(max_reactime/1000)
		else
			pause(max_reactime*rand(1)/1000)
		end
		rt = trialtime - tnow;
	end
	if ontarget,
		toggleobject([target1 target2], 'eventmarker', 36);
		rt = NaN;
		trialerror(1);		%no response
		eventmarker(60);	%end of trial
		putvalue(DAQ.TTL{1}, 0);
		eval(['save ',filesave,' TrialRecord'])
		return
	end

	%
	%look for choice
	if ~fakemonkey
		ontarget = eyejoytrack('acquirefix', [target1 target2], [tgt_radiusL tgt_radiusR], max_sacc_time);
	else
		ontarget = rand(1) < p_nofixbreak;		%1 - prob late response
		ontarget = ceil(ontarget*2*rand(1));	%randomly select targtet
		if ontarget
			pause(max_sacc_time*rand(1)/1000)
		else
			pause(max_sacc_time/1000)
		end
	end
end

%
if rt < sessionParams.trial.toofast
	toggleobject([target1 target2], 'eventmarker', 36);
	trialerror(5);		%too fast
	eventmarker(60);	%end of trial
	putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end
%
if ~ontarget,
	toggleobject([target1 target2], 'eventmarker', 36);
	trialerror(2);		%too slow
	eventmarker(60);	%end of trial
	putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end
%
eventmarker(39);	%choice made
chosenOffer = currTrial.currentOffers(ontarget);

%
%hold target then reward
eyejoytrack('idle', sacctgt_idle);
if ~fakemonkey
	if		(ontarget==1), chtarget = target1; tgt_radius = tgt_radiusR;
	elseif	(ontarget==2), chtarget = target2; tgt_radius = tgt_radiusR;
	end
	ontarget = eyejoytrack('holdfix', chtarget, tgt_radius, tgthold_time-sacctgt_idle);
else
	ontarget = rand(1) < p_nofixbreak;	%1 - prob target fix break
	if ontarget
		pause((tgthold_time-sacctgt_idle)/1000)
	else
		pause((tgthold_time-sacctgt_idle)*rand(1)/1000)
	end
end
if ~ontarget,
	toggleobject([target1 target2], 'eventmarker', 36);
	trialerror(5);		%broke fixation
	eventmarker(60);	%end of trial
	putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

%
%if chosenOffer == null, abort
if ~chosenOffer.quantity && ontarget>0
	toggleobject([target1 target2], 'eventmarker', 36);
	trialerror(6);		%incorrect response
	eventmarker(60);	%end of trial
	putvalue(DAQ.TTL{1}, 0);
	eval(['save ',filesave,' TrialRecord'])
	return
end

trialerror(0); %correct

%
%turn off all visual stimuli
toggleobject([target1 target2], 'eventmarker', 36);

%
%deliver juice
tic;	%start preend_time before juice delivery
juiceline	= chosenOffer.good.juiceline;
ndrops		= chosenOffer.quantity;
quantum		= chosenOffer.good.quantum;
eventmarker(42+abs(juiceline)); %FIXSB
juiceDeliver(DAQ, juiceline, ndrops, quantum);
gotjuice	= 1;

%
%finish preend_time
jdt = round(1000*(toc));	%juice delivery time in ms
if	jdt > preend_time
	disp(['exceeded preend_time: juice delivery took ',num2str(jdt),' ms!'])
else
	% 	idle(preend_time - jdt);
	ontarget = eyejoytrack('idle', preend_time - jdt);
end
eventmarker(60);	%end of trial
putvalue(DAQ.TTL{1}, 0);

%
%update records
trialIndex = currTrial.trialIndex;
nTrialsHalfBlock = size(TrialRecord.User.trialsDoneThisHalfBlock,1);
TrialRecord.User.trialsDoneThisHalfBlock(trialIndex) = 1;
if sum(TrialRecord.User.trialsDoneThisHalfBlock) == nTrialsHalfBlock
	TrialRecord.User.halfBlocksDone = TrialRecord.User.halfBlocksDone + 1;
end
if TrialRecord.User.halfBlocksDone == 2
	TrialRecord.User.blocksDone = TrialRecord.User.blocksDone + 1;
end
%
currTrial.chosenOffer	= chosenOffer;
currTrial.gotjuice		= gotjuice;
TrialRecord.User.currTrial	= currTrial;
TrialRecord.User.trialRecord(trialNumber) = currTrial;
%
%update goodTrials and compute choicePattern
ngoods = TrialRecord.User.ngoods;
currentOffers = currTrial.currentOffers;
thistrial = zeros(1, ngoods+3);
thistrial(1) = trialNumber;
for ioff = 1:2
	goodId = currentOffers(ioff).goodId;
	if goodId
% 		thistrial(1+goodId) = currentOffers(ioff).quantity * sign(currentOffers(ioff).xyposition(1));
		thistrial(1+goodId) = currentOffers(ioff).quantity * sign(currentOffers(ioff).order - 1.5);		% -/+ for offer 1/2
	end
end
goodIds = [TrialRecord.User.currTrial.currentOffers.goodId];
indA = find(goodIds==1); if isempty(indA), indA = find(goodIds==0); end
tgtNumA = TrialRecord.User.currTrial.currentOffers(indA).targetnumber;
thistrial(ngoods+2) = tgtNumA;
thistrial(ngoods+3) = chosenOffer.goodId;
TrialRecord.User.goodTrials_SO = [TrialRecord.User.goodTrials_SO; thistrial];

%
%summary stats
ind = chosenOffer.order;								%order of the offer
TrialRecord.User.summaryStats.choice12(ind) = TrialRecord.User.summaryStats.choice12(ind) + 1;
% ind = (sign(chosenOffer.xyposition(1)) + 3)/2;		%LR side of the offer (obsolete)
ind = (sign(chosenOffer.xytargetposition(1)) + 3)/2;	%LR side of the target
TrialRecord.User.summaryStats.choiceLR(ind) = TrialRecord.User.summaryStats.choiceLR(ind) + 1;
TrialRecord.User.summaryStats.gotjuice(:, chosenOffer.goodId) = TrialRecord.User.summaryStats.gotjuice(:, chosenOffer.goodId) + [1; gotjuice];
TrialRecord.User.summaryStats.gotjuice_hist(trialNumber, chosenOffer.goodId) = gotjuice;
%
%choice pattern
[pairlist, table01_all_SO] = makechoicepattern_SO(TrialRecord.User.goodTrials_SO);
TrialRecord.User.choicePattern_SO = table01_all_SO;

%save TrialRecord to file
eval(['save ',filesave,' TrialRecord'])

end %%% END TWOTASK

