function [] = juiceDeliver(DAQ, juiceline, ndrops, quantum, varargin)
%
% synthax: juiceDeliver(DaqInfo, juiceline, ndrops, quantum, pausetime);
%

% disp(['got to juice deilvery...'])

if isempty(DAQ) 
	disp('manual j'); 
end

%pause time
if isempty(varargin)
% 	pausetime = 40;
% 	pausetime = 25;
	pausetime = 30;
else
	pausetime = varargin{1};
end

%convert to ms
%conversion: new element = old element * (set volume/measured volume)
ml_to_ms_conv = [
%	1	2	3	4	5	6	7	8	9	10		% n drops
	.92	.73	.69	.66	.65	.64	.64	.64	.65 .65		% juice line 1		
   1.16	.95	.91	.85	.86	.85	.84	.83	.82	.81		% juice line 2		%08/20/17 quantum70
	.96	.82	.75	.73	.72	.71	.70	.69	.69 .69];	% juice line 3


% 	.86	.70	.65	.60	.60	.61	.58	.58	.57 .57		% juice line 1		
% 	.81	.66	.61	.58	.58	.56	.63	.57	.55	.54		% juice line 2		%08/20/17 quantum70
% 	.83	.68	.62	.59	.59	.58	.61	.57	.56	.55];	% juice line 3

% 	.79	.65	.61	.59	.58	.57	.64	.56	.64 .56		% juice line 1		
% 	.79	.63	.59	.58	.57	.56	.63	.55	.63	.55		% juice line 2		%04/20/16 quantum70
% 	.77	.62	.58	.57	.55	.54	.65	.52	.64	.52];	% juice line 3

% 	.99	.79	.71	.7	.67	.67	.64	.64	.64 .64		% juice line 1		
% 	.92	.78	.68	.66	.64	.63	.63	.63	.63	.63		% juice line 2		%11/23/15 quantum70
% 	.98	.8	.73	.68	.68	.67	.65	.65	.64	.64];	% juice line 3

% 	.53	.52	.5	.47	.45	.43	.4	.4	.4	.4		% juice line 1			%150424 Quantum = 70	
% 	.54	.53	.51	.48	.46	.44	.4	.42	.4	.41		% juice line 2
% 	.55	.51	.49	.46	.45	.42	.4	.41	.4	.4];	% juice line 3

% ml_to_ms_conv = [
% 	1	2	3	4	5	6	7	8	9	10		% n drops
% 	.89	.81	.73	.67	.65	.63	.63	.62	.60	.59		% juice line 1		
% 	.77	.71	.64	.59	.60	.57	.56	.54	.53	.53		% juice line 2
% 	.68	.61	.57	.55	.54	.53	.53	.54	.53	.53];	% juice line 3

if ndrops > size(ml_to_ms_conv,2)
	duration = ml_to_ms_conv(juiceline, size(ml_to_ms_conv,2)) * quantum;
else
	duration = ml_to_ms_conv(juiceline, ndrops) * quantum;
end

reward_function(duration,'juiceline',juiceline,'numReward',ndrops,'pauseTime',pausetime,'eventmarker',42+juiceline);
