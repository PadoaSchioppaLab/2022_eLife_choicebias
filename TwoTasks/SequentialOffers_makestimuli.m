function offerImage = SequentialOffers_makestimuli(currentOffer, sessionParams, MLConfig, stimulus_type)
%
%
% mod: cps, july 2015
%
modval = 1; %rand(1)*10; %FIXSB  ???????????????

PixelsPerDegree=MLConfig.PixelsPerDegree; %FIXSB
BackgroundColor=MLConfig.SubjectScreenBackground;% FIXSB
ScreenInfo = MLConfig.Screen; %FIXSB

%decode symbol shape
symshape		= currentOffer.symbol;
if		strcmpi(symshape,'square') || strcmpi(symshape,'sqr'),	symshape = 's';
elseif	strcmpi(symshape,'circle') || strcmpi(symshape,'crc'),	symshape = 'o';
elseif	strcmpi(symshape,'diamond'),							symshape = 'd';
elseif	strcmpi(symshape,'cross'),								symshape = '+';
elseif	strcmpi(symshape,'line'),								symshape = 'l';
else	disp('image of visual stimuli not defined!!'),			dummy
end


if		strcmpi(stimulus_type,'target')
    
    % 	rgbcolor		= [1 1 1];						%all targets are white
    rgbcolor		= currentOffer.good.color;		%target color matched offer
    % 	if ~currentOffer.quantity
    % 		rgbcolor = ScreenInfo.BackgroundColor;
    % 	end
    
    symbsize_deg = sessionParams.stimulus.targetSize;
    
    fullImage_sizeDeg = 2*symbsize_deg;
    fullImage_sizePix = modval*ceil(fullImage_sizeDeg * PixelsPerDegree / modval); %FIXSB
    % fullImage_sizePix = ceil(fullImage_sizeDeg * PixelsPerDegree ); %FIXSB
    % offerImage = repmat(ones(fullImage_sizePix),[1 1 3]);%FIXSB
    offerImage = repmat(ones(abs(fullImage_sizePix)),[1 1 3]);%FIXSB
    
    offerImage(:,:,1) = offerImage(:,:,1)*BackgroundColor(1);
    offerImage(:,:,2) = offerImage(:,:,2)*BackgroundColor(2);
    offerImage(:,:,3) = offerImage(:,:,3)*BackgroundColor(3);
    offerImage=im2uint8(offerImage); %FIXSB
    
    symbsize_pix = abs(round(symbsize_deg * PixelsPerDegree));
    
    if		strcmpi(symshape,'s'),	imdata = makesquare(symbsize_pix, rgbcolor, 1,  0, BackgroundColor);
    elseif	strcmpi(symshape,'d'),	imdata = makesquare(symbsize_pix, rgbcolor, 0, 45, BackgroundColor);
    elseif	strcmpi(symshape,'+'),	imdata = makecross(symbsize_pix, rgbcolor, .30, 0, BackgroundColor);
    elseif	strcmpi(symshape,'o'),	imdata = makecircle(symbsize_pix, rgbcolor, 0, BackgroundColor);
    elseif	strcmpi(symshape,'l'),	imdata = makesquare([symbsize_pix, symbsize_pix/4], rgbcolor, 1, 30, BackgroundColor);
    end
    imgsize = size(imdata);
    
    %position imdata in offerImage
    symb_xcenter = fullImage_sizePix/2;
    symb_ycenter = fullImage_sizePix/2;
    %
    sxmin = round(symb_xcenter - imgsize(1)/2);
    symin = round(symb_ycenter - imgsize(2)/2);
    offerImage(symin:symin+imgsize(1)-1, sxmin:sxmin+imgsize(2)-1, 1) = imdata(:,:,1);
    offerImage(symin:symin+imgsize(1)-1, sxmin:sxmin+imgsize(2)-1, 2) = imdata(:,:,2);
    offerImage(symin:symin+imgsize(1)-1, sxmin:sxmin+imgsize(2)-1, 3) = imdata(:,:,3);
    
    
elseif	strcmpi(stimulus_type,'offer')
    
%          keyboard
    
    % params
    offnumber		= currentOffer.quantity;
    rgbcolor		= currentOffer.good.color;
    offposition		= currentOffer.xyposition;
    
    %outer radius for stuimuli
    outerRadius		= sessionParams.stimulus.outerRadius;
    squareSize		= sessionParams.stimulus.squareSize;
    diamondSize		= sessionParams.stimulus.diamondSize;
    crossSize		= sessionParams.stimulus.crossSize;
    circleSize		= sessionParams.stimulus.circleSize;
    lineSize		= sessionParams.stimulus.lineSize;
    symbPos_all		= sessionParams.symbPos_all_SO;
    stimIndices		= sessionParams.stimIndices;
    
    %get the actual positions of symbols and mirros reflect for offers on the left
   
%get the actual positions of symbols and mirros reflect for offers on the left
if offnumber
	symbPos = symbPos_all(stimIndices{offnumber},:);
	if offposition(1)<0
		symbPos(:,1) = -symbPos(:,1);
	end
end

if		strcmpi(symshape,'s') || strcmpi(symshape,'square')		symbsize_deg = squareSize;
elseif	strcmpi(symshape,'d') || strcmpi(symshape,'diamond')	symbsize_deg = diamondSize;
elseif	strcmpi(symshape,'+') || strcmpi(symshape,'cross')		symbsize_deg = crossSize;
elseif	strcmpi(symshape,'o') || strcmpi(symshape,'circle')		symbsize_deg = circleSize;
elseif	strcmpi(symshape,'l') || strcmpi(symshape,'line')		symbsize_deg = lineSize;
else	disp(['image of visual stimuli not defined!!']),		dummy
end

fullImage_sizeDeg = 2*(outerRadius + 2*symbsize_deg);
 fullImage_sizePix = modval*ceil(fullImage_sizeDeg * PixelsPerDegree / modval); %FIXSB
% fullImage_sizePix = ceil(fullImage_sizeDeg * PixelsPerDegree ); %FIXSB

% offerImage = repmat(ones(fullImage_sizePix),[1 1 3]);%FIXSB
offerImage = repmat(ones(abs(fullImage_sizePix)),[1 1 3]);%FIXSB
% offerImage=uint8(offerImage); %FIXSB
offerImage(:,:,1) = offerImage(:,:,1)*BackgroundColor(1);  %FIXSB
offerImage(:,:,2) = offerImage(:,:,2)*BackgroundColor(2);  %FIXSB
offerImage(:,:,3) = offerImage(:,:,3)*BackgroundColor(3);  %FIXSB
offerImage=im2uint8(offerImage); %FIXSB

% offerImage = repmat(ones(fullImage_sizePix),[1 1 4]);
% offerImage(:,:,1) = offerImage(:,:,1)*Screen.BackgroundColor(1);
% offerImage(:,:,2) = offerImage(:,:,2)*Screen.BackgroundColor(2);
% offerImage(:,:,3) = offerImage(:,:,3)*Screen.BackgroundColor(3);
% offerImage(:,:,4) = offerImage(:,:,4)*.5;
if ~offnumber
	return
else
%     keyboard
end

symbsize_pix = abs(round(symbsize_deg * PixelsPerDegree));

for isymb = 1:offnumber
	if		strcmpi(symshape,'s') || strcmpi(symshape,'square')	imdata = makesquare(symbsize_pix, rgbcolor, 1,  0, BackgroundColor);
	elseif	strcmpi(symshape,'d') || strcmpi(symshape,'diamond')imdata = makesquare(symbsize_pix, rgbcolor, 0, 45, BackgroundColor);
% 	elseif	strcmpi(symshape,'+') || strcmpi(symshape,'cross')	imdata = makecross(symbsize_pix, rgbcolor, .15, 0);	%may want to fix the function
	elseif	strcmpi(symshape,'+') || strcmpi(symshape,'cross')	imdata = makecross(symbsize_pix, rgbcolor, .15, 0, BackgroundColor);
	elseif	strcmpi(symshape,'o') || strcmpi(symshape,'circle')	imdata = makecircle(symbsize_pix, rgbcolor, 0, BackgroundColor);
	elseif	strcmpi(symshape,'l') || strcmpi(symshape,'line')	imdata = makesquare([symbsize_pix, symbsize_pix/4], rgbcolor, 1, 30, BackgroundColor);
	else	disp(['image of visual stimuli not defined!!']),	dummy
	end
	imgsize = size(imdata);

	%position imdata in offerImage
	symb_xcenter = fullImage_sizePix/2 + symbPos(isymb,1) * PixelsPerDegree;
	symb_ycenter = fullImage_sizePix/2 + symbPos(isymb,2) * PixelsPerDegree;
	%
	sxmin = abs(round(symb_xcenter - imgsize(1)/2));% FIXSB
	symin = abs(round(symb_ycenter - imgsize(2)/2));%FIXSB
	offerImage(symin:symin+imgsize(1)-1, sxmin:sxmin+imgsize(2)-1, 1) = imdata(:,:,1);
	offerImage(symin:symin+imgsize(1)-1, sxmin:sxmin+imgsize(2)-1, 2) = imdata(:,:,2);
	offerImage(symin:symin+imgsize(1)-1, sxmin:sxmin+imgsize(2)-1, 3) = imdata(:,:,3);
end
else
    disp('Unknown stim type!')
end