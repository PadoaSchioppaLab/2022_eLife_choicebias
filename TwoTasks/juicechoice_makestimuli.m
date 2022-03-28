function offerImage = juicechoice_makestimuli(currentOffer, sessionParams, mlconfig)
%  close all
%  keyboard
%
%
%
%disp(Screen) %FIXSB 11/16
% modval = ScreenInfo.ModVal;
modval = 1; %rand(1)*10; %FIXSB  ???????????????

PixelsPerDegree=mlconfig.PixelsPerDegree; %FIXSB
BackgroundColor=mlconfig.SubjectScreenBackground;% FIXSB
Screen = mlconfig.Screen; %FIXSB



%disp(currentOffer)
% params
%
offnumber		= currentOffer.quantity;
rgbcolor		= currentOffer.good.color;
try
symshape		= currentOffer.good.symbol;
catch
 symshape		= currentOffer.symbol;
end
offposition		= currentOffer.xyposition;
% %
% outerRadius	= sessionParams.outerRadius;		%outer radius for stimuli
% squareSize	= sessionParams.squareSize;
% diamondSize	= sessionParams.diamondSize;
% crossSize	= sessionParams.crossSize;
% circleSize	= sessionParams.circleSize;
% lineSize	= sessionParams.lineSize;
% symbPos_all	= sessionParams.symbPos_all;
% stimIndices	= sessionParams.stimIndices;
%
outerRadius	= sessionParams.stimulus.outerRadius;		%outer radius for stimuli
squareSize	= sessionParams.stimulus.squareSize;
diamondSize	= sessionParams.stimulus.diamondSize;
crossSize	= sessionParams.stimulus.crossSize;
circleSize	= sessionParams.stimulus.circleSize;
lineSize	= sessionParams.stimulus.lineSize;
symbPos_all	= sessionParams.symbPos_all_JC;
stimIndices	= sessionParams.stimIndices;
%
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





% close all
% keyboard

% subplot(1,2,1);  imshow(offerImage) ; subplot(1,2,2); mglimage(offerImage); axis square
 



