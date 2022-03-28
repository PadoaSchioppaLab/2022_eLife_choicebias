function crs = makecross(sz, rgb, varargin)
%SYNTAX:
%        sqr = makecross(size, rbg, [thickness], [rotation], [bgcolor])
%
% Size is in pixels. Thickness is optional and it is a fraction of the
% size. Rotation is optional and is in degrees. Default is rotation = 0,
% which gives a '+'. Bgcolor is optional and is specified as an RGB
% triplet. Must specify thickness and rotation to specify bgcolor.
%

if length(sz) == 1,
	xs = round(sz);
	ys = round(sz);
else
	xs = round(sz(1));
	ys = round(sz(2));
	% if different parity reduce larger size (necessary so that 
	% the 2 arms of the cross have the same thickness)
	if floor((xs + ys)/2) < (xs + ys)/2
		if xs < ys, ys = ys-1;
		elseif xs > ys, xs = xs - 1;
		end
	end
end

crs = zeros(ys, xs);

if ~isempty(varargin),
	if length(varargin) == 1
		linwidth = varargin{1};
		theta = 0;
		bgcolor = [0 0 0];
	elseif length(varargin) == 2
		linwidth = varargin{1};
		theta = varargin{2};
		bgcolor = [0 0 0];
	elseif length(varargin) == 3
		linwidth = varargin{1};
		theta = varargin{2};
		bgcolor = varargin{3};
	end
else
	linwidth = .15;
	theta = 0;
	bgcolor = [0 0 0];
end

% we make sure that the line width in pixels is even/odd depending on
% whether xs and ys are even/odd
if floor(xs/2) == xs/2
	lw = 2 * round(linwidth * min([xs ys])/2);
elseif floor(xs/2) < xs/2
	lw = floor(linwidth * min([xs ys]));
	if floor(lw/2) == lw/2
		lw = lw + 1;
	end
end
for (i = ys/2-lw/2+1 : ys/2+lw/2)	crs(i,:) = 1;	end
for (j = xs/2-lw/2+1 : xs/2+lw/2)	crs(:,j) = 1;	end

if theta > 0,
	mtheta = mod(theta, 90);
	if mtheta > 45,
		mtheta = abs(mtheta - 90);
	end
	gfactor = 1 + (abs(tan(pi*mtheta/180)));
	bxs = round(gfactor*xs);
	bys = round(gfactor*ys);
	x1 = round(xs*(gfactor-1)/2);
	y1 = round(ys*(gfactor-1)/2);
	bcrs = zeros(bys, bxs);
	bcrs(y1+1:y1+ys, x1+1:x1+xs) = crs;
	crs = imrotate(bcrs, theta, 'crop');
end

if max(max(crs)) > 0,
	crs = crs./max(max(crs));
end
crs = repmat(crs, [1 1 3]);
crs(:, :, 1) = crs(:, :, 1).*(rgb(1) - bgcolor(1)) + bgcolor(1);
crs(:, :, 2) = crs(:, :, 2).*(rgb(2) - bgcolor(2)) + bgcolor(2);
crs(:, :, 3) = crs(:, :, 3).*(rgb(3) - bgcolor(3)) + bgcolor(3);

