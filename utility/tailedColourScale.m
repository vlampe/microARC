function cmap = tailedColourScale(lo, hilo, lohi, hi, varargin)% Modify Neil's tailedrainbow function to create custom tails on any% standard colourmaps - I want to create tailed versions of the viridis% colour maps...% cmap = tailedrainbow(lo, hilo, lohi, hi);%          makes a blue-yellow-red colormap (default length 100) such that if you set%          caxis limits to [lo hi], there's a proper rainbow over [hilo lohi] and%          gradually darkening tails beyond that.%%		   e.g. tailedrainbow(0, 24, 30, 35, 0.5) will highlight the range 24..30%               but give a subtler gradient from 0..24 and 30..35, and contains%               just the colors that will be used if you put one contour every 0.5.%% cmap = tailedrainbow(lo, hilo, lohi, hi, n);%          sets the length n of the colormap. Default is 100.%% cmap = tailedrainbow(n);%          just makes the rainbow, no tails.%% cmap = tailedrainbow(vec);%			where vec = [lo hilo lohi hi] or [lo hilo lohi hi n] is okay too.%% Neil Banas, 2014% centerValue = 0.87; % value as in brightness or whiteness% endValue = 0.65;% tailEndValue = 0.3;%% blueHue = 0.55;% yellowHue = 1/6;% redHue = 0.04;%% greenHue = 1/3;% greenDesatWidth = 0.67;% greenDesatAmount = 0.5;% EXTEND THIS TO ALLOW FOR EXTRA FUNCTION ARGUMENTS LATER... SEE COMMENTED% STUFF BELOW...failed = false;if nargin == 4 % default settings    n = 100;    colourMap = 'viridis';elseif nargin == 5    if any(strcmp(varargin{1},{'viridis' 'inferno' 'magma' 'plasma'})) % specified colour map        colourMap = varargin{1};        n = 100;    else % specified interval (number of colours)        colourMap = 'viridis';        interval = varargin{1};        n = (hi-lo) / interval;    endelseif nargin == 6    if any(strcmp(varargin{1},{'viridis' 'inferno' 'magma' 'plasma'}))        colourMap = varargin{1};        interval = varargin{2};        n = (hi-lo) / interval;    elseif any(strcmp(varargin{2},{'viridis' 'inferno' 'magma' 'plasma'}))        colourMap = varargin{2};        interval = varargin{1};        n = (hi-lo) / interval;    endelse    failed = true;    disp('Invalid input arguments to "tailedColourScale"')end% if length(lo) >= 4% 	if length(lo)==5% 		n = (lo(4)-lo(1))/lo(5);% 	else% 		n = 100;% 	end% 	hi = lo(4);% 	lohi = lo(3);% 	hilo = lo(2);% 	lo = lo(1);% else% 	if nargin==5% 		n = round((hi-lo) / interval);% 	elseif nargin==4% 		n = 100;% 	elseif nargin==1% 		n = lo;% 		lo = 0;% 		hilo = 0;% 		lohi = 1;% 		hi = 1;% 	end% end%% if any(sort([lo hilo lohi hi]) ~= [lo hilo lohi hi])% 	disp('tailedrainbow: color ranges out of order. Returning something sensible instead.');% 	cmap = tailedrainbow(100);% 	return% endif failed    cmap = nan;else            nmain = round((lohi-hilo) / (hi-lo) * n);    cmapMain = feval(colourMap,nmain); % focal colours        % make tails ----------------------------------        nLowTail = round((hilo-lo) / (hi-lo) * n);    if nLowTail >= 1        lowCol = rgb2hsv(cmapMain(1,:));        hue = repmat(lowCol(1), [nLowTail 1]);        sat = repmat(lowCol(2), [nLowTail 1]);        val = linspace(4/5*lowCol(3), lowCol(3), nLowTail+1);        val = val(1:end-1);        lowTail = hsv2rgb([hue(:) sat(:) val(:)]);    else        lowTail = [];    end        nHighTail = round((hi-lohi) / (hi-lo) * n);    if nHighTail >= 1        highCol = rgb2hsv(cmapMain(end,:));        hue = repmat(highCol(1), [nHighTail 1]);        sat = repmat(highCol(2), [nHighTail 1]);        val = linspace(highCol(3), 4/5*highCol(3), nHighTail+1);        val = val(2:end);        highTail = hsv2rgb([hue(:) sat(:) val(:)]);    else        highTail = [];    end        cmap = [lowTail; cmapMain; highTail];    nmap = size(cmap,1);        if nmap ~= n % if the rounding has left the map the wrong length, interp it down to the right length        cmap = interp1(linspace(0,1,nmap)', cmap, linspace(0,1,n)');    end    end