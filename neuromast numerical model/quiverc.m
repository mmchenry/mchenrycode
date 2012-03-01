function varargout = quiverc(varargin)
% function quiverc(x,y,u,v,...)
%	or     quiverc(u,v,...)
%	or	   quiverc(...,col,...)
%
% Plots a vector field at coordinates (x,y) with vector length (u,v).  Or
% plots on an evenly spaced grid the vector matrices (u,v).  Can color the
% vectors according to another matrix, col, which must be the same size as
% u and v.  You can also pass a single character color as in other plotting
% routines.  Currently does not handle RGB colors correctly.  Automatically
% scale the vectors so that the 95% percentile vector length corresponds to
% the average spacing between vectors.  Vector head size scales with length,
% also.
%
% Takes a lot of options to customize the display of the vectors.  Options can
% come in any order, and usually take a numeric argument after the option
% name.  All options have a long name and a short name shown in parenthesis.
%
%	'ScaleFactor' ('s'),sf - Multiplies the autoscale length by sf.  Default: 1
%	'RelScale' ('rs'),rs - Changes which percentile corresponds to the average
%		spacing.  Percentiles are from 0 to 1.  Default: 0.95
%	'AbsScale' ('as'),scale - Absolute scaling value.  Vector display lengths are
%		created by multiplying the vector length by this value.
%	'ScaleRange' ('sr'),[lo hi] - Changes the range over which vector display length
%		scales by true length.  Expressed in percentiles from 0 to 1.  Lengths
%		below lo are shown as lo in length, and lengths above hi are no longer.
%	'AbsScaleRange' ('asr'),[lo hi] - Same as ScaleRange, but expressed in absolute
%		length, not percentiles.
%	'HeadSize' ('hs'),hsz - Changes the head size as a percentage from 0 to 1 of 
%		the total length of the vector.  Default: 0.4.
%	'HeadRange' ('hr'),[hlo hhi] - Like 'ScaleRange', but only affects the scaling
%		of the head.  A typical use is to have small heads remain even when
%		the vectors get really small, so that you can still see the direction:
%		'HeadRange',[0.2 1], which means that all vectors shorter than the 20th
%		percentile will still have a head.
%	'AbsHeadRange' ('ahr'),[hlo hhi] - Same as 'HeadRange', but takes absolute lengths,
%		not percentiles.
%	'Truncate' ('t'),maxlen - Old option provided for backwards compatibility.
%		Cuts off scaling at maxlen percentile.  Use ScaleRange instead.
%	'Show' ('sh'),fac - Shows fac times fewer vectors (fac < 1).  Mostly use
%		multiples of 1/4, which corresponds to half the number of vectors in
%		both x and y.
%	'NoHeads' ('nh') - Doesn't show vector heads.
%	'NoTails' ('nt') - Doesn't show vector tails.
%
% Copyright (c) 2003 - Eric Tytell and Harvard University

% get the size of all the arguments
sz(:,1) = cellfun('size',varargin,1)';
sz(:,2) = cellfun('size',varargin,2)';

% reject any 3D matrices
if (any(cellfun('size',varargin,3) > 1)),
	error('Cannot handle 3D matrices.');
end;

% check for various combinations of x,y,u,v
if (all(sz(1,:) == sz(2,:)) & all(sz(1,:) == sz(3,:)) & all(sz(1,:) == sz(4,:))),
	% quiverc(x,y,u,v,...)
	x = varargin{1};
	y = varargin{2};
	u = varargin{3};
	v = varargin{4};
	p = 5;
elseif (all(sz(1,:) == sz(2,:))),
	% quiverc(u,v,...)
	u = varargin{1};
	v = varargin{2};
	[x,y] = meshgrid(1:size(u,2),1:size(u,1));
	
	p = 3;
elseif (all(sz(3,:) == sz(4,:)) & (prod(sz(1,:)) == sz(3,1)) & (prod(sz(2,:)) == sz(3,2))),
	% quiverc(x(1,:),y(:,1),u,v,...)
	[x,y] = meshgrid(varargin{1},varargin{2});
	u = varargin{3};
	v = varargin{4};
	p = 5;
end;

col = [];
if ((p <= nargin) & (all(sz(p,:) == size(u)))),
	col = varargin{p};
	p = p+1;
	isColMatrix = 1;
end;

% process all the options
options = varargin(p:end);
[scaleFac,options,m1] = matchOption(options,'ScaleFactor','s',1);
[relScale,options,m2] = matchOption(options,'RelScale','rs',0.95);
[absScale,options,m3] = matchOption(options,'AbsScale','as',[]);
if (m3 & (m1 | m2)),
	warning('AbsScale option overrides ScaleFactor or RelScale option.');
end;
[scaleRange,options,m1] = matchOption(options,'ScaleRange','sr',[0 1]);
[absScaleRange,options,m2] = matchOption(options,'AbsScaleRange','asr',[]);
if (m1 & m2),
	warning('AbsScaleRange option overrides ScaleRange option.');
end;
[headSize,options] = matchOption(options,'HeadSize','hs',0.4);
[headRange,options,m1] = matchOption(options,'HeadRange','hr',[0 1]);
[absHeadRange,options,m2] = matchOption(options,'AbsHeadRange','ahr',[]);
if (m1 & m2),
	warning('AbsHeadRange option overrides HeadRange option.');
end;

[trunc,options,m1] = matchOption(options,'Truncate','t',1);
if (m1),
	warning('Obsolete option ''truncate''.  Use ''ScaleRange'' instead.');
	scaleRange = [0 trunc];
end;
[shw,options] = matchOption(options,'Show','sh',1);

[noheads,options] = matchOption(options,'NoHeads','nh');
[notails,options] = matchOption(options,'NoTails','nt');
if (noheads & notails),
	error('You must either display heads or tails.');
end;

% Check for color matrix or color string
if (~isempty(options) & ischar(options{1})),
	if (~isempty(col)),
		warning('Color option overrides color variable.');
	end;
	col = options{1};
	options = options(2:end);
	isColMatrix = 0;
	p = p+1;
end;
if (isempty(col)),
	col = 'b';
	isColMatrix = 0;
end;

% error if there's anything left
if (~isempty(options)),
	isc = cellfun('isclass',options,'char');
	str = sprintf('%s, ',options{isc});
	
	error(['Unrecognized options: ' str '.']);
end;

switch get(gca,'NextPlot'),
case 'replace',
	cla reset;
case 'replacechildren',
	cla;
end;

% definition of the arrow shape
AR    = range(y(:))./range(x(:));
headx = [0.85 	-0.15	0		-0.15 	0.85]';
heady = [0		-0.4	0		0.4		0]' .*AR;
tailx = [0		1]';
taily = [0		0]';

% reduce the number of vectors according to the Show option
if (any(size(x) == 1)),
	% random vector positions
	N = length(x);
	
	k = round(linspace(1,N,shw*N));
	x = x(k);
	y = y(k);
	u = u(k);
	v = v(k);
	if (isColMatrix),
		col = col(k);
	end;
	
	m = sqrt(length(x));
	n = m;
else,
	% gridded vector positions
	[m,n] = size(x);
	
	i = round(linspace(1,m,sqrt(shw)*m));
	j = round(linspace(1,n,sqrt(shw)*n));
	x = x(i,j);
	y = y(i,j);
	u = u(i,j);
	v = v(i,j);
	if (isColMatrix),
		col = col(i,j);
	end;
	
	[m,n] = size(x);
end;

%%%%%% Set up scaling

% first get the spacing between points in x and y
dx = diff([min(x(:)) max(x(:))]);
dy = diff([min(y(:)) max(y(:))]);
d = (dx + dy)/(n+m);	% weight spacing by the number of points in each direction

% Remove any NaNs
k = find(isfinite(x) & isfinite(y) & isfinite(u) & isfinite(v));
x = shiftdim(x(k))';
y = shiftdim(y(k))';
u = shiftdim(u(k))';
v = shiftdim(v(k))';
if (isColMatrix),
	col = shiftdim(col(k))';
end;
N = length(k);

% Set up the vectors so that they have a length (len) and a direction (u,v)
len = sqrt(u.^2 + v.^2);
u = u./len;
v = v./len;

% Autoscale
if (~isempty(absScale)),
	scale = absScale;
else,
	long = prctile(len(:),relScale*100);
	scale = d/long * scaleFac;
end;

slen = len;

% Define scale range
if (~isempty(absScaleRange)),
	isShort = find(len < absScaleRange(1));
	isLong = find(len > absScaleRange(2));
	slen(isShort) = absScaleRange(1);
	slen(isLong) = absScaleRange(2);
else
	short = prctile(len(:),scaleRange(1)*100);
	long = prctile(len(:),scaleRange(2)*100);
	isShort = find(len < short);
	isLong = find(len > long);
	slen(isShort) = short;
	slen(isLong) = long;
end;

headlen = len*headSize;
if (~isempty(absHeadRange)),
	isShortHead = find(len < absHeadRange(1));
	isLongHead = find(len > absHeadRange(2));
	headlen(isShortHead) = absHeadRange(1)*headSize;
	headlen(isLongHead) = absHeadRange(2)*headSize;
else
	short = prctile(len(:),headRange(1)*100);
	long = prctile(len(:),headRange(2)*100);
	
	isShortHead = find(len < short);
	isLongHead = find(len > long);
	headlen(isShortHead) = short*headSize;
	headlen(isLongHead) = long*headSize;
end;

taillen = slen - headlen;
taillen(taillen < 0) = 0;
nonZeroTail = find(taillen > 0);

if (scale ~= 0),
	taillen = taillen*scale;
	headlen = headlen*scale;
else
	taillen = repmat(d*(1-headSize),size(len));
	headlen = repmat(d*headSize,size(len));
end;

%%%% Construct vectors

if (noheads),
	ax = repmat(tailx,[1 N]);
	ay = repmat(taily,[1 N]);
	
	ax = ax.*repmat((headlen + taillen),[size(ax,1) 1]);
	ay = ay.*repmat((headlen + taillen),[size(ax,1) 1]);
	
	tailind = 1:length(tailx);
elseif (notails),
	ax = repmat(headx,[1 N]);
	ay = repmat(heady,[1 N]);
	
	ax = ax.*repmat((headlen + taillen),[size(ax,1) 1]);
	ay = ay.*repmat((headlen + taillen),[size(ax,1) 1]);
	
	headind = 1:length(headx);
else,
	hx = repmat(headx,[1 N]);
	hy = repmat(heady,[1 N]);
	tx = repmat(tailx,[1 N]);
	ty = repmat(taily,[1 N]);
	
	tx = tx.*repmat(taillen,[size(tx,1) 1]);
	ty = ty.*repmat(taillen,[size(tx,1) 1]);
	hx = hx.*repmat(headlen,[size(hx,1) 1]) + repmat(taillen,[size(hx,1) 1]);
	hy = hy.*repmat(headlen,[size(hx,1) 1]);
	
	ax = [tx; hx];
	ay = [ty; hy];
	
	tailind = 1:length(tailx);
	headind = (1:length(headx)) + length(tailx);
end;

n = size(ax,1);

ptx = (repmat(u,[n 1]).*ax - repmat(v,[n 1]).*ay) + repmat(x,[n 1]);
pty = (repmat(v,[n 1]).*ax + repmat(u,[n 1]).*ay) + repmat(y,[n 1]);

if (isColMatrix),
	col = repmat(col,[size(ptx,1) 1]);
end;

if (~noheads),
	if (isColMatrix),
		h(1) = patch(ptx(headind,:),pty(headind,:),col(headind,:),'EdgeColor','none');
	else
		h(1) = patch(ptx(headind,:),pty(headind,:),col,'EdgeColor','none');
	end;		
end;
if (~notails),
	if (isColMatrix),
		h(2) = patch(ptx(tailind,nonZeroTail),pty(tailind,nonZeroTail),'b','FaceColor','none','FaceVertexCData',flatten(col(tailind,nonZeroTail)),'EdgeColor','flat');
	else
		h(2) = patch(ptx(tailind,nonZeroTail),pty(tailind,nonZeroTail),'b','FaceColor','none','EdgeColor',col(1));
	end;
end;

if (nargout == 1),
	varargout = {h};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% matchOption
function [val,opt,ismatch] = matchOption(opt,name,nick,def)

if (nargin == 3),
	isdef = 0;
else
	isdef = 1;
end;

q = find(cellfun('isclass',opt, 'char'));
ind = strmatch(lower(name),lower(opt(q)),'exact');
if (isempty(ind)),
	ind = strmatch(lower(nick),lower(opt(q)),'exact');
end;
ind = q(ind);

if (isdef),
	if (isempty(ind)),
		val = def;
		ismatch = 0;
	else
		if ((ind == length(opt)) | ~isnumeric(opt{ind+1})),
			error(sprintf('Option %s requires a numeric parameter.',name));
		end;
		val = opt{ind+1};
		ismatch = 1;
		
		opt = opt([1:ind-1 ind+2:end]);
	end;
else
	if (isempty(ind)),
		val = 0;
		ismatch = 0;
	else
		val = 1;
		opt = opt([1:ind-1 ind+1:end]);
	end;
end;
