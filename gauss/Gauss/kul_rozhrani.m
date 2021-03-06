function fig = kul_rozhrani()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.

load kul_rozhrani
global jazyk scrsz
[X,m]=imread('kul_rozhr.jpg');
texty=[{'Kulov? rozhran?','Ohniskov? vzd?lenost [m]:','D?lka prostoru za soustavou [m]:','Index lomu 1. prost?ed?:','Index lomu 2. prost?ed?:','Pokra?uj'};
{'Ball boundary','Focus length [m]:','Length of the space behind the system [m]:','Index of refraction, 1st environment:','Index of refraction, 2nd environment:','Continue'}];

h0 = figure('Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'FileName','kul_rozhrani.m', ...
	'MenuBar','none', ...
	'Name',texty{jazyk,1}, ...
	'NumberTitle','off', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'Position',[scrsz(3)/2-155 scrsz(4)/2-125 310 250], ...
	'Tag','kulove rozhrani', ...
	'ToolBar','none');

h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[15 165 120 12], ...
	'String',texty{jazyk,2}, ...
	'Style','text', ...
	'Tag','ohniskot');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[180 165 35 15], ...
	'String','5', ...
	'Style','edit', ...
	'Tag','ohnisko');

h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[15 140 160 12], ...
	'String',texty{jazyk,3}, ...
	'Style','text', ...
	'Tag','delkatextzat');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[180 140 35 15], ...
	'String','5', ...
	'Style','edit', ...
	'Tag','delkaza');

h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[15 115 150 12], ...
	'String',texty{jazyk,4}, ...
	'Style','text', ...
	'Tag','index1t');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[180 115 35 15], ...
	'String','1', ...
	'Style','edit', ...
	'Tag','index1');

h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[15 90 150 12], ...
	'String',texty{jazyk,5}, ...
	'Style','text', ...
	'Tag','index2t');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[180 90 35 15], ...
	'String','1.5', ...
	'Style','edit', ...
	'Tag','index2');

h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'Callback','kul_rozhrani_vypocet', ...
	'ListboxTop',0, ...
	'Position',[140 40 60 20], ...
	'String',texty{jazyk,6}, ...
	'Tag','Pushbutton1');

image(X);
set(gca,'Visible','off','Position',[0.1 0.1 0.33 0.37]); 

h1 = axes('Parent',h0, ...
	'Box','on', ...
	'CameraUpVector',[0 -1 0], ...
	'Color',[1 1 1], ...
	'ColorOrder',mat1, ...
	'Layer','top', ...
	'Tag','Axes1', ...
	'Visible','off', ...
	'XColor',[0 0 0], ...
	'XLim',[0.5 230.5], ...
	'XLimMode','manual', ...
	'YColor',[0 0 0], ...
	'YDir','reverse', ...
	'YLim',[0.5 100.5], ...
	'YLimMode','manual', ...
	'ZColor',[0 0 0]);
h2 = image('Parent',h1, ...
	'CData',mat2, ...
	'Tag','Axes1Image1', ...
	'XData',[1 230], ...
	'YData',[1 100]);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',mat3, ...
	'Tag','Axes1Text4', ...
	'VerticalAlignment','cap', ...
	'Visible','off');
set(get(h2,'Parent'),'XLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[-30.29497907949791 51.10240963855422 9.160254037844386], ...
	'Rotation',90, ...
	'Tag','Axes1Text3', ...
	'VerticalAlignment','baseline', ...
	'Visible','off');
set(get(h2,'Parent'),'YLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','right', ...
	'Position',[-38.95606694560669 -8.335341365461829 9.160254037844386], ...
	'Tag','Axes1Text2', ...
	'Visible','off');
set(get(h2,'Parent'),'ZLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',mat4, ...
	'Tag','Axes1Text1', ...
	'VerticalAlignment','bottom', ...
	'Visible','off');
set(get(h2,'Parent'),'Title',h2);
if nargout > 0, fig = h0; end
