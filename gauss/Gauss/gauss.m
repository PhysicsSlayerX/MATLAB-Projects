function fig = gauss()
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

load gauss
%uvodni obrazovka
global scrsz;
scrsz = get(0,'ScreenSize');
close all;

h0 = figure('Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'FileName','gauss.m', ...
	'MenuBar','none', ...
	'Name','Gaussian beam Gauss�v svazek 1.0 ', ...
	'NumberTitle','off', ...
	'PaperPosition',[20 200 550 400], ...
	'PaperUnits','points', ...
	'Position',[scrsz(3)/2-250 scrsz(4)/2-135 500 270], ...
	'Tag','main', ...
	'ToolBar','none');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'Callback','inform(2)', ...
	'ListboxTop',0, ...
	'Position',[93.75 15 60 20.25], ...
	'String','Continue', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'FontSize',24, ...
	'ListboxTop',0, ...
	'Position',[45 120 300 60], ...
	'String','Transition of Gaussian beam trough an optical system', ...
	'Style','text', ...
	'Tag','Nazev');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'FontSize',24, ...
	'ListboxTop',0, ...
	'Position',[45 45 300 60], ...
	'String','Pr�chod Gaussova svazku optickou soustavou', ...
	'Style','text', ...
	'Tag','Nazev');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'Callback','inform(1)', ...
	'ListboxTop',0, ...
	'Position',[206.25 15 60 19.5], ...
	'String','Pokra�uj ', ...
	'Tag','Pushbutton1');
if nargout > 0, fig = h0; end