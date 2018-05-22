% function [lat] = GetStepLatencies(line,fs,time)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% line           = line-nanmean(line);
% line           = line/nanstd(line);

if size(line,2)<size(line,1)
    line            = line';
end
line                = double(line);
line(isnan(line))   = nanmean(line);

if size(line_oppo,2)<size(line_oppo,1)
    line_oppo            = line_oppo';
end
line_oppo           = double(line_oppo);
line(isnan(line_oppo))   = nanmean(line_oppo);

BT_PreProcess_GetStepLatencies_firstplot;

%% Create Sliders and Edit Boxes
%  for STATELEVEL 1 (SL)
tsl1 = uicontrol('style','text', 'str','LowLevel (amplitude):', 'fontsize',14, 'unit','norm');
set(tsl1,'pos',[.12 .17 .1 .03])

slsl1 = uicontrol('style','slider','Tag','sl1-slide', ...
    'Min', min(line), 'Max', -0.05, 'Value', -2, 'SliderStep',[0.01 0.01], ...
    'unit','norm');
set(slsl1,'pos',[.25 .165 .5 .03])

edsl1 = uicontrol('style','edit','Tag','sl1-edit', 'str','-2','fontsize',14,'unit','norm');
set(edsl1,'pos',[.8 .17 .05 .03])

%  for STATELEVEL 2 (SL2)
tsl2 = uicontrol('style','text', 'str','HighLevel (amp):', 'fontsize',14, 'unit','norm');
set(tsl2,'pos',[.12 .135 .1 .03])

slsl2 = uicontrol('style','slider','Tag','sl2-slide', ...
    'Min', 0.05, 'Max', max(line), 'Value', 2, 'SliderStep',[0.01 0.01], ...
    'unit','norm');
set(slsl2,'pos',[.25 .13 .5 .03])

edsl2 = uicontrol('style','edit','Tag','sl2-edit', 'str','2','fontsize',14,'unit','norm');
set(edsl2,'pos',[.8 .132 .05 .03])

% for ARTEFACTTHRESHOLD SLOW FALL
tsf = uicontrol('style','text', 'str','Thresh SlowFall:', 'fontsize',14, 'unit','norm');
set(tsf,'pos',[.12 .095 .1 .03])

slsf = uicontrol('style','slider','Tag','sf-slide', ...
    'Min', 0.01, 'Max', .7, 'Value', .2, 'SliderStep',[0.01 0.01], ...
    'unit','norm');
set(slsf,'pos',[.25 .09 .15 .03])

edsf = uicontrol('style','edit','Tag','sf-edit', 'str','.2','fontsize',14,'unit','norm');
set(edsf,'pos',[.405 .095 .05 .03])

% for ARTEFACTTHRESHOLD STEP REPETITION
tsr = uicontrol('style','text', 'str','Thresh StepRep:', 'fontsize',14, 'unit','norm');
set(tsr,'pos',[.455 .095 .1 .03])

slsr = uicontrol('style','slider','Tag','sr-slide', ...
    'Min', 0.01, 'Max', 1, 'Value', .3, 'SliderStep',[0.01 0.01], ...
    'unit','norm');
set(slsr,'pos',[.55 .09 .12 .03])

edsr = uicontrol('style','edit','Tag','sr-edit', 'str','.3','fontsize',14,'unit','norm');
set(edsr,'pos',[.8 .095 .05 .03])

tsrp = uicontrol('style','text', 'str','Keep 2d:', 'fontsize',14, 'unit','norm');
set(tsrp,'pos',[.67 .095 .05 .03])

slsrp = uicontrol('style','checkbox','Tag','srp-checkbox','Value',1);
set(slsrp,'pos',[1040 75 60 20])


% for ADAPGT RANGE
%fig = figure('unit','norm','pos',[.01 .05 .98 .85]); 

tar = uicontrol('style','text', 'str','Adapt Range:', 'fontsize',14, 'unit','norm');
set(tar,'pos',[.12 .05 .1 .03])

slar = uicontrol('style','slider','Tag','ar-slide', ...
    'Min', 0, 'Max', 40, 'Value', 0, 'SliderStep',[0.025 0.25], ...
    'unit','norm');
set(slar,'pos',[.25 .05 .15 .03])

edar = uicontrol('style','edit','Tag','ar-edit', 'str','0','fontsize',14,'unit','norm');
set(edar,'pos',[.405 .055 .05 .03])

% for ADAPGT RANGE
get_ar_edit_handle  = 'har=findobj(''tag'',''ar-edit''); ';
update_ar_variable  = 'ar=str2num(get(har,''str'')); ';

update_edit_box    = 'set(har,''str'',num2str(get(gcbo,''value''))); ';
update_plot        = 'BT_PreProcess_GetStepLatencies_updateplot;';

set(edar,'call',[get_ar_edit_handle update_ar_variable update_plot])
set(slar,'call',[get_ar_edit_handle update_edit_box update_ar_variable update_plot])






% Manually deleate a value
sldelete = uicontrol('style','pushbutton','Tag','sldelete-button','String','Delete Step');
set(sldelete,'pos',[1200 600 60 20])

sladd = uicontrol('style','pushbutton','Tag','sladd-button','String','Add Step');
set(sladd,'pos',[1200 550 60 20])

% Finish Button
slfinish = uicontrol('style','pushbutton','Tag','finish-button','String','Finished');


%% Callbacks
%  for STATELEVEL 1 (SL1)
get_sl1_edit_handle  = 'hsl1=findobj(''tag'',''sl1-edit''); ';
update_sl1_variable  = 'sl1=str2num(get(hsl1,''str'')); ';

update_edit_box    = 'set(hsl1,''str'',num2str(get(gcbo,''value''))); ';
update_plot        = 'BT_PreProcess_GetStepLatencies_updateplot;';

set(edsl1,'call',[get_sl1_edit_handle update_sl1_variable update_plot])
set(slsl1,'call',[get_sl1_edit_handle update_edit_box update_sl1_variable update_plot])

%  for STATELEVEL 2 (SL2)
get_sl2_edit_handle  = 'hsl2=findobj(''tag'',''sl2-edit''); ';
update_sl2_variable  = 'sl2=str2num(get(hsl2,''str'')); ';

update_edit_box    = 'set(hsl2,''str'',num2str(get(gcbo,''value''))); ';
update_plot        = 'BT_PreProcess_GetStepLatencies_updateplot;';

set(edsl2,'call',[get_sl2_edit_handle update_sl2_variable update_plot])
set(slsl2,'call',[get_sl2_edit_handle update_edit_box update_sl2_variable update_plot])

% for ARTEFACTTHRESHOLD SLOW FALL
get_sf_edit_handle  = 'hsf=findobj(''tag'',''sf-edit''); ';
update_sf_variable  = 'sf=str2num(get(hsf,''str'')); ';

update_edit_box    = 'set(hsf,''str'',num2str(get(gcbo,''value''))); ';
update_plot        = 'BT_PreProcess_GetStepLatencies_updateplot;';

set(edsf,'call',[get_sf_edit_handle update_sf_variable update_plot])
set(slsf,'call',[get_sf_edit_handle update_edit_box update_sf_variable update_plot])


% for ARTEFACTTHRESHOLD STEP REPETITION
get_sr_edit_handle  = 'hsr=findobj(''tag'',''sr-edit''); ';
update_sr_variable  = 'sr=str2num(get(hsr,''str'')); ';

update_edit_box    = 'set(hsr,''str'',num2str(get(gcbo,''value''))); ';
update_plot        = 'BT_PreProcess_GetStepLatencies_updateplot;';

set(edsr,'call',[get_sr_edit_handle update_sr_variable update_plot])
set(slsr,'call',[get_sr_edit_handle update_edit_box update_sr_variable update_plot])


get_srp_edit_handle     = 'hsrp=findobj(''tag'',''srp-checkbox''); ';
update_srp_variable     = 'srp=get(hsrp,''Value''); ';
update_plot             = 'BT_PreProcess_GetStepLatencies_updateplot;';
set(slsrp,'call',[get_srp_edit_handle update_srp_variable update_plot])

% Manually delete a value
get_sladd_edit_handle     = 'hsadd=findobj(''tag'',''sladd-button''); ';
update_sadd_variable     = 'button=get(hsadd,''String''); ';
update_plot             = 'BT_PreProcess_GetStepLatencies_updateplot_manual;';
set(sladd,'call',[get_sladd_edit_handle update_sadd_variable update_plot])


get_sldelete_edit_handle     = 'hsdelete=findobj(''tag'',''sldelete-button''); ';
update_sdelete_variable     = 'button=get(hsdelete,''String''); ';
update_plot             = 'BT_PreProcess_GetStepLatencies_updateplot_manual;';
set(sldelete,'call',[get_sldelete_edit_handle update_sdelete_variable update_plot])

% Finish
get_finish_edit_handle      = 'hfinish=findobj(''tag'',''finish-button''); ';
update_sfinish_variable     = 'lat=steps_t; ';
close_plot                  = 'close all';
set(slfinish,'call',[get_finish_edit_handle update_sfinish_variable close_plot])

waitfor(slfinish,'call');
