function varargout = Savef_v2(varargin)
% SAVEF_V2 MATLAB code for Savef_v2.fig
%      SAVEF_V2, by itself, creates a new SAVEF_V2 or raises the existing
%      singleton*.
%
%      H = SAVEF_V2 returns the handle to a new SAVEF_V2 or the handle to
%      the existing singleton*.
%
%      SAVEF_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SAVEF_V2.M with the given input arguments.
%
%      SAVEF_V2('Property','Value',...) creates a new SAVEF_V2 or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Savef_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Savef_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Savef_v2

% Last Modified by GUIDE v2.5 12-Sep-2016 12:48:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Savef_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @Savef_v2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Savef_v2 is made visible.
function Savef_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Savef_v2 (see VARARGIN)

% Choose default command line output for Savef_v2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Denis Code starts here..................
global resultp ctrlu guival guivaltag guistr guistrtag;
global DataA nvar Nonecol erreur; %#ok<NUSED>
global vgroup vship vcruiseid vtype vlat vlong vyday vflag ;
global indf1 indf2 finalok workingok wfn ffn dateok expot expod pi_names save_opts clear_savew;

if clear_savew,  guival=''; guivaltag=''; guistr=''; guistrtag='';clear_savew=0;end

% scrsz = get(0,'ScreenSize');
% figpos=get(handles.figure1,'Position');
%set(handles.figure1,'Position',[(scrsz(3)-figpos(3))/2 (scrsz(4)-figpos(4))/2 figpos(3) figpos(4) ]);
figi=handles.figure1;
CenterWindow(handles,figi,varargin{1});
save_opts_ctrl=[handles.sac_expocode,handles.sac_group,handles.sac_ship,handles.save_flag_4,handles.exclude_atm];

save_opts(1)=1;%forces expocode to be selected
for i=1:length(save_opts_ctrl),set(save_opts_ctrl(i),'Value',save_opts(i));end

if isempty(pi_names), pi_names='Dupont,D. ; Dupond,T.';end
vgroup=strtrim(vgroup);vship=strtrim(vship);pi_names=strtrim(pi_names);vcruiseid=strtrim(vcruiseid);
if ~isempty(expod), expod=strtrim(expod);end
set(handles.expot,'String', expot,'ToolTipString',expot);
set(handles.pin,'String', pi_names,'ToolTipString',pi_names);

info_mem_update(handles);
if isempty(guival) && isempty(guistr)    
    info_memory_ButtonDownFcn(-1,-1,handles);
else
    %restore the state of the GUI
    if ~isempty(guival)
        guival(find(strcmp(guivaltag,'sac_expocode')))=1;%forces expocode to be selected
        for i=1:length(guival)
            if ~isempty(strtrim(char(guivaltag(i)))) set(findobj('Tag',char(guivaltag(i))),'value',guival(i));end
        end
    end
    if ~isempty(guistr)
        for i=1:size(guistr,1) 
            if ~isempty(strtrim(char(guistrtag(i)))) set(findobj('Tag',char(guistrtag(i))),'string',strtrim(guistr(i,:)),'tooltipstring',strtrim(guistr(i,:))); end
        end
    end
end
dirn_Callback(-1, -1, handles);
wname_Callback(-1, -1, handles);
fname_Callback(-1, -1, handles);

ctrlu=[handles.cruiseidnu, handles.dirnu, handles.wnameu, handles.fnameu];
set(ctrlu,'Visible','off');

cruiseidn_Callback(-1, -1, handles);

%indf1=[];indf2=[];indf=[];
%save_flag_4_Callback(-1, -1, handles);

equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};
equ_atm=cat(1,equ,atm);
    
if save_opts(5),   indf1=find(ismember(DataA(:,vtype),equ));else    indf1=find(ismember(DataA(:,vtype),equ_atm));end
if save_opts(4)
    indf2=find(strcmp(DataA(2:end,vlat),'-999')~=1 & strcmp(DataA(2:end,vlong),'-999')~=1 & strcmp(DataA(2:end,vyday),'-999')~=1)+1;
else
    indf2=find(strcmp(DataA(2:end,vlat),'-999')~=1 & strcmp(DataA(2:end,vlong),'-999')~=1 & strcmp(DataA(2:end,vyday),'-999')~=1  & strcmp(DataA(2:end,vflag),'4')~=1)+1;
end

mess2={'N0 good EQU and/or ATM measurements!';'No good Lat, Long, YDay or Flag 2/3 found!'};mess='Final File: ';
if isempty(indf1), mess=[mess mess2{1,:}]; end
if isempty(indf2), if isempty(indf1), mess=['Final File: No Good Data to Save']; else mess=[mess mess2{2,:}];end, end
set(handles.fnamet,'String',mess,'ToolTipString',mess);
set([handles.fname,handles.chkfover,handles.chksameaswf],'enable','on');set([handles.fname2],'enable','inactive');
finalok=1;
if isempty(indf1) | isempty(indf2)
    set(handles.fnamet,'foregroundcolor','r');
    set([handles.fname,handles.fname2,handles.chkfover,handles.chksameaswf],'enable','off');
    finalok=0;
end

fs=char(filesep);

workingok=1;

% ffn=[get(handles.dirn,'String') char(fs) get(handles.fname,'String') get(handles.fname2,'String')] ;
% wfn=[get(handles.dirn,'String') char(fs) get(handles.wname,'String') get(handles.wname2,'String')] ;

uiwait(handles.figure1);


function info_mem_update(handles)
global expot expod vgroup vship vcruiseid resultp wfn ffn pi_names;

 fs=char(filesep);
 if isempty(expot), expot='00XX';end
 if isempty(expod), expod=datestr(now,'yyyymmdd');end
 if isempty(pi_names), pi_names='Dupont,D. ; Dupond,T.';end
 if isempty(wfn), wfn=[resultp char(fs) vcruiseid get(handles.wname2,'String')] ;end
 if isempty(ffn), ffn=[resultp char(fs) vcruiseid get(handles.fname2,'String')] ;end

 [~,owfn,~]=fileparts(wfn);owfn=regexprep(owfn,'_Working','');
 [~,offn,~]=fileparts(ffn);offn=regexprep(offn,'_Final','');

 infomem=sprintf('Expocode: %s %s\nGroup: %s\nShip: %s\nInvestigators: %s\nCruise ID: %s\nWorking File: %s\nFinal File: %s\nResult Files Path: %s',...
     expot,expod,vgroup,vship,pi_names,vcruiseid,owfn,offn,resultp);
 set(handles.info_memory,'String',infomem,'ToolTipString',infomem);

% --------------------------------------------------------------------

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over info_memory.
function info_memory_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to info_memory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global expot expod vgroup vship vcruiseid resultp wfn ffn;

fs=char(filesep);
if isempty(expot), expot='00XX';end
if isempty(expod), expod=datestr(now,'yyyymmdd');end
if isempty(wfn), wfn=[resultp char(fs) vcruiseid get(handles.wname2,'String')] ;end
if isempty(ffn), ffn=[resultp char(fs) vcruiseid get(handles.fname2,'String')] ;end


set(handles.expot,'String',expot,'ToolTipString',expot);set(handles.expod,'String',expod,'ToolTipString',expod);
expod_Callback(-1,-1,handles);
set(handles.groupn,'String',char(vgroup),'ToolTipString',char(vgroup));
set(handles.shipn,'String',char(vship),'ToolTipString',char(vship));
set(handles.cruiseidn,'String',char(vcruiseid),'ToolTipString',char(vcruiseid));
%Check if directory exists
set(handles.dirn,'String',resultp,'ToolTipString',resultp);
[~,owfn,~]=fileparts(wfn);owfn=regexprep(owfn,'_Working','');set(handles.wname,'String',owfn,'ToolTipString',owfn);
[~,offn,~]=fileparts(ffn);offn=regexprep(offn,'_Final','');set(handles.fname,'String',offn,'ToolTipString',offn);
dirn_Callback(-1, -1, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over changedir.
function changedir_ButtonDownFcn(hObject, eventdata, handles)
global resultp;

sp = uigetdir(resultp,'Select Dir to Save Files');
if ~isnumeric(sp)
    set(handles.dirn,'String',sp,'ToolTipString',sp);end
dirn_Callback(-1,-1,handles);




function expod_Callback(hObject, eventdata, handles)


nexpod=get(handles.expod,'String');
if (~isempty(regexp(nexpod,'[^0-9]')) | length(nexpod)~= 8)
    msgbox(sprintf('Expocode date should be composed\n\nof only 8 number with YYYYMMDD format.\n'));
    nexpod='';set(handles.expod,'String',nexpod,'ToolTipString',nexpod);
end




function dirn_Callback(hObject, eventdata, handles)
global ctrlu;
if ~exist(get(handles.dirn,'String'),'dir')
    set(handles.dirn,'ForegroundColor','red');
    set(handles.save_btn,'enable','off');
    set(handles.savebtnu,'visible','on');
else
    set(handles.save_btn,'enable','on');
    set(handles.dirn,'ForegroundColor','k');
    wname_Callback(-1,-1,handles);
end
if (hObject~=-1), set(ctrlu,'Visible','off');end



function cruiseidn_Callback(hObject, eventdata, handles)
global ctrlu;

if get(handles.chksameascrid,'Value')
    set(handles.wname,'String',get(handles.cruiseidn,'String'),'ToolTipString',get(handles.cruiseidn,'String'));
    wname_Callback(-1, -1, handles);
    fname_Callback(-1, -1, handles);
end
if sum(get(handles.dirn,'ForegroundColor')) %if red, sum is 1
    set(handles.savebtnu,'visible','on');
else
    set(handles.save_btn,'enable','on');
    set(handles.savebtnu,'visible','off');
end
if (hObject~=-1), set(ctrlu,'Visible','off');end




function wname_Callback(hObject, eventdata, handles)
%working file name control
global  ctrlu;

fs=char(filesep);
twfn=[get(handles.dirn,'String') char(fs) get(handles.wname,'String') get(handles.wname2,'String')] ;
set(handles.chkwover,'Value',1,'Visible','off');
d=dir(twfn);
if ~isempty(d)
    if ismember({[get(handles.wname,'String') get(handles.wname2,'String')]},{d.name})
      set(handles.chkwover,'Value',1,'Visible','on');
    end
end
overwrite_Callback(handles.chkwover, eventdata, handles);
if get(handles.chksameaswf,'Value') & strcmp(get(handles.chksameaswf,'Enable'),'on')
    set(handles.fname,'String',get(handles.wname,'String'),'ToolTipString',get(handles.wname,'String'));
    fname_Callback(-1, -1, handles);
end
if sum(get(handles.dirn,'ForegroundColor')) %if red, sum is 1
    set(handles.savebtnu,'visible','on');
else
    set(handles.save_btn,'enable','on');
    set(handles.savebtnu,'visible','off');
end
if (hObject~=-1), set(ctrlu,'Visible','off');end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over wname2.
function wname2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to wname2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global resultp ;
persistent chk
if isempty(chk)
      chk = 1;
      pause(0.5); %Add a delay to distinguish single click from a double click
      if chk == 1
           %fprintf(1,'\nWorking single-click.\n\n');
          chk = [];
      end
else
      chk = [];
       %fprintf(1,'\nWorking double-click.\n\n');
      if ~exist(get(handles.dirn,'String'),'dir')
          nfname=uigetfile('*_Working.csv','Select Result File...',[resultp filesep get(handles.wname,'string') '*_Working.csv']);
      else
          nfname=uigetfile('*_Working.csv','Select Result File...',[get(handles.dirn,'String') filesep get(handles.wname,'string') '*_Working.csv']);
      end
      if nfname
          nfname=regexprep(nfname,'_Working.csv','');
          set(handles.wname,'String',nfname,'ToolTipString',nfname);
          wname_Callback(-1, -1, handles);
      end
end



function fname_Callback(hObject, eventdata, handles)
%final file name control
global   ctrlu;

fs=char(filesep);
tffn=[get(handles.dirn,'String') char(fs) get(handles.fname,'String') get(handles.fname2,'String')] ;
set(handles.chkfover,'Value',1,'Visible','off');
d=dir(tffn);
if ~isempty(d)
    if ismember({[get(handles.fname,'String') get(handles.fname2,'String')]},{d.name})
      set(handles.chkfover,'Value',1,'Visible','on');
    end
end
overwrite_Callback(handles.chkfover, eventdata, handles);
if sum(get(handles.dirn,'ForegroundColor')) %if red, sum is 1
    set(handles.savebtnu,'visible','on');
else
    set(handles.save_btn,'enable','on');
    set(handles.savebtnu,'visible','off');
end
if (hObject~=-1), set(ctrlu,'Visible','off');end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over fname2.
function fname2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to fname2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global resultp ;
persistent chk
if isempty(chk)
      chk = 1;
      pause(0.5); %Add a delay to distinguish single click from a double click
      if chk == 1
          %fprintf(1,'\nFinal a single-click.\n\n');
          chk = [];
      end
else
      chk = [];
     % fprintf(1,'\nFinal a double-click.\n\n');
      
      if ~exist(get(handles.dirn,'String'),'dir')
          nfname=uigetfile('*_Final.csv','Select Result File...',[resultp filesep get(handles.fname,'string') '*_Final.csv']);
      else
          nfname=uigetfile('*_Final.csv','Select Result File...',[get(handles.dirn,'String') filesep get(handles.fname,'string') '*_Final.csv']);
      end
      
      if nfname
          nfname=regexprep(nfname,'_Final.csv','');
          set(handles.fname,'String',nfname,'ToolTipString',nfname);
          fname_Callback(-1, -1, handles);
      end
end


function KeyPressFcn(hObject, eventdata, handles)
% --- Executes on key press on cruiseidn and other controls.
global ctrlu;

set(handles.save_btn,'enable','off');
set(ctrlu(find([handles.cruiseidn handles.dirn handles.wname handles.fname]==hObject)),'Visible','on');



% --- Executes on button press in chkwover and chkfover.
function overwrite_Callback(hObject, eventdata, handles)
% hObject    handle to chkwover or chkfover (see GCBO)
global finalok workingok;
%finalok and workingok: ~enable--> 0  ~vis-->1  if vis, then value. ==>  enable and ((not vis) OR (vis&val))
%if no good data, control would not be enabled
if (strcmp(get(hObject,'Enable'),'on'))
    if ~get(hObject,'Value')
        set(hObject,'String','File will NOT be saved!','ToolTipString','File will NOT be saved!','foregroundcolor','r');
    else
        set(hObject,'String','Overwrite?','ToolTipString','Overwrite?','foregroundcolor','r');
    end
    ok= ~(strcmp(get(hObject,'Visible'),'on') & ~get(hObject,'Value'));
else
    ok=~(strcmp(get(hObject,'Enable'),'off'));
end
if (hObject==handles.chkfover), finalok=ok;end
if (hObject==handles.chkwover), workingok=ok;end


% --- Executes on button press in chksameascrid.
function chksameascrid_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    set(handles.wname,'String',get(handles.cruiseidn,'String'),'TooltipString',get(handles.cruiseidn,'String'));
    drawnow update;
    wname_Callback(-1,-1,handles);
end




% --- Executes on button press in chksameaswf.
function chksameaswf_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    set(handles.fname,'String',get(handles.wname,'String'),'TooltipString',get(handles.wname,'String'));
    drawnow update;
    fname_Callback(-1,-1,handles);
end



% --- Executes on button press in cancel_btn.
function cancel_btn_Callback(hObject, eventdata, handles)
global saveok guival guivaltag guistr guistrtag;

%remember the state of the GUI
set(hObject, 'value',0);
guictrl=findobj('-property','value');
guivaltag=get(guictrl,'Tag');
guival=get(guictrl,'value');guival=cell2mat(guival);
guictrl=findobj('style','edit');
guistrtag=get(guictrl,'tag');
guistr=get(guictrl,'string');guistr=strtrim(char(guistr));

saveok=0;
close all;


% --- Executes on button press in savef_v2.
function save_btn_Callback(hObject, eventdata, handles)
% hObject    handle to savef_v2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global resultp saveok vgroup vship vcruiseid expod savehandles ffn wfn;
global  guival guivaltag guistr guistrtag pi_names save_opts;


global DataA vtype vlat vlong vyday vflag indf indf1 indf2 indf3;

save_opts_ctrl=[handles.sac_expocode,handles.sac_group,handles.sac_ship,handles.save_flag_4,handles.exclude_atm];
for i=1:length(save_opts_ctrl),save_opts(i)=get(save_opts_ctrl(i),'Value');end

equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};
equ_atm=cat(1,equ,atm);
    
if save_opts(5),   indf1=find(ismember(DataA(:,vtype),equ)); else    indf1=find(ismember(DataA(:,vtype),equ_atm));end
if save_opts(4)
    indf2=find(strcmp(DataA(2:end,vlat),'-999')~=1 & strcmp(DataA(2:end,vlong),'-999')~=1 & strcmp(DataA(2:end,vyday),'-999')~=1)+1;
else
    indf2=find(strcmp(DataA(2:end,vlat),'-999')~=1 & strcmp(DataA(2:end,vlong),'-999')~=1 & strcmp(DataA(2:end,vyday),'-999')~=1  & strcmp(DataA(2:end,vflag),'4')~=1)+1;
end

%do not save points with flag other than 2,3 or 4 (points removed from display have flag = 50 + original flag
indf3=find(strcmp(DataA(2:end,vflag),'2')==1 | strcmp(DataA(2:end,vflag),'3')==1 | strcmp(DataA(2:end,vflag),'4')==1)+1;

indf=intersect(intersect(indf1,indf2),indf3);

%remember the state of the GUI
set(hObject, 'value',0);
guictrl=findobj('-property','value');
guivaltag=get(guictrl,'Tag');
guival=get(guictrl,'value');guival=cell2mat(guival);
guictrl=findobj('style','edit');
guistrtag=get(guictrl,'tag');
guistr=get(guictrl,'string');guistr=strtrim(char(guistr));

fs=char(filesep);
resultp=get(handles.dirn,'String');
ffn=[get(handles.dirn,'String') char(fs) get(handles.fname,'String') get(handles.fname2,'String')] ;
wfn=[get(handles.dirn,'String') char(fs) get(handles.wname,'String') get(handles.wname2,'String')] ;
expod=get(handles.expod,'String');
vgroup=get(handles.groupn,'String');
vship=get(handles.shipn,'String');
pi_names=get(handles.pin,'String');
vcruiseid=get(handles.cruiseidn,'String');
for i=1:length(save_opts_ctrl),save_opts(i)=get(save_opts_ctrl(i),'Value');end
savehandles=handles;
saveok=1;
close all;

% --- Outputs from this function are returned to the command line.
function varargout = Savef_v2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = 0;


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over GetExpod.
function GetExpod_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to GetExpod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DataA vyday vdutc dateok;

if dateok
    nexpod=char(DataA(find((str2num(char(DataA(2:end,vyday)))>0 & str2num(char(DataA(2:end,vyday)))<370),1,'first')+1,vdutc)); %1st good date looking at Year Days
    [yexpo, mexpo, dexpo, ~, ~, ~]=datevec(nexpod,'mm/dd/yy');
%    yexpo=year(expod,'mm/dd/yy');[mexpo,~]=month(expod,'mm/dd/yy');dexpo=day(expod,'mm/dd/yy');
    nexpod=[num2str(yexpo,'%04i') num2str(mexpo,'%02i') num2str(dexpo,'%02i')];
    set(handles.expod,'String',nexpod,'ToolTipString',nexpod);
end


function CenterWindow(handles,wh,ph)
global monitors;
%wh is handle of window to center
%ph is handle of parent or main window
pd=10^30;
m=1;w=getpixelposition(wh);p=getpixelposition(ph);
if size(monitors,1)>1
    for mi=1:size(monitors,1)
        if (p(1)- monitors(mi,1))>0 && (p(1) - monitors(mi,1))<pd
            m=mi ;pd=(p(1) - monitors(mi,1));
        end
    end
end
setpixelposition(wh,[monitors(m,1)+monitors(m,3)/2-w(3)/2 monitors(m,4)/2-w(4)/2 w(3) w(4)]);
