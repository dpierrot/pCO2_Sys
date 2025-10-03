function varargout = pCO2_Sys_v140(varargin) 
%  PCO2_SYS_V140 Application M-file for pco2_sys_v140.fig
%    FIG = PCO2_SYS_V140 launch pco2_sys_v140 GUI.
%    PCO2_SYS_V140('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 16-Jan-2025 10:37:31

clearvars -except varargin;
global f1 f2 hp hp2 zoomon motion cmh fig off_on mapfn map_loaded shorelines;
global sub_flag vartitlesc vartitles titles inihead nvar dateok posok;
global projDir datap  headp resultp hfile  sysinip sysinif fversionn where_ini def_folder map_folder where_m where_exe where_xml;
global v_compatible;
v_compatible ={'125','126 127 130 140'};

mapfn = ['gshhs_l.b';'gshhs_h.b';'gshhs_f.b'];

global  dotsize;dotsize=8;
global monitors;monitors=get(0,'MonitorPositions');
global Dmonitor;
% p=strsplit(path,';');
% fprintf('%s\n',p{:})

if nargin == 0  % LAUNCH GUI

%---------------------------set default monitor-----------------------------------------------------------------
desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
desktopMainFrame = desktop.getMainFrame;
% Set default figure position to left hand monitor
Lmonitor=find(monitors(:,1)==min(monitors(:,1)));
Rmonitor=find(monitors(:,1)==max(monitors(:,1)));
set(0,'DefaultFigurePosition',monitors(Lmonitor,1:4));

Dmonitor=1;
if ~isdeployed % Running from MATLAB.
    if size(monitors,1)>1 % If there are more than one monitor
        if desktopMainFrame.getX >= max(monitors(:,1)) % If the command window is in the right hand monitor
            % Set default figure position to left hand monitor
            set(0,'DefaultFigurePosition',monitors(Lmonitor,1:4));
            Dmonitor=Lmonitor;
        else % Then the command window must be on the left hand
            % Set default figure position to right hand monitor
            set(0,'DefaultFigurePosition',monitors(Rmonitor,1:4));
            Dmonitor=Rmonitor;
        end
    end
end
% msgbox(['not isdeployed: ' num2str(~isdeployed)]);      

clear desktop desktopMainFrame 
%------------------------------------------------------------------------------------------------------------
    
    sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';...
        'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';...
        'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'interpolated using questionable SST'};
    vartitlesc={'x-axis Time';'Date UTC';'Time UTC';'Latitude';...
        'Longitude';'EQU Pressure';'EQU Temp';'H2O flow';'gas flow';...
        'LICOR Cavity';'Ambient Pressure';'Sample Type';'Salinity';...
        'in situ Temp';'Std Value';'LICOR XCO2';'LICOR H2O';'ATM Pressure';'cpu Temp';...
        'Deck Box T';'Deck Box P';'Condenser T';'QC';'SubFlag';'User SubFlag'};
    inihead=strvcat('Default Data File Path','Default Result File Path','Default Header File Path',...                                                    
             'Default Header File','Default System Configuration File Path',...                                      
             'Default System Configuration File','Default xml Data Path');                             
    off_on=strvcat('off','on');
    
    fversionn=(strrep(mfilename,'pCO2_Sys_v','')); %#ok<*ST2NM>
    fversion=sprintf('Version\n%s.%s%s',num2str(fversionn(1)),num2str(fversionn(2:3)),fversionn(4:end));
    options.Interpreter='tex';

    vartitles=char(vartitlesc);
    titles=[];
    nvar=zeros(1,51);%Ini variables
    where_m='';where_exe='';where_xml='';
    if (isdeployed & ~ispc) | ~isdeployed  % deployed on mac or Running from matlab
        where_m=which([mfilename, '.m'], '-all');mess={'No m file found.';'Several m files found.'};
        needed_files_mac={['pCO2_Sys_v' fversionn '.m'];['pCO2_Sys_v' fversionn '.fig'];'Edit_xml.m';'Edit_xml.fig';'Savef_v2.m';'Savef_v2.fig';'merge.m';'findjobj.m';};

        if size(where_m,2)==1
            def_folder = strrep([mfilename('fullpath') '.m'],[filesep mfilename '.m'],'');
            projDir=def_folder;
            where_m = def_folder;if isdeployed, where_exe=ctfroot;end
        else
            mmm=msgbox(sprintf('pCO2_Sys m file error:\n\n%s\n\nProgram aborted.',mess{(size(where_m,2)>1)+1}),'path error','error','modal');
            uiwait(mmm);
            return;
        end
        nf=dir(def_folder);
        mf=ismember(needed_files_mac,{nf.name});
        if sum(mf)~=size(needed_files_mac,1)
            mess=sprintf('%s\n',needed_files_mac{~mf});
            mmm=msgbox(sprintf('Files missing:\n\n%s\nProgram function not guaranteed.\n',mess),'file error','error','modal');
            uiwait(mmm);
        end            
    else % deployed on pc
        [status, result] = system('set PATH');
        def_folder = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
        projDir=def_folder;where_exe=def_folder;
        nf=dir(def_folder);
    end

   %check presence of necessary files
   if ~ismember([mfilename '.ini'],{nf.name}),  create_ini_file(def_folder,[mfilename '.ini']);   end
   if ~ismember('Initial System Config.ini',{nf.name}),  create_ini_file(def_folder,'Initial System Config.ini');   end
   if ~ismember('Initial Data Config.csv',{nf.name}),  create_ini_file(def_folder,'Initial Data Config.csv');   end
   where_ini=[def_folder filesep mfilename '.ini'];
  if exist(where_ini,'file')~=2
       mmm=msgbox(sprintf('could not locate or create main ini file.\n\n%s\n Try moving the program to a folder\nwhere you have write privileges.\n\nProgram aborted',[mfilename, '.ini']),'ini error','modal');
       uiwait(mmm);
       return;
   end

%    'Initial Data Config.csv';'Initial System Config.ini';'gshhs_l.b';'gshhs_f.b';'gshhs_h.b'};
   
    %close(findobj(allchild(0), 'flat', 'FileName', [projDir filesep 'pCO2_Sys_v140.fig']));
    close(findobj(allchild(0), 'flat'));
    fig = openfig(mfilename);
 %   movegui(fig,'north');
    
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'DefaultUicontrolBackgroundColor'));

    % Generate a structure of handles to pass to callbacks, and store it.
    cmh=[];%context menu handles for variable assignment
    handles = guihandles(fig);
    guidata(fig, handles);
    set(handles.figure1,'Toolbar','figure');
    set(handles.list_text,'Visible','off','String','','TooltipString','');

    
    map_folder='';
    for  i=1:size(mapfn,1)
        where_map=which(mapfn(i,:), '-all');
        if size(where_map,1)>1
            if ~isempty(find(~cellfun(@isempty,strfind(where_map,def_folder))))
                map_folder=def_folder;
                break;
            else
                mmm=msgbox(sprintf('map file error:\n\nMultiple map files (%s) found.\n\nPut a copy in the project directory:\n\n%s\n\nProgram aborted.','gshhs_*.b',def_folder),'map error','error','modal');
                CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                return;
            end
        elseif size(where_map,1)==1
            map_folder=strrep(char(where_map),[filesep mapfn(i,:)],'');
            break;
        end
    end
    if isempty(map_folder)
        if exist('map','dir')==7
            mmm=msgbox('no gshhs_*.b Files found but Matlab map dir found instead','plotting map', 'none', 'modal');
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        else
            junk=sprintf(['NO map Files or dir found.\n\nMap cannot be plotted.\n\nRequirements are either:'...
                '\n\n   1 - gshhs_*.b files (see manual)'...
                '\n\n   2 - OR the map toolbox provided by MathWorks.'...
                '\n\nYou are missing:\n']);
            
            if isempty(dir('gshhs_*.b'))
                junk=sprintf('%s\n%s',junk,'One of the map files (gshhs_*.b)');
            end
            if exist('map','dir')~=7
                junk=sprintf('%s\n%s',junk,'map toolbox');
            end
            mmm=msgbox(junk,'Plotting Error','error','modal');
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        end
    end
    
    
    
    activity(1,handles); 
    set(handles.versiont,'String',fversion,'TooltipString',fversion);
    drawnow;
    fid=fopen(where_ini);

    if fid==-1, datap=def_folder;headp=def_folder;sysinip=def_folder;resultp=def_folder;where_xml=def_folder;
    else %reads ini file
        junk='';junkl=deblank(fgetl(fid));
        i=1;
        while ischar(junkl)
           % junk=strvcat(junk,junkl);  %#ok<REMFF1>
            junka(i)={junkl};i=i+1;
            junkl=fgetl(fid); if ischar(junkl), junk1=deblank(junkl);end
        end
        fclose(fid);
        %if problem while reading ini file, paths set to projDir and filenames to ''
        %         pb=0;
        i=find(ismember(junka, strtrim(inihead(1,:)))==1);%Default Data file path
        if ~isempty(i), datap=strtrim(char(junka{i+1}));end%Default Data file path
        i=find(ismember(junka, strtrim(inihead(2,:)))==1);%Result file path
        if ~isempty(i), resultp=strtrim(char(junka{i+1}));end%Result file path
        i=find(ismember(junka, strtrim(inihead(3,:)))==1);%Header file path
        if ~isempty(i), headp=strtrim(char(junka{i+1}));end%Header file path
        i=find(ismember(junka, strtrim(inihead(4,:)))==1);%Header File
        if ~isempty(i), hfile=strtrim(char(junka{i+1}));end%Header File
        i=find(ismember(junka, strtrim(inihead(5,:)))==1);%System conf. file path
        if ~isempty(i), sysinip=strtrim(char(junka{i+1}));end%System conf. file path
        i=find(ismember(junka, strtrim(inihead(6,:)))==1);%System conf. File
        if ~isempty(i), sysinif=strtrim(char(junka{i+1}));end%System conf. File
        i=find(ismember(junka, strtrim(inihead(7,:)))==1);%xml data path
        if ~isempty(i), where_xml=strtrim(char(junka{i+1}));end%xml data path
        pb=[size(junka,2)<5; exist(datap,'dir')~=7; exist(resultp,'dir')~=7;...
            exist([headp, filesep, hfile],'file')~=2;...
            exist([sysinip, filesep, sysinif],'file')~=2;...
            exist([where_xml, filesep, 'xml.tsv.txt'],'file')~=2];
        todo={'datap=def_folder;';'resultp=def_folder;';'headp=def_folder;hfile='''';';'sysinip=def_folder;sysinif='''';';'where_xml=def_folder;'};
        for i=2:6
            if pb(i)| pb(1),eval(todo{i-1});end
        end
    end

    dateok=0;posok=0;

    Config_Sys_Load(-1, 0, handles);
    Config_Data_Load(-1, 0, handles, 0);


    hp=[];   hp2=[];   zoomon=0;   motion=0;
    set(handles.axes1,'XColor','r','YColor','r', 'XAxisLocation','bottom', 'YAxisLocation','left','NextPlot','replace','Color','none');
   % f1=   get(handles.axes1,'ButtondownFcn');
    c=get(handles.figure1,'Color');set(handles.Graphc_s,'BackgroundColor',c);
    set(handles.axes2,'Position',get(handles.axes1,'Position'),'XColor','K','YColor','K','XAxisLocation','top','YAxisLocation','right', 'Xtick',[],...
        'NextPlot','replace','Color',c);
  %  f2=   get(handles.axes2,'ButtondownFcn');

    listg=['Graph 1';'Graph 2'];
    set(handles.popGselect, 'String',listg,'Value',1);
        
%     set(handles.Flagc12,'BackgroundColor','r'); set(handles.Flagc13,'BackgroundColor','y');set(handles.Flagc14,'BackgroundColor','w');
%     set(handles.Flagc22,'BackgroundColor','k'); set(handles.Flagc23,'BackgroundColor','b'); set(handles.Flagc24,'BackgroundColor','c');

    set(handles.rdoFlag2,'Value',1); set(handles.rdoFlag3,'Value',1); set(handles.rdoFlag4,'Value',1);
    set(handles.rdomv,'Value',0);
    
    rdos=[handles.rdoAtm handles.rdoEqu];
    set(rdos,'Value',1);
    junk=get(handles.tblstdl,'Data');
    junk(:,1)={true};
    set(handles.tblstdl,'Data',junk); 
    % Sets TooltipString of controls to their label and sets fontunits to Normalized.     
  %  Sty={'text','pushbutton','radiobutton','popupmenu','edit'};
    textb=findobj(handles.figure1,'Style','text','-or','Style','pushbutton','-or','Style','radiobutton','-or','Style','edit','-or','Style','popupmenu');
    ts=get(textb,'String');
    tts=get(textb,'TooltipString');
    for i=1:length(tts)
        if ~strcmp('popupmenu',get(textb(i), 'Style'))
            if iscell(ts{i}) & size(ts{i},1)<3, ts{i}=sprintf('%s %s',char(ts{i}(1)),char(ts{i}(2)));end
            if strcmp(tts{i},'')>0,set(textb(i),'TooltipString',ts{i});end
            %ts{i}
        end
        set(textb(i),'FontUnits','normalized');
    end
    set(handles.axes1,'FontUnits','normalized');set(handles.axes2,'FontUnits','normalized');
%    set(handles.tblrange,'FontUnits','normalized');
    %Position is x y w h
    lbp=get(handles.lb,'Position');lbm=get(handles.lbmagic,'Position');
    mid=lbp(1,1)+lbp(1,3)/2;set(handles.lbmagic,'Position',[mid-lbm(1,3)/2 lbm(1,2) lbm(1,3) lbm(1,4)]);
    uistack(handles.lbmagic,'bottom');
    clear lbp lpm mid;
    
    if ~isempty(shorelines) 
        map_loaded = 2-(length(shorelines)<50000)+(length(shorelines)>180000);handles.map_res.Value=map_loaded;
        map_res_Callback(handles.map_res, 0, handles);
    end
    %     pushb=findobj(handles.figure1,'Style','pushbutton');

    if nargout > 0
        varargout{1} = fig;
    end

    activity(0,handles);
%     set(handles.figure1,'Units', 'normalized','Position',[0 0 0.9 0.9]);
%     drawnow;
    set(handles.figure1,'Units', 'pixels','OuterPosition',monitors(Dmonitor,1:4).*[1 1 0.9 0.9]);
    set(handles.figure1,'Units', 'normalized');
    drawnow;
    movegui(handles.figure1,'center');
    set(handles.figure1,'WindowButtonMotionFcn',{@figure1_WindowButtonMotionFcn,handles});
    set(handles.figure1,'WindowButtonDownFcn',{@figure1_WindowButtonDownFcn,handles});
    set(handles.figure1,'WindowButtonUpFcn',{@figure1_WindowButtonUpFcn,handles});
    f1=   get(handles.axes1,'ButtondownFcn');
    f2=   get(handles.axes2,'ButtondownFcn');

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    try
        [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    catch

        disp(lasterr);            
        if strcmp(lasterr,'Attempt to reference field of non-structure array.')~=1
            handles=guihandles(fig);
            mmm=msgbox(sprintf('Error has occurred\n(%s)',lasterr),'Error','error','modal');
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        end
%         fig = openfig(mfilename,'reuse');
        activity(0,guihandles(fig));
        drawnow;
    end


end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and
%| sets objects' callback properties to call them through the FEVAL
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the
%| callback type separated by '_', e.g. 'map_zoom_out_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.map_zoom_out. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function  btnPlot_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.
global DataA ind1 ind2 listf Currentplot vspeed vspeedok stdv;
global Nonecol vtype  vstdx  vflag  vstdi;
global vxco2dry vlicorxcorr vxco2icorr vlicorxw vlicorxa vfco2w vfco2a vfco2i vdfco2;
% global vflagcol;


activity(1,handles); 
drawnow;

radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4];
tradioh=[handles.rdoAtm handles.rdoEqu];
typea={};
for i=1:length(stdv)
    typea=[typea;{['STD' num2str(i)]}];
end
typea=[typea;{'ATM'};{'EQU'}];
for i=1:size(handles.tblunkl.Data,1)
    typea=[typea;{['UNK' num2str(i)]}];
end
% typea={'STD1';'STD2';'STD3';'STD4';'ATM5';'EQU';...
%     'STD1-DRAIN';'STD2-DRAIN';'STD3-DRAIN';'STD4-DRAIN';'ATM5-DRAIN';'EQU-DRAIN'};
% typea={'STD1';'STD2';'STD3';'STD4';'ATM5';'EQU';...
%     'STD1-DRAIN';'STD2-DRAIN';'STD3-DRAIN';'STD4-DRAIN';'ATM5-DRAIN';'EQU-DRAIN';...
%     'STD1z';'STD2z';'STD3z';'STD4z';'NOTHING';'NOTHING';...
%     'STD1s';'STD2s';'STD3s';'STD4s';'NOTHING';'NOTHING'};

if strcmp(get(handles.RightY,'visible'),'off')==1
    set(handles.RightY,'value',Nonecol);
end

%if StdXC02, Type or not "plotable"
testh=[handles.CurrentX.Value handles.LeftY.Value handles.RightY.Value];

comph=[vstdx  vtype];
lsth=[vstdi vxco2dry vlicorxcorr vxco2icorr vlicorxw vlicorxa vfco2w vfco2a vfco2i vdfco2];

mmm=0;
if mmm~=0,return;end

for i=1:length(comph)%cannot plot what's in comph
    if ~isempty(find(testh==comph(i),1)) & comph(i)~= Nonecol  
        mess=sprintf('\"%s\" cannot be plotted',char(listf(comph(i))));
        mmm=msgbox(mess,'Plotting Error','error','modal');
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    end
    if mmm~=0,break;end
end
if mmm~=0,return;end
% end
if ~isempty(find(testh==vspeed,1)) & ~vspeedok  
    mess=sprintf('Speed not calculated.\nMake sure YDay is OK and Re-calculate Fields');
    mmm=msgbox(mess,'Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
junk=get(handles.tblstdl,'Data');stdl=junk(:,1);
tradiov=[handles.tblstdl.Data{:,1} tradioh.Value]; if size(handles.tblunkl.Data,2)>0,tradiov=[tradiov handles.tblunkl.Data{:,1}];end
stdn=sum([handles.tblstdl.Data{:,1}]);
if (sum(tradiov)==0)
    mmm=msgbox('No Data selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

if strcmp(listf(handles.LeftY.Value),'none')==0
    ind11={[];[];[]};
    for i=1:length(tradiov)  %Select Type
        if tradiov(i)
            if (i>length(stdl)) | (isempty(find(lsth==handles.LeftY.Value,1))) %type not STD OR leftY can be plotted for the type "std"
                for j=1:3
                    if radio(j).Value
                        if verLessThan('matlab','9.1'),ind11{j}=[ind11{j};find(~cellfun(@isempty,strfind(DataA(:,vtype),typea{i})) & strcmp(DataA(:,vflag),num2str(j+1))==1)];
                        else, ind11{j}=[ind11{j};find(contains(DataA(:,vtype),typea{i}) & strcmp(DataA(:,vflag),num2str(j+1))==1)];end
                            ind11{j}=sort(cat(1,ind11{j}));    
                    else
                        ind11{j}=[ind11{j};[]];
                    end
                end
            end
        end
    end
   % ind11=sort(cat(1,indf{:}));
else
    ind11{1}=[]; ind11{2}=[]; ind11{3}=[];
end % leftY

if strcmp(listf(handles.RightY.Value),'none')==0 & strcmp(get(handles.RightY,'visible'),'off')==0
    ind22={[];[];[]};
    for i=1:length(tradiov)  %Select Type
        if tradiov(i)
            if (i>length(stdl)) | (isempty(find(lsth==handles.RightY.Value,1))) %type not STD OR RightY can be plotted for the type "std"
                for j=1:3
                    if radio(j).Value
                        if verLessThan('matlab','9.1'),ind22{j}=[ind22{j};find(~cellfun(@isempty,strfind(DataA(:,vtype),typea{i})) & strcmp(DataA(:,vflag),num2str(j+1))==1)];
                        else, ind22{j}=[ind22{j};find(contains(DataA(:,vtype),typea{i}) & strcmp(DataA(:,vflag),num2str(j+1))==1)];end
                            ind22{j}=sort(cat(1,ind22{j}));    
                    else
                        ind22{j}=[ind22{j};[]];
                    end
                end
            end
        end
    end
    %ind22=sort(cat(1,indf{:}));
else
    ind22{1}=[]; ind22{2}=[]; ind22{3}=[];
end % RightY

if double(get(handles.rdomv,'Value'))==0  %don't plot missing values
    indmv=find(strcmp(DataA(:,handles.CurrentX.Value),'-999')==0 & strcmp(DataA(:,handles.LeftY.Value),'-999')==0);
    for i=1:3     ind11{i}=intersect(indmv,ind11{i}); end
    indmv=find(strcmp(DataA(:,handles.CurrentX.Value),'-999')==0 & strcmp(DataA(:,handles.RightY.Value),'-999')==0);
    for i=1:3     ind22{i}=intersect(indmv,ind22{i}); end
end


mmm=0;
if strcmp(listf(handles.LeftY.Value),'none')==0 & strcmp(listf(handles.RightY.Value),'none')==0
    listg=[(DataA(1,handles.LeftY.Value));(DataA(1,handles.RightY.Value))];
%     set(handles.popGselect, 'String',listg,'Value',2);
elseif strcmp(listf(handles.LeftY.Value),'none')==0
    listg=[(DataA(1,handles.LeftY.Value));'Graph 2'];
    set(handles.popGselect, 'String',listg,'Value',1);
elseif strcmp(listf(handles.RightY.Value),'none')==0
    listg=['Graph 1';(DataA(1,handles.RightY.Value))];
    set(handles.popGselect, 'String',listg,'Value',2);
else
    listg=['Graph 1';'Graph 2'];
%     set(handles.popGselect, 'String',listg,'Value',1);
    set(handles.popGselect, 'String',listg);
    mmm=msgbox('No Data Selected','Plotting Error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
end
if mmm~=0, return;end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@btnPlot_Callback;
clear ind11 ind22 indmv indf;
plotdata_double(handles,ind1, ind2, handles.CurrentX.Value, handles.LeftY.Value, handles.RightY.Value,h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end
activity(0,handles); 


% --------------------------------------------------------------------
function Config_Data_Load(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.loadhc.
global datap resultp  headp hfile sysinip sysinif projDir def_folder where_xml;
global titles vartitles dateok;


if dateok==1
    mmm=questdlg('This will reset the data processing.','Data Config Load','Reset','Cancel','Cancel');
    if strcmp(mmm,'Cancel')==1, return; end
end


if isobject(h) , activity(1,handles); end
drawnow;

if exist([headp,filesep,hfile],'file')~=2,  headp=[def_folder , filesep, 'System Configurations'];end
if exist([headp,filesep,hfile],'file')~=2,  hfile='Initial Data Config.csv';headp=def_folder;   end
if exist([headp,filesep,hfile],'file')~=2,  if sum(ismember('WIN',computer))==3, headp='c:'; elseif sum(ismember('MACI',computer))==4, headp=''; end,    end

if isobject(h) | exist([headp filesep hfile],'file')~=2
%     ccd=pwd;cd(headp);
    mmm=msgbox(sprintf('%s\n\n%s\n\n','Select a DATA configuration file','in the next popup window!'),'Data Config','modal');
    CenterWindow(handles,mmm,handles.figure1); uiwait(mmm);
    [fname,fpath] = uigetfile('*.csv','Select Data Config File',[headp filesep hfile]);
    if (fpath==0), return; end
    hfile=fname;
    if strcmp(fpath(end),filesep)>0,fpath=fpath(1:end-1);end
    headp=fpath;
    set(handles.loadhct,'String',hfile,'TooltipString',hfile);
end


if exist([headp filesep hfile],'file')
    titles={0};nvartitles={0};%initializes titles as cell

    fid=fopen([headp filesep hfile]);

    line1=fgetl(fid);  lines=fgetl(fid);
    if isempty(line1) | isempty(lines),  return;end;
    
    %for compatibility with previous versions with licors 840,6262,7000
    line1=strrep(line1,'LICOR Temp','LICOR Cavity');
    line1=strrep(line1,'LICOR Pressure','Ambient Pressure');
    
    nvartitles = textscan(line1,'%s','delimiter',',');
    titles = textscan(lines,'%s','delimiter',',');

    titles=strtrim(char(titles{:,1}));  nvartitles=strtrim(char(nvartitles{:,1}));
    nvartitles=strvcat(nvartitles,vartitles(end-2,:));%#ok<VCAT> % adds 'QC'
    nvartitles=strvcat(nvartitles,vartitles(end-1,:));%#ok<VCAT> % adds 'SubFlag'
    nvartitles=strvcat(nvartitles,vartitles(end,:));%#ok<VCAT> % adds 'User SubFlag'

    fclose(fid);

    if strcmp(nvartitles,vartitles)~=1
        hh=msgbox(sprintf('First line different than expected.\n%s possibly compromised.\nNeed new configuration.',hfile),'Header Config Error','error','modal');
        CenterWindow(handles,hh,handles.figure1)
        uiwait(hh);
    end
    set(handles.loadhct,'String',hfile,'TooltipString',hfile); 
    Save_Ini(datap,resultp,headp,hfile,sysinip,sysinif,where_xml);
else
    set(handles.loadhct,'String','nothing.yet','TooltipString','nothing.yet');
end

if isobject(h), Reset_Reduction(-1, 0, handles);end

if isobject(h) , activity(0,handles); end


% --------------------------------------------------------------------
function Reset_Reduction(hObject, eventdata, handles)
    
    global DataA f1 f2 listf listf_file listf_calc npts titles nvar ncol dateok ;
    global datap  headp hfile sysinip sysinif ;
    global Nonecol vtype hp hp2 ;
    global stdok atmok airiok vfco2ok xco2corrok vspeedok;
    global stdv ;
    global vartitles motion off_on   nrec ;

    


    try

        %save conditions in case of fail
        DataO=DataA; stdok0=stdok;atmok0=atmok;vfco2ok0=vfco2ok; airiok0=airiok;vspeedok0=vspeedok;xco2corrok0=xco2corrok;motion0=motion;
        datap0=datap;headp0=headp;hfile0=hfile;sysinip0=sysinip;sysinif0=sysinif;

   
        stdok=0;atmok=0;vfco2ok=0; airiok=0;vspeedok=0;xco2corrok=0;motion=0;dateok=0;clearvars -global rfh;nrec=0;
        off_list=[handles.air_interpolation handles.latlong_interpolation handles.btnCalcYDay...
            handles.correctxco2   handles.cfco2  handles.btnSave handles.btnXML  handles.SaveConfig handles.btnplot...
            handles.OSet handles.OApply handles.OCheck];

        off_r_list=[handles.ydayok  handles.posok  handles.stdok  handles.xco2aok];
        fields_list=[handles.ssti_f handles.stdo_f handles.dt_f handles.dp_f handles.xco2d_f handles.sp_f];

        set(off_list,'Visible','off');
        set(off_r_list,'Visible','off','BackgroundColor','r');
        set(fields_list,'TooltipString','','BackgroundColor','r');
        %if  get(handles.rdodryxco2,'value')==0 cc=get(handles.figure1,'Color');set(handles.xco2d_f,'TooltipString','xCO2 dry N/A','BackgroundColor',cc);end
        set([handles.xco2d_f,handles.xco2d_ft],'Visible',off_on(get(handles.rdodryxco2,'value')+1,:));
        % set(handles.xco2d_ft,'Visible',off_on(get(handles.rdodryxco2,'value')+1,:));

        set(handles.ydayokt,'String','Y/C Day - Not OK');
        set(handles.posokt,'String','Positions - Not Flagged');
        set(handles.stdokt,'String','Standards - Not Flagged');
        set(handles.xco2aokt,'String','xCO2(atm) - Not Flagged');

        clear off_list off_r_list fields_list;


        if ~isempty(listf_file) % if data files are in memory already

            if vtype~= Nonecol
                tblunk_modify(handles.tblunkmenu.Children(1), eventdata, handles) ;%adjust number of UNK in table from file info
            end

            %Reset calculated Fields headers and data
            stdtit='';for i=1:length(stdv),stdtit{1,i}=['Std' num2str(i) ' values'];end

            fieldtitles=[{'Speed(knots)'},{'Tin interp'},{'Teq-Tin'},{'Peq-Patm'},{'EQU Pressure (calc)'},{'Std Offset'},...
                stdtit,...
                {'xCO2 dry'},{'xCO2 corr'},{'xCO2 interp corr'},{'xCO2W Corrected'},{'xCO2A Corrected'},{...
                'fCO2 water'},{'fCO2 air'},{'fCO2 interp'},{'dfCO2(w-a)'}];
            qctitles=[{'QC'},{'SubFlag'},{'SubFlag_User'},{'YDay Calc'},{'Cruise Day'},{'none'}];

            %***** Data Config File Use *************
            %Assigns variable from titles to column number in file (listf_file) according to config file
            nvar=zeros(1,51);%Ini variables
            if hObject~=-1,Config_Data_Load(-1, 0, handles, 0);end % if not from Data_Config_Load
            for i=1:size(titles,1)%Ini variables
                nv=find(strcmpi(strtrim(listf_file),strtrim(titles(i,:)))==1);
                if length(nv)>0,   nvar(i)=nv(1);    end
            end

            ncol=size(listf_file,2);
            numvar=length(vartitles);%these are the variables taken from the file + flag, subflag, subflaguser (See vartitles)
            %j=numvar+1;%+1 to leave a space in nvar for vsubfu
            j=numvar;
            DataA=strtrim(DataA);
            DataA(:,ncol+1:end)=[];
            %reset calc fields
            for i=1:length(fieldtitles) %19
                nvar(j+i)=ncol+i;
                DataA(1,ncol+i)=fieldtitles(i);
                DataA(2:npts,ncol+i)={'-999'};
            end
            i=i+1;
            DataA(1,ncol+i)={'QC'}; DataA(2:npts,ncol+i)={'2'}; nvar(numvar-2)=ncol+i; i=i+1;%Sets Flags to 2
            DataA(1,ncol+i)={'SubFlag'}; DataA(2:npts,ncol+i)={''}; nvar(numvar-1)=ncol+i; i=i+1;%Sets SubFlags to ''
            DataA(1,ncol+i)={'SubFlag_User'}; DataA(2:npts,ncol+i)={''}; nvar(numvar)=ncol+i; i=i+1;%Sets SubFlags User to ''
            DataA(1,ncol+i)={'YDay Calc'}; DataA(2:npts,ncol+i)={'-999'}; nvar(numvar+length(fieldtitles)+1)=ncol+i; i=i+1;%Sets YDay Calc column to -999
            DataA(1,ncol+i)={'Cruise Day'}; DataA(2:npts,ncol+i)={'-999'}; nvar(numvar+length(fieldtitles)+2)=ncol+i; i=i+1;%Sets Cruise Day column to -999
            DataA(1,ncol+i)={'none'}; Nonecol=ncol+i; DataA(2:npts,Nonecol)={'-999'}; %Sets "None" column to -999

            nvar(nvar==0)=Nonecol; %replaces 0 by Nonecol


            listf=DataA(1,1:end);
            listf=strtrim(listf);
            listf_calc=DataA(1,ncol+1:end-1);
            listf_calc=strtrim(listf_calc);

            Update_List(handles);
            Update_Var(handles);

            %Check #STD and UNK in file and advise
            mmm2=Check_nSTD_from_file;
            if strcmp(mmm2,'1 - Program fix'),  tblstd_modify(handles.tblstdmenu.Children(1), eventdata, handles) ;end


            %Check for duplicate columns
            ss='';samei_all=[];samei_dup=[];
            for i=1:size(DataA,2)
                samei=find(strcmpi(DataA(1,i),DataA(1,:))==1);
                if size(samei,2)>1
                    if max(samei)>ncol %one duplicate is a calculated column - should be last one (not renumbered)
                        for j=1:size(samei,2)-1,   DataA(1,samei(j)) = {[DataA{1,samei(j)} '_' strtrim(num2str(j+1))]};  end  %1st col renumbered starting at _2
                        ss='It could conflict with a column created by the program.';
                    else
                        for j=2:size(samei,2),  DataA(1,samei(j)) = {[DataA{1,samei(j)} '_' strtrim(num2str(j))]};  end   %2nd col renumbered starting at _2
                    end
                    samei_dup=unique([samei_dup i]);
                end
            end
            if ~isempty(samei_dup)
                dupcol = cellfun(@(x) [x],DataA(1,samei_dup),'UniformOutput',false);
                dupcoltxt=sprintf('%s\n',dupcol{:});
                if ~isempty(ss)
                    mmm=msgbox(sprintf('The following are duplicate header(s):\n\n%s\n%s\n\nThey have been modified\n\nto differentiate them.',dupcoltxt,ss),'Data Import...','Warn','modal');
                else
                    mmm=msgbox(sprintf('The following are duplicate header(s):\n\n%s\nThey have been modified\n\nto differentiate them.',dupcoltxt),'Data Import...','Warn','modal');
                end
                CenterWindow(handles,mmm,handles.figure1)
                uiwait(mmm);
                listf=DataA(1,1:end);
                listf=strtrim(listf);
                listf_file=DataA(1,1:ncol);
                listf_calc=DataA(1,ncol+1:end-1);
                listf_calc=strtrim(listf_calc);

                Update_List(handles);
                Update_Var(handles);
            end

            Calc_Fields_Callback(handles);
            activity(1,handles);%Calc_Fields above turns it off...
        end
        cla(handles.axes1);cla(handles.axes2);
        c=get(handles.Graphc_s,'BackgroundColor');
        set(handles.axes1,'XColor','r','YColor','r', 'XAxisLocation','bottom', 'YAxisLocation','left', 'Color','none');
        set(handles.axes2,'Position',get(handles.axes1,'Position'), 'XColor','k','YColor','k', 'XAxisLocation','top', 'YAxisLocation','right', 'Xtick',[],'Color',c);
        set(handles.axes1,'ButtonDownFcn',f1);set(handles.axes2,'ButtonDownFcn',f2);
        hp=[];hp2=[];
        set(handles.axes2,'XLim',[0 1],'YLim', [0 1]);

    
        on_list=[handles.merge_butt  handles.delete_col_butt  handles.btnCalcYDay  handles.btnSave  handles.btnXML...
            handles.SaveConfig  handles.lr handles.btnplot handles.OSet handles.OApply handles.OCheck];
        set(on_list,'Visible','on');drawnow;

        % cd(ccd);
        Check_for_999(handles);

        activity(0,handles);

        clear mess mess1 mess2 iData
        activity(1,handles);
        drawnow;

    catch

        mmm=msgbox(lasterr,'Error');
        CenterWindow(handles,mmm,handles.figure1)
        uiwait(mmm);
        %restore conditions in case of fail
        if exist('DataO','var'), DataA=DataO; end
        stdok=stdok0;atmok=atmok0;vfco2ok=vfco2ok0; airiok=airiok0;vspeedok=vspeedok0;xco2corrok=xco2corrok0;motion=motion0;
        datap=datap0;headp=headp0;hfile=hfile0;sysinip=sysinip0;sysinif=sysinif0;

        on_list=[handles.merge_butt  handles.delete_col_butt  handles.btnCalcYDay  handles.btnSave handles.btnXML...
            handles.SaveConfig  handles.lr handles.btnplot handles.OSet handles.OApply handles.OCheck];
        set(on_list,'Visible','on');drawnow;
        activity(0,handles);

        clear  DataO stdok0 vfco2ok0 airiok0 vspeedok0 xco2corrok0 motion0 datap0;
        clear headp0 hfile0 sysinip0 sysinif0;

    end



% --------------------------------------------------------------------
function Data_Load(hObject, eventdata, handles)
% Stub for Callback of the uicontrol handles.loadf.
global DataA f1 f2 listf listf_file listf_calc npts titles nvar ncol dateok clear_savew;
global datap resultp headp hfile sysinip sysinif where_xml expot dryxco2;
global vflag  Nonecol vtype hp hp2;
global vxaxis vyday vydayi vcrday vlat vlong vpeq vteq vwflo vgflo  vlicorcav  vpamb  vsal  vtin  vlicorx  vlicorw  vpatm vpdbox;
global vxaxiscok posok stdok atmok airiok vfco2ok xco2corrok vspeedok;
global stdv usestdv use0 GSID baroh equpdiff fversionn fversions vgroup vship pi_names vcruiseid save_opts;
global vartitles motion off_on drag rfh nrec v_compatible hdllohistr hdltype;

try

    %save conditions in case of fail
    DataO=DataA; stdok0=stdok;atmok0=atmok;vfco2ok0=vfco2ok; airiok0=airiok;vspeedok0=vspeedok;xco2corrok0=xco2corrok;motion0=motion;
    datap0=datap;headp0=headp;hfile0=hfile;sysinip0=sysinip;sysinif0=sysinif;
    
    if exist(datap,'dir')~=7, if sum(ismember('WIN',computer))==3, datap='c:'; elseif sum(ismember('MACI',computer))==4, datap=''; end,    end
    v  = version ('-release'); 
    if str2num(v(1:4))<2013, lf=['*dat.txt; *.csv; *.xls; *.mat'];
    else lf=['*dat.txt; *.csv; *.xls; *.xlsx; *.mat'];end
    [fname,fpath] = uigetfile(lf,'Select Data File(s)',[datap filesep],'MultiSelect','on');
    if (fpath==0),return;end
    if strcmp(fpath(end),filesep)>0,fpath=fpath(1:end-1);end
    datap=fpath; %addpath (datap); cd(datap);
    if exist(datap,'dir')==7,    Save_Ini(datap,resultp,headp,hfile,sysinip,sysinif,where_xml); end
    
    
    fname2='';

    DataA={0}; stdok=0;atmok=0;vfco2ok=0; airiok=0;vspeedok=0;xco2corrok=0;motion=0;dateok=0;clearvars -global rfh;nrec=0;
    off_list=[handles.air_interpolation handles.latlong_interpolation handles.btnCalcYDay ...
        handles.correctxco2   handles.cfco2  handles.btnSave handles.btnXML  handles.SaveConfig handles.btnplot...
        handles.OSet handles.OApply handles.OCheck];
    
    off_r_list=[handles.ydayok  handles.posok  handles.stdok  handles.xco2aok];
    fields_list=[handles.ssti_f handles.stdo_f handles.dt_f handles.dp_f handles.xco2d_f handles.sp_f];
    
    set(off_list,'Visible','off');
    set(off_r_list,'Visible','off','BackgroundColor','r');
    set(fields_list,'TooltipString','','BackgroundColor','r');
    %if  get(handles.rdodryxco2,'value')==0 cc=get(handles.figure1,'Color');set(handles.xco2d_f,'TooltipString','xCO2 dry N/A','BackgroundColor',cc);end
    set([handles.xco2d_f,handles.xco2d_ft],'Visible',off_on(get(handles.rdodryxco2,'value')+1,:));
   % set(handles.xco2d_ft,'Visible',off_on(get(handles.rdodryxco2,'value')+1,:));

    set(handles.ydayokt,'String','Y/C Day - Not OK');
    set(handles.posokt,'String','Positions - Not Flagged');
    set(handles.stdokt,'String','Standards - Not Flagged');
    set(handles.xco2aokt,'String','xCO2(atm) - Not Flagged');
    
    clear off_list off_r_list fields_list;
    
    activity(1,handles);
    drawnow;


    nl=0;lname=1;fnameo=fname;
    
    if iscell(fname), cond=~isempty(cell2mat(strfind(fname,'dat.txt')));cond1=~isempty(cell2mat(strfind(fname,'.csv')));
    else cond=~isempty(findstr(fname,'dat.txt'));cond1=~isempty(findstr(fname,'.csv')); end
    
    % ************************************** dat.txt or csv files*****************************************************************
    
    if cond | cond1 %data file is a dat.txt file. Gets rid of shutdown lines and imports data in iData
        if iscell(fname)
            dattxtnum=size(cell2mat(strfind(fname,'dat.txt')),2);
            if dattxtnum~=size(fname,2) & dattxtnum>0
                mmm=msgbox(sprintf(['Several files have been selected.\n\nOnly ' char(39)  '*dat.txt' char(39) ' files will be loaded.']),'File Error','error','modal');
                CenterWindow(handles,mmm,handles.figure1)
                uiwait(mmm);
                %return;
            end
        end
        
        if ischar(fname)% 1 file selected
            lname=1;
        else
            lname=length(fname);     fname=sort(fname);
        end
        
        fname2='temp_dat.txt';    [fid2]=fopen([datap filesep fname2],'w');
        
        j=0;delines=zeros(lname,1);k=0;
        for i=1:lname
            kk=0;%nb of lines loaded per file
            if lname>1
                set(handles.lr,'Visible','on');drawnow;
                if ~isempty(strfind(fname{i},'dat.txt'))
                    fid= fopen([datap filesep fname{i}],'r');k=k+1;
                    set(handles.lr,'String',sprintf('loading file... %d/%d', k, dattxtnum));
                    drawnow;
                else continue;end
            else
                fid=fopen([datap filesep fname],'r');
                set(handles.lr,'String','loading file... 1/1');drawnow;
            end
            
            lines=fgetl(fid); %reads header of 1st file
            if  i==1 & j==0 %prints header of 1st file
                if cond1 %csv file
                    ncol0=length(regexp(lines,','));
                    lines=sprintf('%s, o, 1',lines);%TRICK - adds 2 extra text column to import all in 'textdata'
                    lines=regexprep(lines,',','\t');
                else
                    ncol0=length(regexp(lines,'\t'));
                    lines=sprintf('%s\to\t1',lines);%TRICK - adds 2 extra text column to import all in 'textdata'
                end
                fprintf(fid2,'%s\n',lines);
            end
            while ~feof(fid)%Input of data
                lines=fgetl(fid);j=j+1;nl=nl+1;
                if contains(lines,'06/07/24')>0 & contains(lines,'19:31:44')>0
                    lines
                end
                if cond1 %csv file
                    lines=sprintf('%s, o, 1',lines);%TRICK - adds 2 extra text column to import all in 'textdata'
                    %Checks for missing values and replaces them by -9: missing values detected by 2 commas next to each other or separated by one space
                    lines=strrep(strrep(strrep(lines,',,',',-9,'),', ,',',-9,'),',,',',');% the 3rd strrep is to get rid of ",," left by strrep when holes consecutive
                    lines=regexprep(lines,',','\t');
                    ncol1=length(regexp(lines,','));
                    if ncol1>ncol0+2  %2 extra column added earlier
                        if lname>1 fnn=fname{i};else fnn=fname;end
                        mmm=msgbox(sprintf('Input ABORTED. Fix and restart.\n\nToo many columns in this line of file: %s.\n\n%s',fnn,lines),'File Error','error','modal');
                        CenterWindow(handles,mmm,handles.figure1)
                        uiwait(mmm);
                        set(handles.lr,'Visible','off');drawnow;
                        return;
                    end
                else
                    lines=sprintf('%s\to\t1',lines);%TRICK - adds 2 extra text column to import all in 'textdata'
                    %fixed common bad strings for some ships (this section
                    %can later be pulling regexp from a file...
                    %EQ - search for similar:1*60	$WIMWV	C
                    lines=regexprep(lines,'\<[a-zA-Z0-9]\*\w*\t\<\$\w*\t\w*','-9\t-9');%replaces "any word start with 1* tab any word start with $ tab any word tab"
                    lines=regexprep(lines,'[a-zA-Z]\t[a-zA-Z]','-9\t-9');%replaces "any 2 letters following each other like H R"
                    lines=regexprep(lines,'\t[a-np-zA-Z]\t','\t-9\t');%replaces "any single letter  like M (except o) - see EQNX sometimes"
                    
                    %IslanderII - search for number with 2 dots (in SST)
                    lines=regexprep(lines,'\t[0-9]*\.[0-9]*\.[0-9]*','\t-9');%replaces "any number with 2 dots - see ISL2 sometimes"
                     %Armstrong - remove H when atm P is less than 1000
                     lines=regexprep(lines,'\t(\d\d\d\.\d)H','\t$1');%replaces "997.5H" with "997.5"
                    
                    %Checks for missing values at end of line and replaces them by -9
                    %must be done before checking for other missing values
                    %may not be needed anymore now that \to\t1 is added before the checks so no missing data at end of line anymore
                    % ncol1=length(regexp(lines,'\t'));
                    % if ncol1<ncol0 %assumes last columns are missing data so adds \t-9 for each
                    %     lines=sprintf(['%s' repmat('\t-9',1,ncol0-ncol1)],lines);
                    % end
                    % %another way to check for missing values at end of line
                    % lines=sprintf('%s%s',lines(1:end-1),regexprep(lines(end),'\t','\t-9'));% if missing data is last on line

                    %Checks for missing values and replaces them by -9
                    lines=strrep(strrep(lines,'\t\t','\t-9\t'),'\t\t','\t');% the 2nd strrep is to get rid of "\t\t" left by strrep when holes consecutive
                    lines=regexprep(regexprep(lines,'\t\t','\t-9\t'),'\t\t','\t-9\t'); % trying regexprep in case strrep didn't work (happened on some files)(regexprep workks differently)
                    lines=regexprep(lines,'\t[^0-9\.\\\+:]+\t','\t-9\t');%replaces "any character different than [0-9.\:+] between 2 tabs by -9
                    ncol1=length(regexp(lines,'\t'));
                    if ncol1>ncol0+2  %2 extra column added earlier
                        if lname>1 fnn=fname{i};else fnn=fname;end
                        mmm=msgbox(sprintf('Input ABORTED. Fix and restart.\n\nToo many columns in this line of file: %s.\n\n%s',fnn,lines),'File Error','error','modal');
                        CenterWindow(handles,mmm,handles.figure1)
                        uiwait(mmm);
                        set(handles.lr,'Visible','off');drawnow;
                        return;
                    end
                end
                mot=regexp(lines,{'FILTER' 'SHUT DOWN' 'SLEEP' 'WAKE UP' 'EMERGENCY STOP'});
                if isempty([mot{:}]) & length(lines)>0 
                    type=regexp(lines,{'STD' 'UNK' 'EQU' 'ATM'});
                    if isempty([type{:}]) & length(lines)>0 
                        if lname>1 fnn=fname{i};else fnn=fname;end
                        mmm=msgbox(sprintf('Line skipped. Check file.\n\nNo ''Type'' in this line of file: %s.\n\n%s',fnn,lines),'File Error','error','modal');
                        CenterWindow(handles,mmm,handles.figure1)
                        uiwait(mmm);
                        nl=nl-1;
                        delines(i)=delines(i)+1; 
                    else
                        fprintf(fid2,'%s\n',lines);kk=kk+1;
                        lastline=lines; %for debugging
                    end
                else
                    nl=nl-1;
                    if isempty([mot{:}]) & ncol0~=ncol1, delines(i)=delines(i)+1; end
                end
            end
            fclose(fid);
            if kk==0
                mmm=msgbox(sprintf(['No data imported from ' char(39)  fname{i} char(39) '.\n\nCheck file or its number of columns.']),'File Error','error','modal');
                CenterWindow(handles,mmm,handles.figure1)
                uiwait(mmm);
            end

        end
        
        fclose(fid2);
        iData=importdata([datap filesep fname2],'\t');
        if size(iData.textdata,2)<5  %less than 5 columns imported
                mmm=msgbox(sprintf('Input ABORTED. Fix and restart.\n\nLess than 5 columns imported.\n\nCheck data files.'),'File Error','error','modal');
                CenterWindow(handles,mmm,handles.figure1)
                uiwait(mmm);
                set(handles.lr,'Visible','off');drawnow;
                return;
        end
        % Create new variables in the base workspace from those fields.
        vars = fieldnames(iData);
        for i = 1:length(vars),  assignin('base', vars{i}, iData.(vars{i}));  end
        iData.textdata=iData.textdata(:,1:end -1);
        delete ([datap filesep fname2]);
        if lname==1 | dattxtnum==1
            set(handles.lr,'String',sprintf('file : %s', fname ));
        else
            set(handles.lr,'String',sprintf('%d files loaded', dattxtnum ));
            fname=fname2;
        end
        iData.textdata(find(strcmp(iData.textdata,'-9')==1))={'-999'};
        iData.textdata=strrep(iData.textdata,'NaN','-999');
        DataA=iData.textdata;
        DataA=regexprep(DataA,['[' char(0) '-' char(31) ']'],char(32));% Replaces invisible characters (ascii<32) by space before stripping
        DataA=strip(DataA);
        
        listf_file=DataA(1,:);
        listf_file=strtrim(listf_file);
        
        % ************************************** XLS files *****************************************************************
        
    elseif ~isempty(strfind(fname,'.xls')) | ~isempty(strfind(fname,'.xlsx'))
        [~,sheets] = xlsfinfo([fpath  filesep fname]);
        xlnum=1;
        if numel(sheets)>1, xlnum=listdlg('ListString',sheets, 'SelectionMode', 'single','Name','Excel Data Import','PromptString', 'Select Sheet to Import.');end
        if (isempty(xlnum) | xlnum>numel(sheets)),return;end
        %iData=xlsread([fpath filesep fname],sheets{xlnum});
        [~,~,DataA]=xlsread([fpath filesep fname],sheets{xlnum},'','basic'); %DataA contains all, numbers AND text
        %DataA=readcell([fpath filesep fname]);
        
        nh=1;
        for i=size(DataA,2):-1:1
            if isnan(DataA{1,i}) & i==size(DataA,2) %Only delete last empty columns.
                if sum(cellfun(@isnan,DataA(:,i),'UniformOutput',true))==size(DataA,1),DataA(:,i)=[];end % the whole col (incl. header) is NaN
            else 
                if isnan(DataA{1,i}), DataA{1,i}={['No_Header_' strtrim(numstr(nh))]}; nh=nh+1; end
                if sum(cellfun(@isnumeric,DataA(2:end,i)))>0,DataA(2:end,i)=cellstr(cellfun(@num2str,DataA(2:end,i),'uniformoutput',false));end
            end
        end
      
        DataA=regexprep(DataA,['[' char(0) '-' char(31) ']'],char(32));% Replaces invisible characters (ascii<32) by space before stripping
        DataA=strip(DataA);
        
        listf_file=DataA(1,:);
        listf_file=strtrim(listf_file);
        
        % ************************************** mat files *****************************************************************
        
    elseif ~isempty(strfind(fname,'.mat'))
        
        lname=1;    load([datap filesep fname]);
        if exist('DataO','var'), DataA=DataO;clear DataO;end
        DataA(find(strcmp(DataA,'-9')==1))={'-999'};
        DataA=regexprep(DataA,['[' char(0) '-' char(31) ']'],char(32));% Replaces invisible characters (ascii<32) by space before stripping
        DataA=strip(DataA);
        listf=DataA(1,:);    listf=strtrim(listf);
        if isnumeric(fversions), fversions=num2str(fversions),end
        compatibleok=0;
        for i=1:size(v_compatible,2)
            if ~isempty(strfind(v_compatible{i},fversionn)) && ~isempty(strfind(v_compatible{i},fversions)),compatibleok=1; end
        end
        if compatibleok==0 & ~isempty(fversions)
            mmm=msgbox(sprintf('This File Was Created\nBy a Different Version [v%1.2f]\nThan the Current One [v%1.2f].\nThe Compatibility Is Not Guaranteed',str2num(fversions)/100, str2num(fversionn)/100),'NOTE','warn','modal');
            CenterWindow(handles,mmm,handles.figure1)
            uiwait(mmm);
        end
        
        if exist('stdval2') 
            if ~isempty(stdval2)
                stdv=reshape(stdval2,[length(stdval2),1]);
            elseif exist('stds')
                stdv=str2num(char(stds));
            end
            clear stdval2;
        end %old mat file
 
        Update_List(handles);
        Update_Var(handles);
        set(handles.ydayok,'visible','on'); %   set(handles.posok,'visible','on');
        
        if ~exist('dryxco2','var'),dryxco2=0;end
        if ~exist('save_opts','var'),save_opts=[1,1,0,0,0];end
        if size(save_opts,1)>1,save_opts=reshape(save_opts,[1,length(save_opts)]);end
        if length(save_opts)<5, while length(save_opts)<5 , save_opts=[save_opts ,0]; end, end
        if isempty(pi_names), pi_names='';end
        set(handles.rdodryxco2,'value',dryxco2);  set([handles.xco2d_f,handles.xco2d_ft],'Visible',off_on(dryxco2+1,:));
        %ctrla=[handles.std1v,handles.std2v,handles.std3v,handles.std4v];
        if ~isempty(expot), set(handles.expot,'string',expot,'TooltipString',expot);end
        %if ~isempty(stds),for i=1:length(ctrla), set(ctrla(i),'String',strtrim(char(stds(i,:))),'TooltipString',strtrim(char(stds(i,:))));end, end
        if ~isempty(stdv)
                    set(handles.tblstdv,'Data',stdv,'ColumnFormat',({'bank'}));
                    tts=''; for i=1:length(stdv), tts= [tts 'STD ' num2str(i) sprintf(' : %0.2f ppm',stdv(i))];if i<length(stdv), tts=[tts sprintf('\n')]; end, end
                    set(handles.tblstdv,'TooltipString',tts);
                    stdl={};
                    for i=1:length(stdv)  stdl(i,:)={true,['STD ' num2str(i,0)]};end
                    tts=''; for i=1:length(stdv), tts= [tts 'STD ' num2str(i)];if i<length(stdv), tts=[tts sprintf('\n')]; end, end
                    set(handles.tblstdl,'Data',stdl,'TooltipString',tts);
        end
        if ~isempty(equpdiff),set(handles.rdoequpdiff,'value',equpdiff);end
        if ~isempty(usestdv),set(handles.rdostdv,'value',usestdv);end
        if ~isempty(use0),set(handles.rdoUse0,'value',use0);end
        if ~isempty(baroh),set(handles.dbaroh,'string',baroh,'TooltipString',baroh);end
        if ~isempty(GSID),vgroup=GSID(1,:);vship=GSID(2,:);vcruiseid=GSID(3,:);end
        
        if str2num(fversionn)<140
            ctrla=[handles.gflmin,handles.gflmax,handles.wflmin,handles.wflmax,handles.sstmin,handles.sstmax,handles.sssmin,handles.sssmax,...
                handles.condsrmin,handles.condsrmax,handles.dtmin,handles.dtmax];
            if ~isempty(SysC1),for i=1:length(ctrla), set(ctrla(i),'String',strtrim(char(SysC1(i,:))),'TooltipString',strtrim(char(SysC1(i,:))));end, end
            ctrla=[handles.peqmin,handles.peqmax,handles.teqmin,handles.teqmax,handles.plimin,handles.plimax,handles.tlimin,handles.tlimax,...
                handles.pdeckmin,handles.pdeckmax,handles.tdeckmin,handles.tdeckmax,];
            if ~isempty(SysC2),for i=1:length(ctrla), set(ctrla(i),'String',strtrim(char(SysC2(i,:))),'TooltipString',strtrim(char(SysC2(i,:))));end, end
            ctrla=[handles.sstdpp,handles.sssdpp];
            if ~isempty(SysC3),for i=1:length(ctrla), set(ctrla(i),'String',strtrim(char(SysC3(i,:))),'TooltipString',strtrim(char(SysC3(i,:))));end, end
            ctrla=[handles.dspksssdelfl, handles.dspksstdelfl, handles.dspkdecktfl, handles.dspklicortfl,handles.dspkequtfl, handles.dspkdeckpfl,...
                handles.dspklicorpfl,handles.dspkequpfl, handles.dspkdtfl,handles.dspkcondfl,handles.dspksssfl, handles.dspksstfl,...
                handles.dspkgflfl, handles.dspkwflfl];
            if ~isempty(dspkfla),for i=1:length(ctrla), set(ctrla(i),'value',7-dspkfla{i});despikefl(ctrla(i), 1, handles);end, end
        else %version>140
            if exist('hdllohistr') & exist('hdltype')
                if ~isempty(hdllohistr) & ~isempty(hdltype)
                    for i=1:10
                        handles.("atm" + num2str(i)).Value=0;handles.("equ" + num2str(i)).Value=0;
                        if ismember('ATM',hdltype{i}), handles.("atm" + num2str(i)).Value=true;end
                        if ismember('EQU',hdltype{i}), handles.("equ" + num2str(i)).Value=true;end
                        handles.("lo" + num2str(i)).String = hdllohistr{i,1};
                        handles.("hi" + num2str(i)).String = hdllohistr{i,2};
                    end
                end
            end
        end
        
        if exist('OTin','var'), set(handles.AOvalue,'string',num2str(OTin,'%10.2f'),'TooltipString',num2str(OTin,'%10.2f'));end
        
        %Checks Data Config file
        fnotok=1;
        while fnotok
            if exist([fpath filesep hfile],'file')==2%checks if config file exists in same folder as .mat file
                headp=fpath;
                mess=sprintf(['NOW USING DATA CONFIG FILE (''%s'')\n\n' ...
                    'LOCATED IN THE FOLDER BELOW!\n\n%s\n\n'],hfile,headp);
                mmm=msgbox(mess,'DATA Config File ');
                CenterWindow(handles,mmm,handles.figure1)
                uiwait(mmm);
            elseif exist([headp filesep hfile],'file')==2 %checks if config file exists where .mat file says it is
                %options.Interpreter='tex';
                mess=sprintf(['LOADED DATA CONFIG FILE (''%s'')\n\n' ...
                    'IS ALREADY ON THIS COMPUTER!\n\n%s\n\n'],hfile,[headp filesep hfile]);
                mmm=questdlg(mess,'DATA Config File EXISTS','Overwrite','Save Copy','Save Copy');
                if strcmp(mmm,'Save Copy'), Config_Data_Save(handles.SaveConfig, eventdata, handles);end
            else Config_Data_Save(handles.SaveConfig, eventdata, handles);
            end
            fnotok=(exist([headp filesep hfile],'file')~=2);
        end
        set(handles.loadhct,'String',hfile,'ToolTipString',hfile);drawnow;
        clear_savew=1;
        
        %Checks System Config ini file
        fnotok=1;
        while fnotok
            if exist([fpath filesep sysinif],'file')==2%checks if config file exists in same folder as .mat file
                sysinip=fpath;
                mess=sprintf(['NOW USING SYSTEM CONFIG FILE (''%s'')\n\n' ...
                    'LOCATED IN THE FOLDER BELOW!\n\n%s\n\n'],sysinif,sysinip);
                mmm=msgbox(mess,'SYSTEM Config File ');
                CenterWindow(handles,mmm,handles.figure1)
                uiwait(mmm);
            elseif exist([sysinip filesep sysinif],'file')==2
                mess=sprintf(['LOADED SYSTEM CONFIG FILE (''%s'')\n\n' ...
                    'IS ALREADY ON THIS COMPUTER!\n\n%s\n\n'],sysinif,[sysinip filesep sysinif]);
                mmm=questdlg(mess,'Data SYSTEM File EXISTS','Overwrite','Save Copy','Save Copy');
                if strcmp(mmm,'Save Copy'),Config_Sys_Save(-1, -1, handles);  end
            else Config_Sys_Save(-1, -1, handles);
            end
            fnotok=(exist([sysinip filesep sysinif],'file')~=2);drawnow;
        end
        set(handles.paraframe,'title',sysinif);
              
        if  vxaxiscok==1
            set(handles.ydayok,'visible','on');
            Status(handles.ydayok, eventdata, handles);
            if  posok==1
                set(handles.posok,'visible','on');  
                Status(handles.posok, eventdata, handles);
                if  stdok==1
                    set(handles.stdok,'visible','on');
                    Status(handles.stdok, eventdata, handles);
                    if  xco2corrok==1
                        set(handles.xco2aok,'visible','on');
                        if  atmok==1
                            Status(handles.xco2aok, eventdata, handles);
                            if  airiok==1
                                set(handles.cfco2,'Visible','on');end
                        else airiok=0; end
                        
                    else atmok=0;airiok=0; end
                else xco2corrok=0;atmok=0;airiok=0; end
            else stdok=0;xco2corrok=0;atmok=0;airiok=0; end
        else  posok=0;stdok=0;xco2corrok=0;atmok=0;airiok=0; end
        
        Stat_Update_Callback(0,0,handles);
        listf=DataA(1,1:end);
        listf=strtrim(listf);
        listf_calc=DataA(1,length(listf_file)+1:end-1);
        listf_calc=strtrim(listf_calc);
    else
        mmm=msgbox(sprintf('Invalid File Type.\nValid File Types are: *dat.txt ; *.csv ; *.xls ; *.mat.'),'NOTE','warn','modal');
        CenterWindow(handles,mmm,handles.figure1)
        uiwait(mmm);
        return;
    end % ************************************** end file type *************************************************************

    if vtype~= Nonecol 
        tblunk_modify(handles.tblunkmenu.Children(1), eventdata, handles) ;%adjust number of UNK in table from file info
    end
    
    npts=size(DataA,1);
    DataA=strtrim(DataA);
    
    if ~contains(fname,'.mat')%Add calculated Fields to ALL files except .mat files and files with uploaded fields
        stdtit='';for i=1:length(stdv),stdtit{1,i}=['Std' num2str(i) ' values'];end
        
        fieldtitles=[{'Speed(knots)'},{'Tin interp'},{'Teq-Tin'},{'Peq-Patm'},{'EQU Pressure (calc)'},{'Std Offset'},...
            stdtit,...
            {'xCO2 dry'},{'xCO2 corr'},{'xCO2 interp corr'},{'xCO2W Corrected'},{'xCO2A Corrected'},{...
            'fCO2 water'},{'fCO2 air'},{'fCO2 interp'},{'dfCO2(w-a)'}];
        qctitles=[{'QC'},{'SubFlag'},{'SubFlag_User'},{'YDay Calc'},{'Cruise Day'},{'none'}];
        
        %***** Data Config File Use *************
        %Assigns variable from titles to column number in file (listf_file) according to config file
        nvar=zeros(1,51);%Ini variables
        Config_Data_Load(-1, 0, handles, 0);
        for i=1:size(titles,1)%Ini variables
            nv=find(strcmpi(strtrim(listf_file),strtrim(titles(i,:)))==1);
            if length(nv)>0,   nvar(i)=nv(1);    end
        end
        
        ncol=size(DataA,2);
        numvar=length(vartitles);%these are the variables taken from the file + flag, subflag, subflaguser (See vartitles)
        %j=numvar+1;%+1 to leave a space in nvar for vsubfu
        j=numvar;
        DataA=strtrim(DataA);
        
        if sum(ismember([fieldtitles,qctitles],listf_file))>0 %calc or qc fields detected 
             mmm=questdlg(sprintf('Use the uploaded calculated or QC fields\n\nOR\n\nCreate new ones?'),'Calc Fields Detected','Use Uploaded','Create New','Use Uploaded');
             if strcmp(mmm,'Use Uploaded')
               
                 k=1;
                 for i=1:length(fieldtitles) %19
                     if ismember(fieldtitles(i),listf_file)
                         nvar(j+i)=find(strcmp(fieldtitles(i),listf_file));
                     else
                         nvar(j+i)=ncol+k;DataA(1,ncol+k)=fieldtitles(i);DataA(2:npts,ncol+k)={'-999'};k=k+1;
                     end                
                 end
                 i=i+1;
                 junk={'QC'};if ismember(junk,listf_file),nvar(numvar-2)=find(strcmp(junk,listf_file)); else DataA(1,ncol+k)=junk; DataA(2:npts,ncol+k)={'2'}; nvar(numvar-2)=ncol+k;k=k+1;end, i=i+1;%Sets Flags to 2
                 junk={'SubFlag'};if ismember(junk,listf_file),nvar(numvar-1)=find(strcmp(junk,listf_file)); else DataA(1,ncol+k)=junk; DataA(2:npts,ncol+k)={''}; nvar(numvar-1)=ncol+k;k=k+1;end, i=i+1;%Sets SubFlags to ''
                 junk={'SubFlag_User'};if ismember(junk,listf_file),nvar(numvar)=find(strcmp(junk,listf_file)); else DataA(1,ncol+k)=junk; DataA(2:npts,ncol+k)={''}; nvar(numvar)=ncol+k;k=k+1;end, i=i+1;%Sets SubFlags User to ''
                 junk={'YDay Calc'};if ismember(junk,listf_file),nvar(numvar+length(fieldtitles)+1)=find(strcmp(junk,listf_file)); else DataA(1,ncol+k)=junk; DataA(2:npts,ncol+k)={'-999'}; nvar(numvar+length(fieldtitles)+1)=ncol+k;k=k+1;end, i=i+1;%Sets YDay Calc column to -999
                 junk={'Cruise Day'};if ismember(junk,listf_file),nvar(numvar+length(fieldtitles)+2)=find(strcmp(junk,listf_file)); else DataA(1,ncol+k)=junk; DataA(2:npts,ncol+k)={'-999'}; nvar(numvar+length(fieldtitles)+2)=ncol+k;k=k+1;end, i=i+1;%Sets Cruise Day column to -999
                 junk={'none'};if ismember(junk,listf_file),Nonecol=find(strcmp(junk,listf_file)); else DataA(1,ncol+k)=junk; DataA(2:npts,ncol+k)={'-999'}; Nonecol=ncol+k;k=k+1;end, i=i+1;%Sets "None" column to -999
                 %when loading a "_Working.csv" file, empty values (like SubFlag or SubFlagUser) were replaced by '-999'. Need to be ''
                 DataA(:,find(strcmp('SubFlag',listf_file)))=strrep(DataA(:,find(strcmp('SubFlag',listf_file))),'-999','');
                 DataA(:,find(strcmp('SubFlag_User',listf_file)))=strrep(DataA(:,find(strcmp('SubFlag_User',listf_file))),'-999','');
             end
         end
        
        if (sum(ismember([fieldtitles,qctitles],listf_file))>0 & ~strcmp(mmm,'Use Uploaded')) | sum(ismember([fieldtitles,qctitles],listf_file))==0 % no previous calc fields in original file
            for i=1:length(fieldtitles) %19
                nvar(j+i)=ncol+i;
                DataA(1,ncol+i)=fieldtitles(i);
                DataA(2:npts,ncol+i)={'-999'};
            end
            i=i+1;
            DataA(1,ncol+i)={'QC'}; DataA(2:npts,ncol+i)={'2'}; nvar(numvar-2)=ncol+i; i=i+1;%Sets Flags to 2
            DataA(1,ncol+i)={'SubFlag'}; DataA(2:npts,ncol+i)={''}; nvar(numvar-1)=ncol+i; i=i+1;%Sets SubFlags to ''
            DataA(1,ncol+i)={'SubFlag_User'}; DataA(2:npts,ncol+i)={''}; nvar(numvar)=ncol+i; i=i+1;%Sets SubFlags User to ''
            DataA(1,ncol+i)={'YDay Calc'}; DataA(2:npts,ncol+i)={'-999'}; nvar(numvar+length(fieldtitles)+1)=ncol+i; i=i+1;%Sets YDay Calc column to -999
            DataA(1,ncol+i)={'Cruise Day'}; DataA(2:npts,ncol+i)={'-999'}; nvar(numvar+length(fieldtitles)+2)=ncol+i; i=i+1;%Sets Cruise Day column to -999
            DataA(1,ncol+i)={'none'}; Nonecol=ncol+i; DataA(2:npts,Nonecol)={'-999'}; %Sets "None" column to -999
        end
        nvar(nvar==0)=Nonecol; %replaces 0 by Nonecol
        
        listf=DataA(1,1:end);
        listf=strtrim(listf);
        listf_calc=DataA(1,ncol+1:end-1);
        listf_calc=strtrim(listf_calc);
        
        Update_List(handles);
        Update_Var(handles);
        
       %Check #STD and UNK in file and advise
       mmm2=Check_nSTD_from_file;
       if strcmp(mmm2,'1 - Program fix'),  tblstd_modify(handles.tblstdmenu.Children(1), eventdata, handles) ;end

        if (exist('mmm','var') & strcmp(mmm,'Use Uploaded')) & strcmp(DataA(1, vyday),'YDay Calc') & ~isempty(find(~strcmp(DataA(2:end, vyday),'-999'))) %uploaded YDay column has non -999 values
            dateok=1;
            set(handles.ydayok,'Visible','on');drawnow;nvar(1)=vyday;Update_Var(handles);
        end
    end
    
    ss='';samei_all=[];samei_dup=[];
    for i=1:size(DataA,2)
            samei=find(strcmpi(DataA(1,i),DataA(1,:))==1);
            if size(samei,2)>1 
                if max(samei)>ncol %one duplicate is a calculated column - should be last one (not renumbered)
                    for j=1:size(samei,2)-1,   DataA(1,samei(j)) = {[DataA{1,samei(j)} '_' strtrim(num2str(j+1))]};  end  %1st col renumbered starting at _2
                    ss='It could conflict with a column created by the program.';
                else
                    for j=2:size(samei,2),  DataA(1,samei(j)) = {[DataA{1,samei(j)} '_' strtrim(num2str(j))]};  end   %2nd col renumbered starting at _2
                end
                samei_dup=unique([samei_dup i]);
            end
    end
    if ~isempty(samei_dup)
        dupcol = cellfun(@(x) [x],DataA(1,samei_dup),'UniformOutput',false);
        dupcoltxt=sprintf('%s\n',dupcol{:});
        if ~isempty(ss)
            mmm=msgbox(sprintf('The following are duplicate header(s):\n\n%s\n%s\n\nThey have been modified\n\nto differentiate them.',dupcoltxt,ss),'Data Import...','Warn','modal');
        else
            mmm=msgbox(sprintf('The following are duplicate header(s):\n\n%s\nThey have been modified\n\nto differentiate them.',dupcoltxt),'Data Import...','Warn','modal');
        end
        CenterWindow(handles,mmm,handles.figure1)
        uiwait(mmm);
        listf=DataA(1,1:end);
        listf=strtrim(listf);
        listf_file=DataA(1,1:ncol);
        listf_calc=DataA(1,ncol+1:end-1);
        listf_calc=strtrim(listf_calc);
        
        Update_List(handles);
        Update_Var(handles);
    end
    
    Calc_Fields_Callback(handles);
    activity(1,handles);%Calc_Fields above turns it off...
    
    
    cla(handles.axes1);cla(handles.axes2);
    c=get(handles.Graphc_s,'BackgroundColor');
    set(handles.axes1,'XColor','r','YColor','r', 'XAxisLocation','bottom', 'YAxisLocation','left', 'Color','none');
    set(handles.axes2,'Position',get(handles.axes1,'Position'), 'XColor','k','YColor','k', 'XAxisLocation','top', 'YAxisLocation','right', 'Xtick',[],'Color',c);
    set(handles.axes1,'ButtonDownFcn',f1);set(handles.axes2,'ButtonDownFcn',f2);
    hp=[];hp2=[];
    set(handles.axes2,'XLim',[0 1],'YLim', [0 1]);
    
    on_list=[handles.merge_butt  handles.delete_col_butt  handles.btnCalcYDay  handles.btnSave  handles.btnXML...
        handles.SaveConfig  handles.lr handles.btnplot handles.OSet handles.OApply handles.OCheck];
    set(on_list,'Visible','on');drawnow;
    
    if lname==1
        set(handles.lr,'String',sprintf('file : %s', fname ),'ToolTipString',sprintf('file : %s', fname ));
    end
    
    
    if exist('delines', 'var')
        if ~isempty(delines) & delines>0
            mess=sprintf('number of skipped lines other than SHUTDOWN...etc:\n'); mess1='';
            for i=1:length(delines)
                if ischar(fnameo), mess2=fnameo; else mess2=fnameo{:,i}; end
                if delines(i)>0,   mess1=[mess1 sprintf('\n%d lines in %s',delines(i),mess2)]; end
            end
            mmm=msgbox([mess mess1],'File Info','warn','modal');
            CenterWindow(handles,mmm,handles.figure1)
            uiwait(mmm);
        end
    end
    
    % cd(ccd);
    Check_for_999(handles);
    Check_for_text(handles);
    activity(0,handles);
    
    clear mess mess1 mess2 iData
    
catch
    
    mmm=msgbox(lasterr,'Error');
    CenterWindow(handles,mmm,handles.figure1)
    uiwait(mmm);
    %restore conditions in case of fail
    if exist('DataO','var'), DataA=DataO; end
    stdok=stdok0;atmok=atmok0;vfco2ok=vfco2ok0; airiok=airiok0;vspeedok=vspeedok0;xco2corrok=xco2corrok0;motion=motion0;
    datap=datap0;headp=headp0;hfile=hfile0;sysinip=sysinip0;sysinif=sysinif0;
    
    on_list=[handles.merge_butt  handles.delete_col_butt  handles.btnCalcYDay  handles.btnSave handles.btnXML...
        handles.SaveConfig  handles.lr handles.btnplot handles.OSet handles.OApply handles.OCheck];
    set(on_list,'Visible','on');drawnow;
    activity(0,handles);

    clear  DataO stdok0 vfco2ok0 airiok0 vspeedok0 xco2corrok0 motion0 datap0;
    clear headp0 hfile0 sysinip0 sysinif0;

end



% --------------------------------------------------------------------
function Calc_Speed(handles)
%Calc speed from lat, long, yday...
global DataA vspeedok;
global  vxaxis vspeed vlat vlong Nonecol;
jstr={'Latitude';'Longitude';'x-axis Time'};
jlist=[vlat vlong vxaxis];
mmm=ismember([vlat vlong vxaxis],Nonecol);indx=find(jlist,Nonecol);

mess='Speed OK';cc='g';
% if  sum(ismember(indx0,find(tf)))>0 
%     for ii=indx0(ismember(indx0,find(tf))),   junk=sprintf('%s\n\"%s\"',junk,strtrim(ncstr(ii,:)));  end

if sum(mmm)>0 
    mess='';
    for ii=find(mmm),   mess=sprintf('%s\n\"%s\"',mess,strtrim(char(jstr{ii,:})));  end
    mess=sprintf('Speed not calculated.\n%s \n\nnot assigned yet!\n\nCorrect and Re-calc Fields',mess);
    vspeedok=0;cc='r';
%     mmm=msgbox(mess,'Speed Error','error','modal');
%     uiwait(mmm);
    set(handles.sp_f,'BackgroundColor',cc,'TooltipString',mess);
    return;
end

posdataok=find(strcmp(strtrim(DataA(2:end,vlat)),'-999')==0 & strcmp(strtrim(DataA(2:end,vlong)),'-999')==0 )+1;
ydayy0=str2num(char(DataA(posdataok(1:end-1),vxaxis)));
ydayy1=str2num(char(DataA(posdataok(2:end),vxaxis)));
latt0=str2num(char(DataA(posdataok(1:end-1),vlat)));
lonn0=str2num(char(DataA(posdataok(1:end-1),vlong)));
latt1=str2num(char(DataA(posdataok(2:end),vlat)));
lonn1=str2num(char(DataA(posdataok(2:end),vlong)));

ER=6371; KHtoKN=0.539957; %ER=Earth Radius in km  KHtoKN=km/h to knots factor
% ydayy0=str2num(char(DataA(2:end-1,vxaxis)));
% ydayy1=str2num(char(DataA(3:end,vxaxis)));
% latt0=str2num(char(DataA(2:end-1,vlat)));
% lonn0=str2num(char(DataA(2:end-1,vlong)));
% latt1=str2num(char(DataA(3:end,vlat)));
% lonn1=str2num(char(DataA(3:end,vlong)));

if isempty(latt0) | isempty(lonn0) | isempty(latt1) | isempty(lonn1) 
    mess=sprintf('Problem with Position Data.\nSpeed not calculated.');
    vspeedok=0;cc='r';
%     mmm=msgbox(mess,'Speed Error','error','modal');
%     uiwait(mmm);
    set(handles.sp_f,'BackgroundColor',cc,'TooltipString',mess);
    return;
end

speed=ER*acos(sind(latt0).*sind(latt1)+cosd(latt0).*cosd(latt1).*cosd(lonn1-lonn0));
speed=speed./((ydayy1-ydayy0)*24);
speed=speed*KHtoKN;
speed(find(imag(speed)))=0;

DataA(posdataok(2:end),vspeed)=cellstr(num2str(speed,'%1.1f'));
DataA(posdataok(1),vspeed)=DataA(posdataok(2),vspeed);
vspeedok=1;
set(handles.sp_f,'BackgroundColor',cc,'TooltipString',mess);



% --------------------------------------------------------------------
function minDT=Auto_SST_Offset(OTin)
global DataA  vxaxis vtin  vflag vteq;


        %nnn = indices of good sst and equ t data
        nnn=find(strcmp(DataA(:,vxaxis),'-999')==0 & strcmp(DataA(:,vtin),'-999')==0 & strcmp(DataA(:,vteq),'-999')==0 & strcmp(DataA(:,vflag),'2'));
        
        if length(nnn)>2
            xin=str2num(char(DataA(nnn,vxaxis)));
            yin=str2num(char(DataA(nnn,vtin)));
            xout=str2num(char(DataA(nnn,vxaxis)))-OTin/60/24;
            if (isempty(xin)|isempty(yin))
               minDT=-999;
            else
                yout=interp1(xin,yin,xout,'linear','extrap');
                minDT=sum(abs(yout-str2num(char(DataA(nnn,vteq)))));
                % DataA(nnn,vtini)=cellstr(num2str(yout));
                % vtiniok=1;mess='SST interpolated OK';cc='g';
            end
        else
            minDT=-999;
        end
        clear xin yin xout yout nnn irange;





% --------------------------------------------------------------------
function Calc_Fields_Callback(handles,varargin)
% Stub for Callback of the uicontrols handles.loadf.
global DataA  Nonecol vartitles;
global vdt vdp vtiniok vsubf vsubfu vflag sub_flag;
global  vpeq  vteq vpamb  vtin  vtini vpeqa  vpatm vxaxis xco2corrok;
% global datap  headp hfile  sysinip sysinif;
global stdv vstdo  vtype  vstdx  vlicorx vlicorxcorr vlicorw vxco2dry;

vtiniok=0;
% handles=varargin{1};
activity(1,handles);  drawnow;
% typea={'STD1';'STD2';'STD3';'STD4';...
%     'STD1-DRAIN';'STD2-DRAIN';'STD3-DRAIN';'STD4-DRAIN';...
%     'STD1z';'STD2z';'STD3z';'STD4z';...
%     'STD1s';'STD2s';'STD3s';'STD4s'};
typea={};
for i=1:length(stdv)
    typea=[typea;{['STD' num2str(i)]}];
end
xrange=[2:size(DataA,1)];
%varargin is the range over which to recalculate fields.
if ~isempty(varargin), xrange=[varargin{1}];end
%column needed to calc fields
nc=[vpeq vteq vpamb vtype vtin vlicorx vpatm vxaxis vlicorw];
ncstr=vartitles([6 7 11 12 14 16 18 1 17],:);
tf=ismember(nc,Nonecol);

%initiates Dry and Corrected xCO2 so something gets saved even when xco2 not dry or corrected.
%vlicorxcorr Will be replaced by real corrected one in 'corrxco2_calc'
if ~xco2corrok, DataA(xrange,vlicorxcorr)=DataA(xrange, vlicorx);end
%vxco2dry will be calculated below if necessary data available
DataA(xrange,vxco2dry)=DataA(xrange, vlicorx);

%Insitu T interpolated **********************************************************************************************************************
OTin=str2num(get(handles.AOvalue,'String'));
mess=' ';cc='r';
if  OTin==0,  DataA(xrange,vtini)=DataA(xrange,vtin);vtiniok=1; mess='SST interpolated OK';cc='g';end

indx0=[5 8];junk='';%vtin and vxaxis
if sum(ismember(indx0,find(tf)))>0 
    for ii=indx0(ismember(indx0,find(tf))),   junk=sprintf('%s\n\"%s\"',junk,strtrim(ncstr(ii,:)));  end
    mess=sprintf(['Time Offset cannot be applied to Insitu Temperature.'...
                     '\n%s \n\nnot assigned yet.'...
                     '\n\nCorrect and Re-Calc Fields'],junk);
%     mmm=msgbox(mess,'Calc Fields Error','error','modal');
%     uiwait(mmm);
    vtiniok=0;cc='r';
elseif OTin~=0 
     %increase range considered (xrange) by 10 for the time offset interpolation
        if min(xrange)>10, minx=min(xrange)-10;  else  minx=2;   end
        if max(xrange)<size(DataA,1)-10, maxx=max(xrange)+10;  else  maxx=size(DataA,1);   end
        irange=[minx:maxx];
        nnn=find(strcmp(DataA(irange,vxaxis),'-999')==0 & strcmp(DataA(irange,vtin),'-999')==0);
        if length(nnn)>2
            xin=str2num(char(DataA(irange(nnn),vxaxis)));
            yin=str2num(char(DataA(irange(nnn),vtin)));
            xout=str2num(char(DataA(irange(nnn),vxaxis)))-OTin/60/24;
            if (isempty(xin)|isempty(yin))
                msga={DataA(1,vxaxis);DataA(1,vtin)};
                mmm=msgbox(sprintf('no usable values found in %s.\nPossible cause is TEXT in data.\n Check your data.',char(msga{[isempty(xin) isempty(yin)],1})),'Calc Error','modal');
                CenterWindow(handles,mmm,handles.figure1)
                uiwait(mmm);
            elseif (size(xin,1)>size(yin,1))
                mmm=msgbox(sprintf('Some values in %s could not be converted to numbers.\nPossible cause is TEXT in data.\n Check your data.',char(DataA(1,vtin))),'Calc Error','modal');
                CenterWindow(handles,mmm,handles.figure1)
                uiwait(mmm);
            else
                yout=interp1(xin,yin,xout,'linear','extrap');
                DataA(irange(nnn),vtini)=cellstr(num2str(yout));
                %find xout on either side of each quest. SST values (subflag = 2)
                %gets rid of special characters
                ssearch=regexprep(char(sub_flag(2)),'(','\\('); ssearch=regexprep(ssearch,')','\\)'); ssearch=regexprep(ssearch,'+','\\+');
                qssti=regexp(DataA(irange(nnn),vsubf),ssearch);%find quest. sst subflag
                qssti=~cellfun('isempty',qssti);%0 and 1s
                qsstii=find(qssti);
                qissti=[];% quest. interpolated sst indices
                for i=1:length(qsstii)
                    [~,closesti]=min(abs(xin(qsstii(i)) + OTin/60/24 - xin));
                    qissti=[qissti closesti];
                    %1 quest sst generates 2 quest. interp. sst (before and
                    %after) unless its the first or last value
                    if (xin(qsstii(i)) + OTin/60/24 - xin(closesti))>0
                        if irange(closesti) < maxx, qissti=[qissti closesti+1];end
                    elseif (xin(qsstii(i)) + OTin/60/24 - xin(closesti))<0
                        if irange(closesti) > 2, qissti=[qissti closesti-1];end
                    end
                end
                qissti=sort(unique(qissti));
                if ~isempty(qissti)
                    irangen=irange(nnn);
                    Set_subflag(irangen(qissti),10,handles);Set_subflaguser(irangen(qissti),17,handles);
                end
                vtiniok=1;mess='SST interpolated OK';cc='g';
            end
        end
        clear irangen qsstii qssti qissti closesti ssearch xin yin xout yout nnn irange;
elseif OTin==0
    DataA(2:end,vtini)=DataA(2:end,vtin);
    vtiniok=1;mess='SST interpolated OK';cc='g';
end
set(handles.ssti_f,'BackgroundColor',cc,'TooltipString',mess);
datastr={'Type';'xCO2'; 'xH2O';'equ T';'SST';'equ P';'ambiant P';'atm P'};
types=(char(DataA(xrange, vtype)));
licorx=str2num(char(DataA(xrange, vlicorx)));
licorw=str2num(char(DataA(xrange, vlicorw)));
teq=str2num(char(DataA(xrange, vteq)));
tini=str2num(char(DataA(xrange, vtini)));
peq=str2num(char(DataA(xrange, vpeq)));
pamb=str2num(char(DataA(xrange, vpamb)));
patm=str2num(char(DataA(xrange, vpatm)));
badd=find([size(types,1) size(licorx,1) size(licorw,1) size(teq,1) size(tini,1) size(peq,1) size(pamb,1) size(patm,1)] ~= size(types,1));
if ~isempty(badd)
    %junk=['error in ' sprintf('\n%s',datastr{1:3,:})];
    mmm=msgbox(sprintf('Issue with calculating Fields.\n\nError in:%s\n\nCheck Data',sprintf('\n - %s',datastr{badd,:})),'Calc Fields Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1);
    uiwait(mmm);
    return;
end
%Adds DT, DP, STd Offset columns to DataA

%******************************************* Std Offset *********************************************
indx0=[4 6];junk='';%vtype and vlicorx

if  sum(ismember(indx0,find(tf)))>0 
    for ii=indx0(ismember(indx0,find(tf))),   junk=sprintf('%s\n\"%s\"',junk,strtrim(ncstr(ii,:)));  end
    mess=sprintf('Std Offsets cannot be calculated.\n%s \n\nnot assigned yet.\n\nCorrect and Re-Calc Fields',junk);
    %stdv=[];
else
    rdostdv_Callback(-1, -1, handles);
end

if ~isempty(stdv)
    typea={};
    for i=1:length(stdv)
        typea=[typea;{['STD' num2str(i)]}];
    end

    %Standard Offsets
    DataA(xrange,vstdo)={'-999'};
    for j=1:length(stdv)  %Select std
        if verLessThan('matlab','9.1'),ind=find(~cellfun(@isempty,strfind(DataA(xrange,vtype),typea{j})));
        else, ind=find(contains(DataA(xrange,vtype),typea{j}));end
            DataA(xrange(ind),vstdo)=cellstr(num2str(licorx(ind)-stdv(j),'%-1.3f'));
    end
    cc='g';mess='Std Offset OK';set(handles.stdo_f,'BackgroundColor',cc,'TooltipString',mess)
else
    cc='r';set(handles.stdo_f,'BackgroundColor',cc);
end


 
%******************************************* Delta T *********************************************
%DT
indx0=[2 5];junk='';%teq or tin
if  sum(ismember(indx0,find(tf)))>0 
    for ii=indx0(ismember(indx0,find(tf))),   junk=sprintf('%s\n\"%s\"',junk,strtrim(ncstr(ii,:)));  end
    mess=sprintf('Delta T cannot be calculated.\n%s \n\nnot assigned yet.\n\nCorrect and Re-Calc Fields',junk);
    cc='r';
%     mmm=msgbox(mess,'Calc Fields Error','error','modal');
%     uiwait(mmm);
    DataA(xrange,vdt)={'-999'};
%     return;
else
    DataA(xrange,vdt)={'-999'};
    ind=intersect(find(teq~=-999),find(tini~=-999)) ;
    DataA(xrange(ind),vdt)=cellstr(num2str(teq(ind)-tini(ind),'%-1.3f'));
    mess='Delta T OK';cc='g';
end
set(handles.dt_f,'BackgroundColor',cc,'TooltipString',mess);

%***************************** Adjust P Equ if differential by adding Ambient Pressure
DataA(xrange,vpeqa)={'-999'};
if  get(handles.rdoequpdiff,'value')==1
    ind=intersect(find(peq~=-999),find(pamb~=-999));
    DataA(xrange(ind),vpeqa)=cellstr(num2str(peq(ind)+pamb(ind),'%-1.3f'));
else
    DataA(xrange(ind),vpeqa)=DataA(xrange(ind),vpeq);
end
peqa= str2num(char(DataA(xrange(ind), vpeqa)));   

%******************************************* DP(Peq-Patm) *********************************************
%DP(Peq-Patm)
indx0=[1 3 7];junk='';%peq or pamb or patm
if  sum(ismember(indx0,find(tf)))>0  
    for ii=indx0(ismember(indx0,find(tf))),   junk=sprintf('%s\n\"%s\"',junk,strtrim(ncstr(ii,:)));  end
    mess=sprintf('Delta P cannot be calculated.\n%s \n\nnot assigned yet.\n\nCorrect and Re-Calc Fields',junk);
    cc='r';
%     mmm=msgbox(mess,'Calc Fields Error','error','modal');
%     uiwait(mmm);
    DataA(xrange,vdp)={'-999'};
%     return;
else
    DataA(xrange,vdp)={'-999'};
    ind=intersect(find(peqa~=-999),find(patm~=-999));
    DataA(xrange(ind),vdp)=cellstr(num2str(peqa(ind)-patm(ind),'%-1.3f'));
    mess='Delta P OK';cc='g';
end
set(handles.dp_f,'BackgroundColor',cc,'TooltipString',mess);

% *************************************Correct licor xCO2 for H2O signal if selected
mess='xCO2 dry N/A';cc=get(handles.figure1,'Color');
if  get(handles.rdodryxco2,'value')==1
    indx0=[6 9];junk='';%licorx and licorw
    if  sum(ismember(indx0,find(tf)))>0
        for ii=indx0(ismember(indx0,find(tf))),   junk=sprintf('%s\n\"%s\"',junk,strtrim(ncstr(ii,:)));  end
        mess=sprintf('xCO2 dry cannot be calculated.\n%s \n\nnot assigned yet.\n\nCorrect and Re-Calc Fields',junk);
        cc='r';
%         mmm=msgbox(mess,'Calc Fields Error','error','modal');
%         uiwait(mmm);
        %     return;
    else
        ind=intersect(find(licorx~=-999),find(licorw~=-999));
        DataA(xrange(ind),vxco2dry)=cellstr(num2str(licorx(ind)./(1-licorw(ind)/1000.),'%-1.3f'));
        mess='xCO2 dry OK';cc='g';
    end
end
set(handles.xco2d_f,'BackgroundColor',cc,'TooltipString',mess);


Calc_Speed(handles);
DataA(strcmp(DataA,'-999.000')==1)={'-999'};

activity(0,handles); 

clear types licorx licorw teq tini peq pamb patm peqa mess cc;
if exist('stdx','var'), clear stdx; end
%------------------------------------------------------------------------


function Update_Var(handles)
global listf vartitles nvar stdv;
global vxaxis vdutc vtutc vlat vlong vpeq vteq vwflo  ...
    vgflo  vlicorcav  vpamb  vtype  vsal  vtin  vstdx  vlicorx  ...
    vlicorw  vpatm vtcpu vtdbox vpdbox vtcond vflag vsubf vsubfu vspeed vtini vdt vdp vpeqa...
    vstdo vstdi ...
    vxco2dry vlicorxcorr vxco2icorr  vlicorxw vlicorxa ...
    vfco2w  vfco2a   vfco2i  vdfco2 vyday vydayi vcrday vcrdayi;

i=1; vstdi=[];
vxaxis=nvar(i);i=i+1;vdutc=nvar(i);i=i+1;vtutc=nvar(i);i=i+1;vlat=nvar(i);i=i+1;% 1 to 4
vlong=nvar(i);i=i+1;vpeq=nvar(i);i=i+1;vteq=nvar(i);i=i+1;vwflo=nvar(i);i=i+1;% 5 to 8
vgflo=nvar(i);i=i+1;vlicorcav=nvar(i);i=i+1;vpamb=nvar(i);i=i+1;vtype=nvar(i);i=i+1;% 9 to 12
vsal=nvar(i);i=i+1;vtin=nvar(i);i=i+1;vstdx=nvar(i);i=i+1;vlicorx=nvar(i);i=i+1;% 13 to 16
vlicorw=nvar(i);i=i+1;vpatm=nvar(i);i=i+1;vtcpu=nvar(i);i=i+1;% 17 to 19
vtdbox=nvar(i);i=i+1;vpdbox=nvar(i);i=i+1;vtcond=nvar(i);i=i+1;% 20 to 22
vflag=nvar(i);i=i+1;vsubf=nvar(i);i=i+1;vsubfu=nvar(i);i=i+1;% 23 to 25
%start of Calculated fields
vspeed=nvar(i);i=i+1;vtini=nvar(i);i=i+1;vdt=nvar(i);i=i+1;vdp=nvar(i);i=i+1;% 26 to 29
vpeqa=nvar(i);i=i+1;vstdo=nvar(i);i=i+1;
for j=1:length(stdv), vstdi(j)=nvar(i);i=i+1;end % 30 to 30+ num std
vxco2dry=nvar(i);i=i+1;vlicorxcorr=nvar(i);i=i+1;vxco2icorr=nvar(i);i=i+1;vlicorxw=nvar(i);i=i+1;vlicorxa=nvar(i);i=i+1; % 36 to 40?
vfco2w=nvar(i);i=i+1;vfco2a=nvar(i);i=i+1;vfco2i=nvar(i);i=i+1;vdfco2=nvar(i);i=i+1;% 41 to 44?
vyday=nvar(i);vydayi=i;i=i+1;vcrday=nvar(i);vcrdayi=i;% 45 to 46?

l=zeros(length(vartitles),'int8');%preallocating space
for i=1:length(vartitles)
    l(i)=length(strtrim(char(listf(nvar(i)))));
end
ll=max(l);fs=sprintf('%%10s - %%%ds',ll);
for i=1:length(vartitles)
    nlist(i,:)=sprintf(fs,vartitles(i,:),strtrim(char(listf(nvar(i)))));
end
set(handles.lb, 'String',nlist);



% --------------------------------------------------------------------
function Update_List(handles)
global listf listf_file listf_calc cmh cmhf cmhc cmhn Nonecol;

listf=strtrim(listf);
hv=[handles.CurrentX,handles.LeftY,handles.RightY];
for i=1:size(hv,2)
    olds=get(hv(i), 'String');
    oldv=olds(get(hv(i), 'Value'));nv=find(ismember(listf,oldv));
    set(hv(i), 'String',listf);
    if nv>0,set(hv(i), 'Value',nv);else set(hv(i),'Value',Nonecol);end
end
% set(handles.CurrentX, 'String',listf);
% set(handles.RightY, 'String',listf);
% set(handles.LeftY, 'String',listf);

if cmh~=0
    delete(cmh);cmh=[];
end
if ishandle(cmhf),delete(cmhf);cmhf=0;end
if ishandle(cmhc),delete(cmhc);cmhc=0;end
if ishandle(cmhn),delete(cmhn);cmhn=0;end

cmhf=uimenu('Parent',handles.C_list,'Label','File Columns...','Tag','cmhfile');
cmhc=uimenu('Parent',handles.C_list,'Label','From Calculated Fields...','Tag','cmhcalc');
cmhn=uimenu('Parent',handles.C_list,'Label','none','Callback', @assignvar,...
        'Tag',num2str(Nonecol));
for i=1:length(listf_file)
    cmh(i)=uimenu('Parent',cmhf,'Label',char(listf_file{i}),'Callback', @assignvar,...
        'Tag',num2str(i));
end
for j=1:length(listf_calc)
    cmh(j+length(listf_file))=uimenu('Parent',cmhc,'Label',char(listf_calc{j}),'Callback', @assignvar,...
        'Tag',num2str(j+length(listf_file)));
 end
% --------------------------------------------------------------------
function  plot_xco2_sst_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.xco2_sst.
global DataA ind1 ind2 Currentplot OTin;
global vxaxis vtype vtin vtini vtiniok vxco2dry vlicorxcorr xco2corrok vflag;


activity(1,handles); %drawnow;
equ_atm={'EQU';'ATM';'EQU-DRAIN';'ATM-DRAIN'};
equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};

if xco2corrok==1
    vlicorxplot=vlicorxcorr;
else
    vlicorxplot=vxco2dry;
end

vtinplot=vtini;
if OTin==0 | ~vtiniok,  vtinplot=vtin; end


set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1)
    uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vlicorxplot),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
        ind22{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vtinplot),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@plot_xco2_sst_Callback;
clear ind11 ind22;
plotdata_double(handles,ind1, ind2, vxaxis, vlicorxplot, vtinplot,h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end

activity(0,handles); 


% --------------------------------------------------------------------
function  plot_track_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.xco2_sst.
global DataA ind1 ind2 Currentplot  ;
global vspeed vlat vlong  vflag;


activity(1,handles); %drawnow;


%set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1)
    uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(strcmp(DataA(:,vflag),num2str(i+1))==1);
        ind22{i}=find(strcmp(DataA(:,vflag),num2str(i+1))==1);
        if double(get(handles.rdomv,'Value'))==0  %don't plot missing values
            indmv=find(strcmp(DataA(:,vlong),'-999')==0 & strcmp(DataA(:,vlat),'-999')==0);
            ind11{i}=intersect(ind11{i},indmv);
            indmv=find(strcmp(DataA(:,vlat),'-999')==0 & strcmp(DataA(:,vspeed),'-999')==0);
            ind22{i}=intersect(ind22{i},indmv);
        end
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
if isobject(h),Currentplot=@plot_track_Callback;end
if h==handles.plotposmap, Currentplot=@plot_trackpos_Callback; end
clear ind11 ind22;
plotdata_double(handles,ind1, ind2, vlong, vlat, vspeed, h);
    
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end

activity(0,handles); 

% --------------------------------------------------------------------

function  plot_trackpos_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.xco2_sst.

plot_track_Callback(h, eventdata, handles, varargin);


% --------------------------------------------------------------------
function  plot_T_DT_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.T_DT.
global DataA ind1 ind2 Currentplot  OTin;
global vxaxis  vteq   vtype  vtin vtini vflag  vdt vtiniok;


activity(1,handles); %drawnow;
equ={'EQU';'EQU-DRAIN'};

OTin=str2num(get(handles.AOvalue,'String'));
vtinplot=vtini;
if OTin==0 | vtiniok==0,  vtinplot=vtin; end

set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1)
    uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vtinplot),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
        ind22{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vtinplot),'-999')==0 & strcmp(DataA(:,vteq),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@plot_T_DT_Callback;
plotdata_double(handles,ind1, ind2, vxaxis, vtinplot, vdt,h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end

activity(0,handles); 


% --------------------------------------------------------------------
function  plot_tin_teq_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.T_DT.
global DataA ind1 ind2 Currentplot  OTin;
global vxaxis vteq vtin vtini vflag vtiniok vtype;


activity(1,handles); %drawnow;
equ={'EQU';'EQU-DRAIN'};

if ~isobject(h)
    tinplot=get(handles.plot_tin_teq, 'ForegroundColor');%when red, tinplot will be [1 0 0], not black ([0 0 0])
    if tinplot(1), h=handles.plot_tin_teq ; else h=handles.plot_tini_teq ;end
end

OTin=str2num(get(handles.AOvalue,'String'));
vtinplot=vtini;
if OTin==0 | ~vtiniok | h==handles.plot_tin_teq ,vtinplot=vtin;end
if h==handles.plot_tin_teq
    indt=1:size(DataA,1); %all indices
else
    indt=find(ismember(DataA(:,vtype),equ));

end


set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1)
    uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=intersect(indt,find(strcmp(DataA(:,vtinplot),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1))));
        ind22{i}=intersect(indt,find(strcmp(DataA(:,vteq),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1))));
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@plot_tin_teq_Callback;
plotdata_double(handles,ind1, ind2, vxaxis, vtinplot, vteq,h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end

activity(0,handles); 


% --------------------------------------------------------------------
function  plot_P_DP_Callback(h, eventdata, handles, varargin) %#ok<*INUSL>
% Stub for Callback of the uicontrol handles.P_DP.
global DataA ind1 ind2 Currentplot  vflag ;
global vxaxis  vpamb  vtype  vdp  vpeqa;


activity(1,handles); %drawnow;
equ={'EQU';'EQU-DRAIN'};

set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1);
    uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vpeqa),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
        ind22{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vpamb),'-999')==0 & strcmp(DataA(:,vpeqa),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@plot_P_DP_Callback;
plotdata_double(handles,ind1, ind2, vxaxis, vpeqa, vdp,h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end

activity(0,handles); 


function plot_Peq_Patm_Callback(h, eventdata, handles, varargin)
global DataA ind1 ind2 Currentplot ;
global vxaxis  vpeqa vpatm vflag vtype;


activity(1,handles); %drawnow;
equ={'EQU';'EQU-DRAIN'};

set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vpeqa),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
        ind22{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vpatm),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@plot_Peq_Patm_Callback;
plotdata_double(handles,ind1, ind2, vxaxis, vpeqa, vpatm,h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end

activity(0,handles); 


% --------------------------------------------------------------------
function  plot_xco2_gflow_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.xco2_gflow.
global DataA ind1 ind2 Currentplot ;
global vxaxis vgflo  vtype  vxco2dry vlicorxcorr xco2corrok vflag;


activity(1,handles); %drawnow;
equ={'EQU';'EQU-DRAIN'};

if xco2corrok==1
    vlicorxplot=vlicorxcorr;
else
    vlicorxplot=vxco2dry;
end

set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vlicorxplot),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
        ind22{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vgflo),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@plot_xco2_gflow_Callback;
clear ind11 ind22;
plotdata_double(handles,ind1, ind2, vxaxis, vlicorxplot, vgflo,h);

if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end
activity(0,handles); 

% --------------------------------------------------------------------
function  plot_xco2_h2oflow_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.xco2_gflow.
global DataA ind1 ind2 Currentplot ;
global vxaxis vwflo  vtype  vxco2dry vlicorxcorr xco2corrok vflag;


activity(1,handles); %drawnow;
equ={'EQU';'EQU-DRAIN'};

if xco2corrok==1
    vlicorxplot=vlicorxcorr;
else
    vlicorxplot=vxco2dry;
end

set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vlicorxplot),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
        ind22{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vwflo),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@plot_xco2_h2oflow_Callback;
clear ind11 ind22;
plotdata_double(handles,ind1, ind2, vxaxis, vlicorxplot, vwflo,h);

if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end
activity(0,handles); 



% --------------------------------------------------------------------
function plot_Stdoffset_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Stdoffset.
global DataA ind1  Currentplot ;
global vxaxis vtype  vflag vstdo;


activity(1,handles); %drawnow;

clear  ind11 ind12 ind2 ;
junk=get(handles.tblstdl,'Data');stdl=junk(:,1);
if (sum(cellfun(@sum,stdl))==0)
    mmm=msgbox('No Standards selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
for j=1:length(stdl) % for each std
    for i=1:3 %for each flag
        if stdl{j} & radio(i).Value
            if verLessThan('matlab','9.1'),ind11{j}{i}=find(~cellfun(@isempty,strfind(DataA(:,vtype),['STD' num2str(j)])) & strcmp(DataA(:,vstdo),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1))==1);
            else, ind11{j}{i}=find(contains(DataA(:,vtype),['STD' num2str(j)]) & strcmp(DataA(:,vstdo),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1))==1);end
        else
            ind11{j}{i}=[];
        end
    end
end
%zero and span only
for j=1:length(stdl) % for each std
    for i=1:3 %for each flag
        if stdl{j} & radio(i).Value
            ind12{j}{i}=find((~cellfun(@isempty,regexp(DataA(:,vtype),['^STD' num2str(j) '(?=[sz])']))) & strcmp(DataA(:,vstdo),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1))==1);
            ind12{j}{i}=find((~cellfun(@isempty,regexp(DataA(:,vtype),['^STD' num2str(j) '(?=[sz])']))) & strcmp(DataA(:,vstdo),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1))==1);
        else
            ind12{j}{i}=[];
        end     
    end
end

ind1=[]; %combine and sort all flags per std
for j=1:length(stdl) % for each std
    ind2=[];
    for i=1:3 %for each flag
       ind2=cat(1,ind2,ind11{j}{i});
    end
    ind2=sort(ind2);
    ind1{j,1}=ind2;
end
% add zero/span plots
ind2=[];
for j=1:length(stdl) % for each std
    for i=1:3 %for each flag
       ind2=cat(1,ind2,ind12{j}{i});
    end
end
ind2=sort(ind2);
ind1{j+1,1}=ind2;

legendstr=[];
for i=1:length(stdl) % for each std
    if length(ind1{i})>0
        if isempty(legendstr)
            legendstr=char(['STD' num2str(i) '     ']);
        else
            legendstr=cat(1,legendstr,char(['STD' num2str(i) '     ']));
        end
    end
end
legendstr=cat(1,legendstr,'zero/span');
legendstr={legendstr};
cols=['b';'r';'k';'c'];  lins=['- ';'--';': '];syms=['o';'^';'v'];
colsn=1;linsn=1;symsn=1;clear col;
for i=1:length(stdl) % for each std
    if colsn> length(cols)
        colsn=1;linsn=linsn+1;
        if linsn> length(lins)
            linsn=1;symsn=symsn+1;
            if symsn> length(syms), mmm=msgbox('max # of color-symbol-line reached (too many STDs','Error','error','modal');CenterWindow(handles,mmm,handles.figure1);uiwait(mmm);return;end
        end
    end
    col(i,:)=[lins(linsn,:) cols(colsn,:) syms(symsn,:)];
    colsn=colsn+1;
end
col(i+1,:)='  k*';col={col}; %for zero/span
Currentplot=@plot_Stdoffset_Callback;
clear indti indt indf ind11;
plotdata_single(handles,ind1, vxaxis, vstdo, col, legendstr, h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end
activity(0,handles); 


% --------------------------------------------------------------------
function  plot_SSS_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.SSS.
global DataA ind1  ind2 Currentplot Nonecol;
global vxaxis vtype  vsal  vflag;


activity(1,handles); %drawnow;
equ={'EQU';'EQU-DRAIN'};

clear ind11;

radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
         ind11{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vsal),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={[];[];[]};
clear indmv indf ind11;
Currentplot=@plot_SSS_Callback;

plotdata_double(handles,ind1, ind2, vxaxis, vsal, Nonecol,h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end
activity(0,handles); 


% --------------------------------------------------------------------
function plot_xco2air_atmflow_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.xco2air_atmflow.
global DataA ind1 ind2 Currentplot   ;
global vxaxis vgflo  vtype vxco2dry vlicorxcorr xco2corrok vflag ;


activity(1,handles); %drawnow;
atm={'ATM';'ATM-DRAIN'};

if xco2corrok==1
    vlicorxplot=vlicorxcorr;
else
    vlicorxplot=vxco2dry;
end

set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(ismember(DataA(:,vtype),atm) & strcmp(DataA(:,vlicorxplot),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
        ind22{i}=find(ismember(DataA(:,vtype),atm) & strcmp(DataA(:,vgflo),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@plot_xco2air_atmflow_Callback;
clear ind11 ind22;
plotdata_double(handles,ind1, ind2, vxaxis, vlicorxplot, vgflo,h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end

activity(0,handles); 


function plot_xco2air_Callback(h, eventdata, handles, varargin)
global DataA ind1 ind2 Currentplot Nonecol;
global vxaxis vtype xco2corrok vxco2dry vlicorxcorr vflag ;


activity(1,handles); %drawnow;
atm={'ATM';'ATM-DRAIN'};

if xco2corrok==1
    vlicorxplot=vlicorxcorr;
else
    vlicorxplot=vxco2dry;
end

radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(ismember(DataA(:,vtype),atm) & strcmp(DataA(:,vlicorxplot),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={[];[];[]};
clear indmv indf ind11;
Currentplot=@plot_xco2air_Callback;

plotdata_double(handles,ind1, ind2, vxaxis, vlicorxplot, Nonecol,h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end
activity(0,handles); 



function plot_fco2sw_fco2atm_Callback(h, eventdata, handles, varargin)
global DataA ind1 ind2 Currentplot   ;
global vxaxis  vfco2w vfco2i vflag vfco2ok vtype;


activity(1,handles); %drawnow;
equ={'EQU';'EQU-DRAIN'};

if vfco2ok==0
    mmm=msgbox(sprintf('fCO2 not calculated yet.\nUse the ''Calc fCO2'' button'),'Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vfco2w),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
        ind22{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vfco2i),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@plot_fco2sw_fco2atm_Callback;
clear ind11 ind22;
plotdata_double(handles,ind1, ind2, vxaxis, vfco2w, vfco2i, h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end

activity(0,handles); 



function plot_dfco2_tin_Callback(h, eventdata, handles, varargin)
global DataA ind1 ind2 Currentplot OTin;
global vxaxis vdfco2 vtin vtini vtiniok vflag vfco2ok vtype;


activity(1,handles); %drawnow;
equ={'EQU';'EQU-DRAIN'};

if vfco2ok==0
    mmm=msgbox(sprintf('fCO2 not calculated yet.\nUse the ''Calc fCO2'' button'),'Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

vtinplot=vtini;
if OTin==0 | vtiniok==0,  vtinplot=vtin; end


set(handles.rdomv,'Value',0); %Missing values not plotted
radio=[handles.rdoFlag2 handles.rdoFlag3 handles.rdoFlag4 ];
if (sum([radio.Value])==0)
    mmm=msgbox('No Flags selected','Plot Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
for i=1:3
    if radio(i).Value
        ind11{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vdfco2),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
        ind22{i}=find(ismember(DataA(:,vtype),equ) & strcmp(DataA(:,vtinplot),'-999')==0 & strcmp(DataA(:,vflag),num2str(i+1)));
    else
        ind11{i}=[];
        ind22{i}=[];
    end
end

ind1={ind11{1};ind11{2};ind11{3}};
ind2={ind22{1};ind22{2};ind22{3}};
Currentplot=@plot_dfco2_tin_Callback;
clear ind11 ind22;
plotdata_double(handles,ind1, ind2, vxaxis, vdfco2, vtinplot,h);
if isobject(h)
    set(findobj(handles.figure1,'-regexp','Tag','plot'),'ForegroundColor',[0 0 0]);
    set(h, 'ForegroundColor',[1 0 0]);
end

activity(0,handles); 



% --------------------------------------------------------------------
function axes2_ButtondownFcn(h, eventdata, handles, varargin)
global x y ;
% Stub for ButtondownFcn of the axes handles.axes2.
[x,y]=gpos(handles.figure1,h);
% --------------------------------------------------------------------
function  axes1_ButtondownFcn(h, eventdata, handles, varargin)
global x y ;
% Stub for ButtondownFcn of the axes handles.axes1.

[x,y]=gpos(handles.figure1,handles.axes1);
% --------------------------------------------------------------------
function  figure1_WindowButtonDownFcn(h, eventdata, handles, varargin)
% Stub for WindowButtonDownFcn of the figure handles.figure1.
global x y x0 y0 souris;

% s=get(handles.popGselect,'Value');
% if s==1 & length(char(get(handles.popGselect,'String')))>1
%     [x,y]=gpos(handles.figure1,handles.axes1);
% else
%     [x,y]=gpos(handles.figure1,handles.axes2);
% end
ax=get(handles.figure1,'Currentaxes');
[x, y]=gpos(handles.figure1,ax);
xlimits = get(ax,'XLim');
ylimits = get(ax,'YLim');
xv=[xlimits(1);xlimits(1);xlimits(2);xlimits(2)];
yv=[ylimits(1);ylimits(2);ylimits(2);ylimits(1)];
if inpolygon(x,y,xv,yv)
    souris='down';
    x0=x;y0=y;
end

% --------------------------------------------------------------------
function figure1_WindowButtonUpFcn(h, eventdata, handles, varargin)
% Stub for WindowButtonUpFcn of the figure handles.figure1.
global souris rh rfh nrec drag;

   if drag % rectangle was drawn
       ax=get(handles.figure1,'Currentaxes');
       [x,y]=gpos(handles.figure1,ax);
       xlimits = get(ax,'XLim');
       ylimits = get(ax,'YLim');
       if x<=max(xlimits) & x>=min(xlimits) & y<=max(ylimits) & y>=min(ylimits)
           if isempty(nrec),nrec=0;end
           nrec=nrec+1;
           rfh=[rfh rh];
       end
       drag=0;
   end

   souris='up';

% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DataA hp hp2 ind1 ind2 motion  x0 y0 x2 y2 rh rfh souris vsubfu drag ;

resolution=200;

if motion==0 | (isempty(hp) & isempty(hp2)), return; end
    %     set(handles.figure1,'WindowButtonMotionFcn','');
    %     reset(handles.figure1);

    pti=[];
    ax=get(handles.figure1,'Currentaxes');
    [x,y]=gpos(handles.figure1,ax);   
    xlimits = get(ax,'XLim');
    ylimits = get(ax,'YLim');
    if x<=max(xlimits) & x>=min(xlimits) & y<=max(ylimits) & y>=min(ylimits)
        set(handles.currx,'Visible','on');
        set(handles.curry,'Visible','on');
        set(handles.currx,'String',sprintf('x = %.4f', x  ));
        set(handles.curry,'String',sprintf('y = %.3f', y  ));
        xp = (xlimits(2)-xlimits(1))/resolution;
        yp = (ylimits(2)-ylimits(1))/resolution;
        if ax==handles.axes1
            ph=hp;ind=ind1;
            yax=handles.LeftY;
        end
        if ax==handles.axes2
            ph=hp2;ind=ind2;
            yax=handles.RightY;
        end
        %Draws Rectangle
        if strcmp(char(souris),'down')==1
            drag=1;
            [x2,y2]=gpos(handles.figure1,ax);
            xlim=get(ax,'XLim');
            ylim=get(ax,'YLim');
            cond1=(x0<=xlim(2) & x0 >= xlim(1) & y0 <= ylim(2) & y0 >= ylim(1));
            cond2=(x2<=xlim(2) & x2 >= xlim(1) & y2 <= ylim(2) & y2 >= ylim(1));
            cond3=(x0~=x2 & y0~=y2);
            if cond1 & cond2 & cond3
                rh=findobj(ax,'Type','Rectangle');
                if isobject(rh)
                    modifier=get(handles.figure1,'CurrentModifier');
                    if size(modifier,2)>0 & sum(ismember(modifier,{'shift','alt','command'}))>0
                        delete(rh(~ismember(rh,rfh)));
                    else
                        delete(rh); nrec=0;clearvars -global rfh;
                    end
                end
                rh=rectangle('Position',[min(x0,x2),min(y0,y2),abs(x2-x0),abs(y2-y0)],'Parent',ax);
                set(ax,'XLim',xlim,'YLim',ylim);
            end
        end
        if motion==2
            %Searches for point in data
            for i =1:length(ph)
                pti=find(abs(str2num(char(DataA(ind{i},handles.CurrentX.Value)))-x)<xp & abs(str2num(char(DataA(ind{i},get(yax,'Value'))))-y)<yp);
                if ~isempty(pti), break,end
            end
            xx=[];yy=[];
            if ~isempty(pti)
                xx=str2num(char(DataA(ind{i}(pti),handles.CurrentX.Value)));
                yy=str2num(char(DataA(ind{i}(pti),get(yax,'Value'))));
            end
            if ~isempty(xx) & ~isempty(yy)
                set (handles.figure1,'Pointer','crosshair');
                set(handles.currxy,'Visible','on');
                if length(xx)==1 & length(yy)==1
                    subfstr=char(DataA(ind{i}(pti),vsubfu));
                    if ~length(subfstr)
                        subfstr='None';
                    end
                    set(handles.currxy,'String',sprintf(' x = %-.4f\n y = %-.3f \nSub Flag: %s',[ xx  yy],subfstr));
                else
                    subfstr='N/A';
                    set(handles.currxy,'String',sprintf(' %i points\n selected',length(xx)));
                end
                bs=get(handles.currxy,'Position');
%                 set(handles.currxy,'Position',cat(2,get(handles.figure1,'CurrentPoint')+[0.5*bs(4) -bs(4) ],[bs(3) bs(4)]));
                set(handles.currxy,'Position',cat(2,get(handles.figure1,'CurrentPoint')+[0.5*bs(4) -bs(4) ],[0.06+length(subfstr)/200 bs(4) ]));
            else
                set (handles.figure1,'Pointer','arrow');
                set(handles.currxy,'Visible','off');
            end
        end
    else
        set (handles.figure1,'Pointer','arrow');
        set(handles.currxy,'Visible','off');
        set(handles.currx,'Visible','off');
        set(handles.curry,'Visible','off');
    end



% --------------------------------------------------------------------
function lstGselect_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.lstGselect.
%graph_select
 global VarO omo InferTO; %from data inference

s=get(h,'Value');    c=get(handles.Graphc_s,'BackgroundColor');

if s==1 & size(char(get(h,'String')),1)==1 %for single plots
    rh=findobj(handles.axes2,'Type','Rectangle');
    if ~isempty(rh);delete(rh);end
    set(handles.axes2,'Visible','off');
    axes(handles.axes1);
    uistack(handles.axes1,'top');
    set(handles.axes1,'XColor','r','YColor','r', 'XAxisLocation','bottom', 'YAxisLocation','left', 'Color',c);
    return;
end

if s==1 & size(char(get(h,'String')),1)>1
    rh=findobj(handles.axes2,'Type','Rectangle');
    if ~isempty(rh);delete(rh);end
    axes(handles.axes1);
    uistack(handles.axes1,'top');
    set(handles.axes1,'XColor','r','YColor','r', 'XAxisLocation','bottom', 'YAxisLocation','left', 'Color','none');
    set(handles.axes2,'Position',get(handles.axes1,'Position'),'XColor','K','YColor','K','XAxisLocation','top','YAxisLocation','right', 'Xtick',[],...
        'NextPlot','replace','Color',c);

else
    rh=findobj(handles.axes1,'Type','Rectangle');
    if ~isempty(rh);delete(rh);end
    axes(handles.axes2);
    uistack(handles.axes2,'top');
    set(handles.axes1,'XColor','r','YColor','r', 'XAxisLocation','bottom', 'YAxisLocation','left', 'Color',c);
    set(handles.axes2,'Position',get(handles.axes1,'Position'), 'XColor','K','YColor','K',...
        'XAxisLocation','top', 'YAxisLocation','right', 'Xtick',[], 'Color','none');

end

%updates data inference menu if necessary
if isempty(VarO), omo=1;
elseif (VarO(1)==handles.LeftY.Value),omo=1;elseif (VarO(2)==handles.LeftY.Value),omo=-1;
end
prec=int2str(2*(abs(InferTO)>0.5)+3*(abs(InferTO)<=0.5));
frmt=['%1.' prec 'f'];
mess=sprintf(['Apply Offset(~ ' frmt ') on Selected Data'],omo*InferTO);
set(handles.gd1ao,'Label',mess);
mess=sprintf(['Apply Offset(~ ' frmt ') on Selected Data'],-omo*InferTO);
set(handles.gd2ao,'Label',mess);


% --------------------------------------------------------------------
function plotdata_double(handles,index1, index2, vx, vy1, vy2,h)
% Stub for Callback of the uicontrol handles.CurrentX.
global DataA f1 f2 hp hp2  motion Nonecol xco2corrok airiok vfco2ok Currentplot;
global vstdi  vlicorxcorr   vfco2w  vfco2a vxco2icorr vfco2i vdfco2;
global map_xlim map_ylim  x1lim y1lim;
global vlong vlat dotsize hcb rfh nrec;
global newlong1 newlong2 mapfn map_loaded map_folder shorelines;
% Stub for Callback of the uicontrol handles.xco2_sst.
offon={'off' 'on'};
p3d=get(handles.p3d,'Value');rotate3d(handles.figure1, 'off');
t=0;t1=0;t2=0;tt='';
for i=1:3, t=t+length(index1{i})+length(index2{i});t1=t1+length(index1{i});t2=t2+length(index2{i});end
if ~t, tt='both variables.';
elseif ~t1 && vy1 ~= Nonecol, tt=['' char(DataA(1,vy1)) ''];
elseif ~t2 && vy2 ~= Nonecol, tt=['' char(DataA(1,vy2)) ''];
end
if ~isempty(tt), mmm=msgbox(['No data to plot for ' tt],'Data Error');CenterWindow(handles,mmm,handles.figure1) ;uiwait(mmm);end
if ~t, return;end

coul=[get(handles.Flagc12,'BackgroundColor'); get(handles.Flagc13,'BackgroundColor');get(handles.Flagc14,'BackgroundColor')];
coul2=[get(handles.Flagc22,'BackgroundColor'); get(handles.Flagc23,'BackgroundColor');get(handles.Flagc24,'BackgroundColor')];

%button to flag quest. air values visible or not 
junka=[handles.tQuestAirVal_Button handles.QuestAirVal_Button];
junkba=[handles.plotxco2air_atmflow handles.plot_xco2air];
if isobject(h), set(junka,'Visible',offon{1+ismember(h,junkba)});end

%button to zoom in or out when map is plotted and check box to center map at date line
junka=[handles.map_zoom_out];
junkba=[handles.plotposmap];
if isobject(h), set(junka,'Visible',offon{1+(ismember(h,junkba)|| p3d)});end

%check box to center map at date line
junka=[handles.rdodateline];
if isobject(h), set(junka,'Visible',offon{1+(vx==vlong || p3d)});end

if  vx== Nonecol & h~=handles.btnplot
    mmm=msgbox(sprintf('\"X Data\" column set to \"None\".\nData needs to be assigned \nbefore it can be plotted.'),'Plotting Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
if   vy1== Nonecol & h~=handles.btnplot
    mmm=msgbox(sprintf('\"Left Y Data\" column set to \"None\".\nData needs to be assigned \nbefore it can be plotted.'),'Plotting Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
%deals with missing x values (y values taken care of in individual callbacks
indmvx=find(strcmp(DataA(2:end, vx),'-999')==0)+1;
if   isempty(indmvx)
    mmm=msgbox(sprintf('\"X Data\" (%s) are ALL equal to \"-999\".\nCorrect and Re-plot.',char(DataA(1, vx))),'Plotting Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end% indmvy1=find(strcmp(DataA(:, vy1),'-999')==0);
% indmvy2=[];
% if vy2>0,indmvy2=find(strcmp(DataA(:, vy2),'-999')==0);end

 if double(get(handles.rdomv,'Value'))==0  %don't plot missing values
    for i=1:3
        if ~isempty(index1), index1{i}=intersect(index1{i},indmvx);end%index1{i}=intersect(index1{i},indmvy1);
        if ~isempty(index2), index2{i}=intersect(index2{i},indmvx);end%index2{i}=intersect(index2{i},indmvy2);
    end
 end
%if data part of xco2 correction, check that calc has been done
xco2c=[vstdi  vlicorxcorr];
%tf = 0;indx = 0;
[tf,indx]=ismember(vy1,xco2c);
if  tf==1 & xco2corrok==0
    mmm=msgbox(sprintf('\"%s\" cannot be plotted...Data not calculated yet.\nUse the ''Correct xCO2'' button',char(DataA(1, xco2c(indx)))),'Plotting Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
%if data part of air interpolation, check that calc has been done
xco2ic=vxco2icorr;
%tf = 0;indx = 0;
[tf,indx]=ismember(vy1,xco2ic);
if  tf==1 & airiok==0
    mmm=msgbox(sprintf('\"%s\" cannot be plotted...Data not calculated yet.\nUse the ''Interpolate Air xCO2'' button',char(DataA(1, xco2ic(indx)))),'Plotting Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

%if data part of fco2 calc, check that calc has been done
fco2c=[vfco2w  vfco2a vfco2i vdfco2];
%tf = 0;indx = 0;
[tf,indx]=ismember(vy1,fco2c);
if  tf==1 & vfco2ok==0
    mmm=msgbox(sprintf('\"%s\" cannot be plotted...Data not calculated yet.\nUse the ''Calc fCO2'' button',char(DataA(1, fco2c(indx)))),'Plotting Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

axes(handles.axes1)
set(handles.axes1,  'Visible','on', 'NextPlot','replace');
if ~isobject(h)
    x1lim=get(handles.axes1,'XLim');
    y1lim=get(handles.axes1,'YLim');
end
if length(index1{1})+ length(index1{2})+ length(index1{3})==0
    if ~isempty(hp),cla(handles.axes1,'reset');hp=[];clearvars -global rfh;nrec=0;end
end

hp=[];newlong1=cell(1,3);
for i=1:3
    if length(index1{i})>0 & ~p3d
        if (vx==vlong & get(handles.rdodateline,'value'))
            newlong1{:,i}=str2num(char(DataA(index1{i}, vx)));
            newlong1{1,i}(newlong1{:,i}<-3)=newlong1{1,i}(newlong1{:,i}<-3)+360;
            hp(i)=plot(newlong1{:,i},str2num(char(DataA(index1{i}, vy1))),'Parent',handles.axes1,...
            'Marker','.', 'MarkerSize',dotsize,'MarkerEdgeColor',coul(i,:), 'MarkerFaceColor',coul(i,:),'LineStyle','none');
        else
            hp(i)=plot(str2num(char(DataA(index1{i}, vx))),str2num(char(DataA(index1{i}, vy1))),'Parent',handles.axes1,...
                'Marker','.', 'MarkerSize',dotsize,'MarkerEdgeColor',coul(i,:), 'MarkerFaceColor',coul(i,:),'LineStyle','none');
        end
        %moves plot from axes
        yy=str2num(char(DataA(index1{i}, vy1)));
        if ~exist('ymin','var'), ymin=min(yy);ymax=max(yy);else, ymin=min(ymin,min(yy));ymax=max(ymax,max(yy));end
        xx=str2num(char(DataA(index1{i}, vx)));
        if ~exist('xmin','var'), xmin=min(xx);xmax=max(xx);else, xmin=min(xmin,min(xx));xmax=max(xmax,max(xx));end
        set(hp(i),'HitTest','off');
        if ~isobject(h)
            set(handles.axes1,'XLim',x1lim,'YLim',y1lim);
        end
        set(handles.axes1,  'NextPlot','add');
    end
end
if ~p3d
    y1lim=handles.axes1.YLim;
    if abs(ymin-y1lim(1))< 0.1*(ymax-ymin),   set(handles.axes1,'YLim',[ymin-0.1*(ymax-ymin) y1lim(2)]);    end
    if abs(ymax-y1lim(2))< 0.1*(ymax-ymin),   set(handles.axes1,'YLim',[y1lim(1) ymax+0.1*(ymax-ymin)]);    end
    x1lim=handles.axes1.XLim;
    if abs(xmin-x1lim(1))< 0.1*(xmax-xmin),   set(handles.axes1,'XLim',[xmin-0.1*(xmax-xmin) x1lim(2)]);    end
    if abs(xmax-x1lim(2))< 0.1*(xmax-xmin),   set(handles.axes1,'XLim',[x1lim(1) xmax+0.1*(xmax-xmin)]);    end
end
                                                                                        set(handles.axes1,'HitTest','on');
if ~p3d
    set(get(handles.axes1,'XLabel'),'String',char(DataA(1, vx)))
    set(get(handles.axes1,'YLabel'),'String',char(DataA(1, vy1)))
end
set(handles.LeftY, 'Value', vy1);
set(handles.CurrentX, 'Value', vx);

if length(index1{1})+length(index1{2})+length(index1{3}) > 0
    
    set(handles.axes1,  'Visible','on');
    set(handles.axes1,'ButtonDownFcn',f1);
    set(handles.axes1,'XColor','r','YColor','r', 'XAxisLocation','bottom', 'YAxisLocation','left');
    
    if ~p3d & (h==handles.plotposmap || (~isobject(h) && strcmp(func2str(Currentplot),'plot_trackpos_Callback')))
        
        x1lim=get(handles.axes1,'XLim');
        y1lim=get(handles.axes1,'YLim');
  
        if ~isempty(dir([map_folder filesep 'gshhs_*.b']))
            if isempty(shorelines) || (handles.map_res.Value~=map_loaded)
                map_loaded=handles.map_res.Value;   maps=['l';'h';'f'];
                if isempty(dir([map_folder filesep 'gshhs_' maps(map_loaded,:) '.b']))
                    mmm=msgbox(sprintf('Corresponding Map file (%s) NOT FOUND!\n Select another resolution!',['gshhs_' maps(map_loaded,:) '.b']));
                    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                    return;
                end
                handles.plotposmap.String='Loading map...';drawnow;
                shorelines = gshhs(mapfn(map_loaded,:));
                handles.plotposmap.String='Track w/ map';drawnow;
            end
            map_loaded = 2-(length(shorelines)<50000)+(length(shorelines)>180000);
 
            long=[shorelines.Lon];lat=[shorelines.Lat];
            if get(handles.rdodateline,'value')==1
                long(long<-3)=long(long<-3)+360;
                long(abs(long(1,1:end-1)-long(1,2:end))>300)=NaN;%breaks polygons across dateline
            end
            
        elseif exist('map','dir')==7
            load coast;%loads variables lat and long of coast
        else

            return;
        end
        
        %TO DO: plot using ortho projection when lat >80 or 85
        % axesm ('ortho', 'Frame', 'on', 'Grid', 'on');
        % geoshow(shorelines,'FaceColor',[1 1 .5],'EdgeColor',[.6 .6 .6]);
        % tissot;
        hh=plot(long,lat,'b','Parent',handles.axes1, 'Marker','none', 'MarkerSize',dotsize, 'MarkerEdgeColor','b', 'MarkerFaceColor','none','LineStyle','-');

        map_xlim = get(handles.axes1,'XLim');
        map_ylim = get(handles.axes1,'YLim');
        
        set(handles.axes1,'XLim',x1lim,'YLim',y1lim);
        set(hh,'HitTest','off');
        set(handles.map_zoom_out, 'value',0);
        index2=[];
        
    end    

    if p3d

        index2=[];

        hwb = waitbar(0,'\fontname{technical}Preparing Data...');
        CenterWindow(handles,hwb,handles.figure1);
        xx=[];yy=[];zz=[];
        indmv=find(strcmp(DataA(:, vlong),'-999')==0 & strcmp(DataA(:, vlat),'-999')==0);
        for i=1:3
            if ~isempty(index1{i})
                index1{i}=intersect(index1{i},indmv);
                xx=cat(1,xx,str2num(char(DataA(index1{i}, vlong))));
                yy=cat(1,yy,str2num(char(DataA(index1{i}, vlat))));
                zz=cat(1,zz,str2num(char(DataA(index1{i}, vy1))));
            end
        end
        
        if (get(handles.rdodateline,'value')), xx(xx<-3)=xx(xx<-3)+360;end

        waitbar(1/4,hwb,'Plotting scatter data...');
        mx=max(zz);    zzc=zz/mx;
        scatter3(handles.axes1,xx,yy,zz,4,zzc);
        hold (handles.axes1,'on');

        waitbar(2/4,hwb,'Color Bar...');
        hcb=colorbar('peer',handles.axes1); set(hcb,'location','eastoutside','yaxislocation','right');
        yt=get(hcb,'YTick');ytl='';
        for i=1:length(yt)
            ytl=strvcat(ytl,int2str(yt(i)*mx));
        end
        set(hcb,'YTickLabel',ytl);
        
        x1lim=get(handles.axes1,'XLim');
        y1lim=get(handles.axes1,'YLim');
        waitbar(3/4,hwb,'Plotting Coast...');

        if ~isempty(dir([map_folder filesep 'gshhs_*.b']))
            if isempty(shorelines) || (handles.map_res.Value~=map_loaded)
                map_loaded=handles.map_res.Value;   maps=['l';'h';'f'];
                if isempty(dir([map_folder filesep 'gshhs_' maps(map_loaded,:) '.b']))
                    mmm=msgbox(sprintf('Corresponding Map file (%s) NOT FOUND!\n Select another resolution!',['gshhs_' maps(map_loaded,:) '.b']));
                    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                    return;
                end
                handles.plotposmap.String='Loading map...';drawnow;
                shorelines = gshhs(mapfn(map_loaded,:));
                handles.plotposmap.String='Track w/ map';drawnow;
            end
            map_loaded = 2-(length(shorelines)<50000)+(length(shorelines)>180000);
 
            long=[shorelines.Lon];lat=[shorelines.Lat];
            if get(handles.rdodateline,'value')==1
                long(long<-3)=long(long<-3)+360;
                long(abs(long(1,1:end-1)-long(1,2:end))>300)=NaN;%breaks polygons across dateline
            end
            
        elseif exist('map','dir')==7
            load coast;%loads variables lat and long of coast
        else

            return;
        end
        
        map_height=mean(str2num(char(DataA(index1{1}, vy1)))); % average of plotted data flag 2 only
        hh=plot3(long,lat,zeros(length(long),1)+map_height,'Parent',handles.axes1, 'Marker','none', 'MarkerSize',dotsize, 'MarkerEdgeColor','b', 'MarkerFaceColor','none','LineStyle','-');

        map_xlim = get(handles.axes1,'XLim');
        map_ylim = get(handles.axes1,'YLim');
        
        set(handles.axes1,'XLim',x1lim,'YLim',y1lim);
        set(get(handles.axes1,'XLabel'),'String',char(DataA(1, vlong)))
        set(get(handles.axes1,'YLabel'),'String',char(DataA(1, vlat)))
        set(get(handles.axes1,'ZLabel'),'String',char(DataA(1, vy1)))

        set(hh,'HitTest','off');
        set(handles.map_zoom_out, 'value',0);
        index2=[];
        waitbar(4/4,hwb,'Done...');
        close (hwb);
        rotate3d(handles.figure1, 'on');

    end%end if p3d
    
    if motion==0
        motion=1;
    end
    
else % no data in index1, most likely set to none
    set(handles.axes1,  'Visible','off');
    listg=get(handles.popGselect, 'String');
    listg(1)=[(DataA(1, vy1))]; %#ok<NBRAK>
    set(handles.popGselect, 'String',listg,'Value',2);
    %     set(handles.LeftY, 'Visible','off');
    %     set(handles.leftytxt, 'Visible','off');
end

% if isempty(index2),index2={[];[];[]};end
clear ymin ymax
if vy2==0, vy2=Nonecol; end
set(get(handles.axes2,'XLabel'),'String','');
set(get(handles.axes2,'YLabel'),'String',char(DataA(1, vy2)));
set(handles.rightytxt, 'Visible','on');
set(handles.RightY,'Visible','on', 'Value', vy2);
if ~isempty(index2)
    if length(index2{1})+ length(index2{2})+ length(index2{3})>0
        %if data part of fco2 calc, check that calc has been done
        [tf,indx]=ismember(vy2,fco2c);
        if  tf & vfco2ok==0
            mmm=msgbox(sprintf('\"%s\" cannot be plotted...Data not calculated yet.\nUse the ''Calc fCO2'' button',char(DataA(1, fco2c(indx)))),'Plotting Error','error','modal');
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
            return;
        end
        axes(handles.axes2)
        set(handles.axes2,  'Visible','on', 'NextPlot','replace');
        if ~isobject(h)
            x2lim=get(handles.axes2,'XLim');
            y2lim=get(handles.axes2,'YLim');
        end
        %     if length(index2{1})+ length(index2{2})+ length(index2{3})==0
        %         if ~isempty(hp2),cla(handles.axes2,'reset');hp2=[];end
        %     end
        
        hp2=[];newlong2=cell(1,3);
        for i=1:3
            if length(index2{i})>0
                if (vx==vlong & get(handles.rdodateline,'value'))
                    newlong2{i}=str2num(char(DataA(index2{i}, vx)));                    
                    newlong2{1,i}(newlong2{:,i}<-3)=newlong2{1,i}(newlong2{:,i}<-3)+360;
                    hp2(i)=plot(newlong2{i},str2num(char(DataA(index2{i}, vy2))),'Parent',handles.axes2,...
                        'Marker','.', 'MarkerSize',dotsize,'MarkerEdgeColor',coul2(i,:), 'MarkerFaceColor',coul2(i,:),'LineStyle','none');
                else
                    hp2(i)=plot(str2num(char(DataA(index2{i}, vx))),str2num(char(DataA(index2{i}, vy2))),'Parent',handles.axes2,...
                        'Marker','.', 'MarkerSize',dotsize,'MarkerEdgeColor',coul2(i,:), 'MarkerFaceColor',coul2(i,:),'LineStyle','none');
                    yy=str2num(char(DataA(index2{i}, vy2)));
                    if ~exist('ymin','var'), ymin=min(yy);ymax=max(yy);else, ymin=min(ymin,min(yy));ymax=max(ymax,max(yy));end
                end
                set(hp2(i),'HitTest','off');
                if ~isobject(h)
                    set(handles.axes2,'XLim',x2lim,'YLim',y2lim);
                elseif ~isempty(index1{i})
                    x1lim=get(handles.axes1,'XLim');
                    x2lim=get(handles.axes2,'XLim');
                    set(handles.axes2,'XLim',[min(x1lim(1),x2lim(1)) max(x1lim(2),x2lim(2))]);
                    set(handles.axes1,'XLim',[min(x1lim(1),x2lim(1)) max(x1lim(2),x2lim(2))]);
                end
                set(handles.axes2,  'NextPlot','add');
            end
        end
        y2lim=handles.axes2.YLim;
        if abs(ymin-y2lim(1))< 0.1*(ymax-ymin)
            set(handles.axes2,'YLim',[ymin-0.1*(ymax-ymin) y2lim(2)]);%gives space at bottom of plot
        end
        set(handles.axes2,'Position',get(handles.axes1,'Position'), 'XColor','K','YColor','K',...
            'XAxisLocation','top', 'YAxisLocation','right', 'Xtick',[]);
        set(handles.axes2,'ButtonDownFcn',f2);
        set(handles.axes2,'HitTest','on');
        set(get(handles.axes2,'XLabel'),'String','');
        set(get(handles.axes2,'YLabel'),'String',char(DataA(1, vy2)));
        set(handles.rightytxt, 'Visible','on');
        set(handles.RightY,'Visible','on', 'Value', vy2);
        if isobject(h)
            listg=[(DataA(1, vy1));(DataA(1, vy2))];
            set(handles.popGselect, 'String',listg);
        end
    else
        set(get(handles.axes2,'XLabel'),'String','');
        set(get(handles.axes2,'YLabel'),'String',char(DataA(1, vy2)));
        set(handles.rightytxt, 'Visible','on');
        set(handles.RightY,'Visible','on', 'Value', vy2);
        if ~isempty(hp2), delete(hp2(hp2>0));hp2=[];end;
        set(handles.axes2,  'Visible','off');
        listg=[(DataA(1, vy1))]; %#ok<NBRAK>
        set(handles.popGselect, 'String',listg,'Value',1);
        %     set(handles.RightY, 'Visible','off');
        %     set(handles.rightytxt, 'Visible','off');
    end
else
    if ~isempty(hp2), delete(hp2(hp2>0));hp2=[];end;
    set(handles.axes2,  'Visible','off');
    listg=[(DataA(1, vy1))]; %#ok<NBRAK>
    set(handles.popGselect, 'String',listg,'Value',1);
    
end
lstGselect_Callback(handles.popGselect,1,handles,1);


% --------------------------------------------------------------------
function plotdata_single(handles,indexa, vx, vy1, colors, labels, h)
% Stub for Callback of the uicontrol handles.CurrentX.
global DataA  f1 hp hp2 motion  Nonecol rfh nrec;


axes(handles.axes2);
if isobject(h)
     cla(handles.axes2,'reset');
     hp2=[];clearvars -global rfh;nrec=0;
end
if  vx== Nonecol & h~=handles.btnplot
    mmm=msgbox(sprintf('\"X Data\" column set to \"None\".\nData needs to be assigned \nbefore it can be plotted.'),'Plotting Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
if   vy1== Nonecol & h~=handles.btnplot
    mmm=msgbox(sprintf('\"Left Y Data\" column set to \"None\".\nData needs to be assigned \nbefore it can be plotted.'),'Plotting Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
%deals with missing x values (y values taken care of in individual callbacks
indmvx=find(strcmp(DataA(:, vx),'-999')==0);

if double(get(handles.rdomv,'Value'))==0  %don't plot missing values
    for i=1:size(indexa,1)
        indexa{i}=intersect(indexa{i},indmvx);
    end
end
set(handles.axes2,  'Visible','off');

axes(handles.axes1);
c=get(handles.Graphc_s,'BackgroundColor');
set(handles.axes1,  'Visible','on', 'NextPlot','replace','Color',c);
if ~isobject(h)
    x1lim=get(handles.axes1,'XLim');
    y1lim=get(handles.axes1,'YLim');
end
if isobject(h)
     cla(handles.axes1,'reset');
     hp=[];clearvars -global rfh;nrec=0;
end
j=1;labels2=[];labels=char(labels);
for i=1:length(indexa);
    if length(indexa{i})>0
        hp(i)=plot(str2num(char(DataA(indexa{i}, vx))),str2num(char(DataA(indexa{i}, vy1))),colors{1}(i,:),'Parent',handles.axes1);
        yy=str2num(char(DataA(indexa{i}, vy1)));
        if ~exist('ymin','var'), ymin=min(yy);ymax=max(yy);else, ymin=min(ymin,min(yy));ymax=max(ymax,max(yy));end
        set(hp(i),'HitTest','off');
%         labels2{j}=labels(i,:);j=j+1;
        labels2{j}=labels(j,:);j=j+1;
        if ~isobject(h)
            set(handles.axes1,'XLim',x1lim,'YLim',y1lim);
        end
        set(handles.axes1,  'NextPlot','add','Color',c);
    end
end
y1lim=handles.axes1.YLim;
if abs(ymin-y1lim(1))< 0.1*(ymax-ymin)
    set(handles.axes1,'YLim',[ymin-0.1*(ymax-ymin) y1lim(2)]);
end

if size(labels2,2)==1
    legend(handles.axes1,'off');
else
    legend(handles.axes1,char(labels2),'Location','Best');
end

set(handles.axes1,'Position',get(handles.axes1,'Position'), 'XColor','K','YColor','K',...
    'XAxisLocation','bottom', 'YAxisLocation','left',  'Color',c);
set(handles.axes1,'ButtonDownFcn',f1);
set(handles.axes1,'HitTest','on');
set(get(handles.axes1,'XLabel'),'String',char(DataA(1, vx)))
set(get(handles.axes1,'YLabel'),'String',char(DataA(1, vy1)))
set(handles.CurrentX, 'Value', vx);
% set(handles.RightY, 'Visible','off');
% set(handles.rightytxt, 'Visible','off');
set(handles.LeftY, 'Value', vy1);
set(handles.RightY, 'Value', Nonecol);

if isobject(h)
    listg=(DataA(1, vy1));
     set(handles.popGselect, 'String',listg,'Value',1);
%    set(handles.popGselect, 'String',listg);
    if motion == 0,  motion=1;   end
end
lstGselect_Callback(handles.popGselect,1,handles,1);


function map_zoom_out_Callback(hObject, eventdata, handles)
global map_xlim map_ylim x1lim y1lim;

zout=get(hObject,'Value');
set(handles.axes1,'XLim',x1lim+(map_xlim-x1lim)*zout,'YLim',y1lim+(map_ylim-y1lim)*zout);

% --------------------------------------------------------------------
function  Flag_Color(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Flagc1.
c=uisetcolor([1 1 1],'Color');
set(h,'BackgroundColor',c);
set(h,'ForegroundColor','k');

% --------------------------------------------------------------------
function  Graph_Color(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Flagc1.
c=uisetcolor([1 1 1],'Color');
set(h,'BackgroundColor',c);set(h,'ForegroundColor','k');
lstGselect_Callback(handles.popGselect,1,handles,1);

% --- Executes on slider movement.
function Graphc_s_Callback(hObject, eventdata, handles)
% hObject    handle to Graphc_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
cv=get(hObject,'Value');c=[1-cv 1-cv 1-cv];
set(hObject,'BackgroundColor',c);set(hObject,'ForegroundColor','k');
lstGselect_Callback(handles.popGselect,1,handles,1);

% --------------------------------------------------------------------
function Change_flag(h, eventdata, handles, varargin)
global DataA ind1 ind2 hp hp2 x y Currentplot;
global sub_flag;
global vflag vsubf vsubfu vtype vlong newlong1 newlong2 vtin vtini;

% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';...
%         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';...
%         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';... %         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';... %         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'interpolated using questionable SST'};  variable missing';'interpolated using questionable SST'};

activity(1,handles); %drawnow;
equ_atm={'EQU';'ATM';'EQU-DRAIN';'ATM-DRAIN'};
equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};


gs=get(handles.popGselect,'String');
if strcmp(gs(1,:),'Graph 1')==1 & strcmp(gs(2,:),'Graph 2')==1
    mmm=msgbox('No data plotted','Point flagging','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    activity(0,handles);return;
end


ax_fl=get(h,'Tag');  %Tag is 'fl'+axe # + flag #
axstr=ax_fl(3);fl=ax_fl(4);
subflag=ax_fl(5:end);
if strcmp(axstr,'1')==1
    ax=handles.axes1;
    yax=handles.LeftY;
    ind=ind1;
    chp=hp;
    if exist('newlong1','var') newlong=newlong1; end
else
    ax=handles.axes2;
    yax=handles.RightY;
    ind=ind2;
    chp=hp2;
    if exist('newlong2','var') newlong=newlong2; end
end

xlimits = get(ax,'XLim');
ylimits = get(ax,'YLim');
xlimits1 = get(handles.axes1,'XLim');
ylimits1 = get(handles.axes1,'YLim');
xlimits2 = get(handles.axes2,'XLim');
ylimits2 = get(handles.axes2,'YLim');
xp = (xlimits(2)-xlimits(1))/100;
yp = (ylimits(2)-ylimits(1))/100;
type='';types='';
switch subflag
    case {'3','4','5','6','7'} % EqT, DT, SSS, pressures, low EQ flow
        type=equ;types='EQU';
    case {'0','1','2','9','10'}% no subflag, std range, SST, <3 std corr, Other
        type=equ_atm;types='EQU and ATM';
    case {'8'}% Q? Air value
        type=atm;types='ATM';
end

if get(yax,'Value')==vtin && str2num(fl)==4
    msg1='MEASURED SST cannot be flagged 4 due to the time offset interpolation.';
    msg2='Only INTERPOLATED SST (Tin interp) can be flagged 4.';
    msg3='Either INTERPOLATE or GENERATE from other source and flag 3 for Quest. SST';
    msg4='Or simply flag 3 for Quest. SST';
    mmm=msgbox(sprintf('%s\n%s\n%s\n%s',msg1,msg2,msg3,msg4),'Flag Issue','modal');CenterWindow(handles,mmm,handles.figure1) ; uiwait(mmm);
    return;
end
if get(yax,'Value')==vtini && str2num(fl)==3
    msg1='INTERPOLATED SST cannot be flagged 3 manually.';
    msg2='If MEASURED SST used to interpolate it is not questionable...';
    msg3='then simply flag 4.';
    mmm=msgbox(sprintf('%s\n%s\n%s',msg1,msg2,msg3),'Flag Issue','modal');CenterWindow(handles,mmm,handles.figure1) ; uiwait(mmm);
    return;
end

if ax_fl(1)~='r', return;end   %point only, not rectangle

    
%rectangle

rh=findobj(ax,'Type','Rectangle');
ptaallr=[];rm=0;%all rectangles
if isempty(rh)
    mmm=msgbox('No area selected','Area Flagging','error','modal');
    CenterWindow(handles,mmm,handles.figure1); uiwait(mmm);activity(0,handles);return;
end

for r=1:length(rh)% for each rectangle
    rv=get(rh(r),'Position');
    xv=[rv(1) rv(1) rv(1)+rv(3) rv(1)+rv(3)];
    yv=[rv(2) rv(2)+rv(4) rv(2)+rv(4) rv(2)];
    pta1all=[];pta2all=[];ptaall=[];
    for i =1:length(chp)
        if (handles.CurrentX.Value==vlong & get(handles.rdodateline,'value'))
            pta1{i} = inpolygon(newlong{i},str2num(char(DataA(ind{i},get(yax,'Value')))),xv,yv);
        else
            pta1{i} = inpolygon(str2num(char(DataA(ind{i},handles.CurrentX.Value))),str2num(char(DataA(ind{i},get(yax,'Value')))),xv,yv);
        end
        %pta1 are points in polygon, pta2 are points of right type
        pta1{i}=find(pta1{i});pta1all=union(pta1all,ind{i}(pta1{i}));
        pta2{i}=[];
        if str2num(fl)==3 & ~strcmp(type,''),pta2{i}=find(ismember(DataA(ind{i},vtype),type));
        else pta2{i}=pta1{i};end
        pta2all=union(pta2all,ind{i}(pta2{i}));
        pta{i}=intersect(pta1{i},pta2{i});
    end
    % ptaall are points of right type in rectangle
    ptaall=intersect(pta1all,pta2all);
    % ptaallr are points of right type in all rectangles drawn
    ptaallr=union(ptaallr,ptaall);
    
    %for fl=3 if (all or some) points are of the wrong type (pta1all>ptaall) --> message
    %for fl=4, can flag anything
    if str2num(fl)==3 & ~isempty(pta1all) & (size(pta1all,1)>size(ptaall,1)), rm=1;end
end %end of each rectangle

ptaallr=unique(ptaallr); %to not select same points several times in case rectangles intersect each other
%if points of the wrong type were selected: message
if rm==1, mmm=msgbox(sprintf('Only %s will be flagged for this parameter.',types),'Flag Issue','modal');CenterWindow(handles,mmm,handles.figure1) ; uiwait(mmm);end


if str2num(fl)~=3 %flags 2 and 4
    DataA(ptaallr,vflag)={fl};
    DataA(ptaallr,vsubf)=cellstr('');
    DataA(ptaallr,vsubfu)=cellstr('');
else
    if str2num(subflag)==0, Set_subflag(ptaallr,str2num('-1'),handles);  % Set_subflag sets flag to 4 if subflag=0
    elseif str2num(subflag)<=10,Set_subflag(ptaallr,str2num(subflag),handles);
    elseif str2num(subflag)>10, Set_subflag(ptaallr,10,handles);Set_subflaguser(ptaallr,str2num(subflag),handles);
    end
end

feval(Currentplot,1,1,handles,1);
set(handles.axes1,'XLim',xlimits1);set(handles.axes1,'YLim',ylimits1);
set(handles.axes2,'XLim',xlimits2);set(handles.axes2,'YLim',ylimits2);



clear pti1all pti2all ptiall pti1 pti2 pti;
clear pta1all pta2all ptaall pta1 pta2 pta ptaallr;

Stat_Update_Callback(0,0,handles);
Check_for_999(handles);
activity(0,handles);


% --------------------------------------------------------------------
function Reset_flags(h, eventdata, handles, varargin)
global DataA Currentplot;
global rangeok xco2corrok vfco2ok airiok;
global vflag vsubf vsubfu ;

% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';...
%         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';...
%         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';... %         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';... %         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'interpolated using questionable SST'};  variable missing';'interpolated using questionable SST'};

activity(1,handles); %drawnow;

gs=get(handles.popGselect,'String');
if strcmp(gs(1,:),'Graph 1')==1 & strcmp(gs(2,:),'Graph 2')==1
    mmm=msgbox('No data plotted','Point flagging','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
else
    DataA(2:end,vflag)={'2'};
    DataA(2:end,vsubf)=cellstr('');
    DataA(2:end,vsubfu)=cellstr('');
    feval(Currentplot,1,1,handles,1);
end
rangeok=0;xco2corrok=0;vfco2ok=0; airiok=0;
set(handles.air_interpolation,'Visible','off');
set(handles.correctxco2,'Visible','off');
set(handles.cfco2,'Visible','off');
% set(handles.rangeok,'Visible','off','BackgroundColor','r');
set(handles.stdok,'visible','off','BackgroundColor','r');
% set(handles.xco2rangeok,'visible','off','BackgroundColor','r');
set(handles.xco2aok,'visible','off','BackgroundColor','r');

% set(handles.despiket,'String', 'Check Range - Not xCO2');
% set(handles.rangeokt,'String','Range Check - Not OK');
set(handles.stdokt,'String','Standards - Not Flagged');
% set(handles.xco2rangeokt,'String','xCO2 Range - Not OK');
set(handles.xco2aokt,'String','xCO2(atm) - Not Flagged');
Check_for_999(handles);


function Interpolate(h, eventdata, handles, varargin)
global DataA ind1 ind2 hp hp2 Currentplot dotsize monitors;
global   vteq vtin vtini vsal vpeqa vpatm vpamb listf_file sub_flag;

% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';...
%         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';...
%         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';... %         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';... %         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'interpolated using questionable SST'};  variable missing';'interpolated using questionable SST'};

activity(1,handles); %drawnow;

mmm=0;interpok=0;
gs=get(handles.popGselect,'String');
if strcmp(gs(1,:),'Graph 1')==1 & strcmp(gs(2,:),'Graph 2')==1
    mmm=msgbox('No data plotted','Interpolation','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

ax_fl=get(h,'Tag');  %Tag is 'ip'+axe #+flag to interpolate (3, 4 or 7(3+4))
axstr=ax_fl(3);
%     fl=str2num(ax_fl(4));
if strcmp(axstr,'1')==1
    ax=handles.axes1;
    yax=handles.LeftY;
    ind=ind1;
    chp=hp;
else
    ax=handles.axes2;
    yax=handles.RightY;
    ind=ind2;
    chp=hp2;
end
xlimits = get(ax,'XLim');
ylimits = get(ax,'YLim');
xp = (xlimits(2)-xlimits(1))/100;
yp = (ylimits(2)-ylimits(1))/100;

%     if get(yax,'Value')~=vteq & get(yax,'Value')~=vtin & get(yax,'Value')~=vlicorcav & get(yax,'Value')~=vsal
if get(yax,'Value')>length(listf_file)
    mmm=msgbox(sprintf('Can Only Interpolate Original Data.\n\n%s is Calculated.\n\nAction Cancelled',char(DataA(1,get(yax,'Value')))),'Interpolation Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
rh=findobj(ax,'Type','Rectangle');
ptaallr=[];rm=0;%all rectangles
if isempty(rh)
    mmm=msgbox('No area selected','Area Interpolation','error','modal');
    CenterWindow(handles,mmm,handles.figure1);uiwait(mmm);activity(0,handles);return;
end

for r=1:length(rh)% for each rectangle
    ptatoi=[];
    rv=get(rh(r),'Position');
    xv=[rv(1) rv(1) rv(1)+rv(3) rv(1)+rv(3)];
    yv=[rv(2) rv(2)+rv(4) rv(2)+rv(4) rv(2)];
    %index of points to interp
    ptatoi = inpolygon(str2num(char(DataA(2:end,handles.CurrentX.Value))),str2num(char(DataA(2:end,get(yax,'Value')))),xv,yv);
    ptatoi=find(ptatoi)+1;
    indall=[];
    for i=1:size(ind,1), indall=union(indall,ind{i});end % take only plotted points
    ptatoi=intersect(ptatoi,indall);clear indall;

    % ptaallr are points to interp in all rectangles drawn
    ptaallr=union(ptaallr,ptatoi);

end

if mmm==0  %if area selected, interpolates data
    xout=str2num(char(DataA(ptaallr,handles.CurrentX.Value)));
    yorig=str2num(char(DataA(ptaallr,get(yax,'Value'))));
    % use good point (ind{1} are flag2) less those to interp
    xin=str2num(char(DataA(setdiff(ind{1},ptaallr),handles.CurrentX.Value)));
    yin=str2num(char(DataA(setdiff(ind{1},ptaallr),get(yax,'Value'))));
    yin=yin(xin~=-999);xin=xin(xin~=-999);%gets rid of data with YDay=-999
    if length(xin)~=0 & length(yin)~=0 & length(xout)~=0
        interpok=1;
        try
            yout=interp1(xin,yin,xout,'linear');
        catch
            mmm=msgbox(lasterr,'Error');
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
            return;
        end
        if ~isempty(find(isnan(yout), 1))
            howmany='ALL';
            if length(find(isnan(yout)))<length(ptaallr),howmany=sprintf('%d',length(find(isnan(yout))));end
            mmm=msgbox(sprintf(['%s Points Were Outside the Range of Flag 2 Data\n\nAnd Were Not Extrapolated.'...
                '\n\nPoints Should Be Bracketed by Flag 2 Data to Be Interpolated.\n\nConsider flagging the %s points ''bad''.'],howmany,howmany),'Error');
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
            if strcmp(howmany,'ALL'), return;end
        end
        figi=figure;figip=get(figi,'position');set(figi,'position',[figip(1) figip(2) 0.5*monitors(1,3) 0.5*monitors(1,4)]);CenterWindow(handles,figi,handles.figure1);
        plot(xin,yin,'ok',xout(~isnan(yout)),yout(~isnan(yout)),'or',xout(~isnan(yout)),yorig(~isnan(yout)),'ob','MarkerSize',dotsize);
        set(gca,'XLim',xlimits);set(gca,'YLim',ylimits);
        legend('Good Data','Interpolated Data','Problematic Data','Location','Best');
        axes(gca);
        
        set(gcf,'Toolbar','figure')
        axl=0.15;axb=0.2;axw=0.75;axh=0.75;
        gcfp=get(gcf,'position'); set(gca,'position',[axl axb axw axh]);
        fl=gcfp(1);fb=gcfp(2);fw=gcfp(3);fh=gcfp(4);
        butw=fw/6;buth=0.1*fh;sp=fw/20;clear gcbf;
%         bt(1)=uicontrol('Position',[(fw-butw)/2 2 butw buth],'String','Continue...',...
%             'Callback','uiresume');
        bt(1)=uicontrol('Style','text','Position',[0 -10 butw buth],'String','Select: ');
        bt(2)=uicontrol('Style','radiobutton','Position',[butw 2 butw buth],'String','Auto Flag',...
            'Callback','uiresume','Value',0);
        bt(3)=uicontrol('Style','radiobutton','Position',[2*butw 2 butw buth],'String','Manual Flag',...
            'Callback','uiresume','Value',0);
        bt(4)=uicontrol('Style','radiobutton','Position',[3*butw 2 butw buth],'String','No Flag',...
            'Callback','uiresume','Value',0);
        bt(5)=uicontrol('Style','radiobutton','Position',[4*butw 2 butw buth],'String','Cancel',...
            'Callback','uiresume','Value',0);
        uiwait(figi);
        mmm=find(cell2mat(get(bt,'Value')));
%         mmm=questdlg(sprintf('Accept Interpolation?'),'Interpolation Results','Yes','No','No');
        close(figi);
        
        if mmm<5  %not Cancel
            DataA(ptaallr(~isnan(yout)),get(yax,'Value'))=cellstr(num2str(yout(~isnan(yout))));
            
            tf=ismember([vtin vtini vteq vsal vpeqa vpatm vpamb],get(yax,'Value'));
            if sum(tf)>0
%             mmm=questdlg(sprintf('FLAG DATA?\n\nAUTO = program does it.\nMANUAL = User chooses'),'Interpolation Flag','AUTO (3)', 'MANUAL', 'NO FLAG','AUTO (3)');
              if mmm==2  % strcmp(mmm,'AUTO (3)')==1
                  
                    tf=ismember([vtin vtini vteq vsal],get(yax,'Value'));
                    if sum(tf)>0
                        if sum(tf)>1
                            mmm=questdlg(sprintf(['These Data Are Assigned to at Least 2 of the Following Variables:'...
                                '\n\n\"in situ Temp\"\n\"EQU Temp\"\n\"Salinity\"\n\nHow Do You Want to Flag Them?'...
                                '\n\n2 - %s\n3 - %s\n5 - %s'],sub_flag{2,:},sub_flag{3,:},sub_flag{5,:}),'Interpolation Flag','2', '3', '5', '2');
                            flv=int8(str2num(mmm));
                        else
                            flv=2*(ismember(get(yax,'Value'),[vtin vtini])) + 3*(get(yax,'Value')==vteq)+5*(get(yax,'Value')==vsal);
                        end
                        Set_subflag(ptaallr(~isnan(yout)),flv,handles);
                    end
                    
                    tf=ismember([vpeqa vpatm vpamb],get(yax,'Value'));
                    if sum(tf)>0
                        flv=6; %Questionable/Interpolated P
                        Set_subflag(ptaallr(~isnan(yout)),flv,handles);
                    end                   
                elseif mmm==3 % strcmp(mmm,'MANUAL')==1
                    %--------------------------------------------------------------------------
                    % Creates interface
                    %--------------------------------------------------------------------------
                    [scrsz]=get(0,'ScreenSize');   wd=400;ht=200;spc=50;
                    
                    fig2 = figure('position', [(scrsz(3)-wd)/2 (scrsz(4)-ht)/2 wd ht],'menubar','none','numbertitle','off');
                    uicontrol ('Style','text','position', [(wd-200)/2 ht-30 200 20],'String','Select Subflag to Erase','FontWeight','bold','backgroundcolor',get(fig2,'color'));
                    r(1) = uicontrol ('Style','popupmenu','position', [(wd-250)/2 ht-60 250 20],'String',sub_flag,'backgroundcolor',get(fig2,'color'));
                    ccb=uicontrol ('Style','togglebutton','position', [(wd-100) 20 40 20],'String','Cancel','Callback','uiresume(gcbf)');
                    okb=uicontrol ('Style','togglebutton','position', [(wd-50) 20 30 20],'String','OK','Callback','uiresume(gcbf)');
                    uicontrol(okb);
                    uiwait(fig2);
                    
                    %set(fig2,'visible','off');
                    if ~get(okb,'value'), return;end
                    flv=min(get(r(1),'Value'),10);
                    if flv>10,  Set_subflag(ptaallr(~isnan(yout)),10,handles);    end
                    Set_subflag(ptaallr(~isnan(yout)),flv,handles);
                                        
                    close(fig2);
                end
                
            end
            Calc_Fields_Callback(handles,ptaallr(~isnan(yout)));
            Stat_Update_Callback(1, 1, handles);
            Check_for_999(handles);
            feval(Currentplot,1,1,handles,1);
        end
    end
end
if interpok==0
    mmm=msgbox (sprintf('No Data to interpolate.\n\nCheck Data or Flags'),'Interpolation Error!','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

activity(0,handles);

            


% --------------------------------------------------------------------
function Erase_Areas(h, eventdata, handles, varargin)
global rfh nrec;

mmm=questdlg('Erase last or All rectangles?','Erase','Last','All','All');

ax_fl=get(h,'Tag');  %Tag is 'fl'+axe # + flag #
axstr=ax_fl(3);
if strcmp(axstr,'1')==1
    ax=handles.axes1;
else
    ax=handles.axes2;
end
rh=findobj(ax,'Type','Rectangle');
if ~isempty(rh)
    if strcmp(mmm,'All')>0 | size(rh,1)==1
        delete(rh); nrec=0;clearvars -global rfh;
    else
        delete(rh(rh==rfh(end)));rfh=rfh(1:end-1);
    end
else
    mmm=msgbox('No Area found','Area Deletion','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
end

% --------------------------------------------------------------------
function btnSave_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.btnSave.
global DataA nvar Nonecol erreur datap resultp where_xml; %#ok<NUSED>
global expot expod vgroup vship vcruiseid dryxco2 sub_flag;
global vyday vdutc vtutc vlat vlong vpeqa vteq   vtype  vsal  vtini  vpatm  vflag vsubf vsubfu...
      vxco2icorr  vlicorxw vlicorxa   vfco2w     vfco2i  vdfco2 ;
global headp hfile sysinip sysinif; %#ok<NUSED>
global vxaxiscok posok stdok atmok airiok vfco2ok dateok xco2corrok vtiniok listf_file; %#ok<NUSED>
global stdv usestdv GSID baroh equpdiff fversionn; %#ok<NUSED>
global hdllohistr hdltype vxaxis; %#ok<NUSED>
global indf1 indf2 indf saveok finalok workingok wfn ffn pi_names save_opts OTin;

equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};
equ_atm=cat(1,equ,atm);

saveok=0;

Savef_v2(handles.figure1);

if ~saveok, return;end
% if ~dateok, return;end

activity(1,handles); 

    
Save_Ini(datap,resultp,headp,hfile,sysinip,sysinif,where_xml);

headp0=headp;hfile0=hfile;sysinip0=sysinip;sysinif0=sysinif;
headp=resultp;sysinip=resultp;hfile=[vcruiseid ' Data Config_v' strtrim(num2str(fversionn)) '.csv'];sysinif=[vcruiseid ' System Config_v' strtrim(num2str(fversionn)) '.ini'];
vgroup=strtrim(vgroup);vship=strtrim(vship);pi_names=strtrim(pi_names);vcruiseid=strtrim(vcruiseid);
Config_Sys_Save(-2, -1, handles); 
Config_Data_Save(-2,-1,handles);
    
% *************************************************************************
%if ispc, ext2='.xls';else ext2='.csv';end
ext2='.csv';

 
if finalok
%--------------------------------    % save xls (or csv for mac) final data --------------------------------------------------------------------------

    %puts DataA in DataO for FINAL File
    %Add columns you want to output in outvar
    outtit={'YD_UTC'; 'DATE_UTC__ddmmyyyy'; 'TIME_UTC_hh:mm:ss'; 'LAT_dec_degree'; 'LONG_dec_degree';...
        'xCO2_EQU_ppm'; 'xCO2_ATM_ppm'; 'xCO2_ATM_interpolated_ppm'; 'PRES_EQU_hPa';'PRES_ATM@SSP_hPa';'TEMP_EQU_C';...
        'SST_C'; 'SAL_permil';'fCO2_SW@SST_uatm'; 'fCO2_ATM_interpolated_uatm'; 'dfCO2_uatm';...
        'WOCE_QC_FLAG'; 'QC_SUBFLAG'};
    outvar=[vyday vdutc vtutc vlat vlong vlicorxw vlicorxa vxco2icorr ...
        vpeqa vpatm vteq vtini vsal vfco2w vfco2i vdfco2 vflag vsubf];
    clear DataO;%drawnow;


         height=str2num(get(handles.dbaroh,'string'));
         if isempty(height), height=0;end
         if height<0 | height>100, height=0; end
         set(handles.dbaroh,'String',strtrim(num2str((height))));

         DataO=cell(length(indf)+1,length(outvar));

%          %Add columns for Expocode, group, ship
          extratit={'Expocode';'Group';'Ship'};j=1;
          for i=1:3
             if save_opts(i)
                 DataO(1,j)=extratit(i);
                 if i==1, DataO(2:end,j)={strcat(expot,expod)};
                 elseif i==2,DataO(2:end,j)={vgroup};
                 elseif i==3,DataO(2:end,j)={vship};end
                 j=j+1;
             end
          end
          datao_shift=sum(save_opts(1:3));
          %DataO(2:end,1)={strcat(expot,expod)};DataO(2:end,2)={strcat(vgroup,'_',vship)}; DataO(2:end,3)={vcruiseid};
          %datao_shift=0; %expocode, pi names, ship, group added as lines of text before data for SOCAt dashboard
         

          for i=1:length(outtit)
             if ismember(outvar(i),[vteq vpeqa vpatm vtini vsal]) %all teq, peq, patm,tin, sal for ATM measurements are set to -999.
                 junk=DataA(:,outvar(i));
                 junk(ismember(DataA(:,vtype),atm))={'-999'};
                 DataO(2:end,i+datao_shift)=junk(indf);
                 
                 if outvar(i)==vpatm  %Patm corrected to sea level
%                  junk=DataA(:,outvar(i));
                     indf3=find(ismember(DataA(:,vtype),equ));
                     indf4=intersect(indf3,indf2);
                     patm=str2num(char(junk(indf4)));
                     tk=str2num(char(DataA(indf4,vteq)))+273.15;
                     for j= 1:length(indf4),junk(indf4(j))={num2str(patm(j)+((patm(j)./tk(j)./(8.314/28.97e-3))*9.8*height))};end
                     DataO(2:end,i+datao_shift)=junk(indf);
                 end
             elseif OTin~=0 %if OTin not 0, SST is interpolated and vflag has to be modified
                 if outvar(i)==vflag
                     junk=DataA(:,vflag);
                     junksubf=DataA(:,vsubf);
                     %if ONLY flagged quest. SST, resets flag to 2 and erase subflag text
                     ssearch=regexprep(char(sub_flag(2)),'(','\\('); ssearch=regexprep(ssearch,')','\\)'); ssearch=regexprep(ssearch,'+','\\+');
                     qssti=regexp(DataA(:,vsubf),ssearch);%find quest. sst subflag
                     qssti=~cellfun('isempty',qssti);%0 and 1s
                     onlysst=ismember(DataA(:,vsubf),sub_flag(2));%0s and 1s
                     junk(onlysst)={'2'};
                     junksubf(onlysst)={''};
                     %remove 'quest.sst;' if not alone
                     junksubf(~onlysst & qssti)=cellstr(strcat(regexprep(junksubf(~onlysst & qssti), strcat(sub_flag(2),';'),'')));
                     
                     %if contains 'quest.interp.SST'
                     %if no other flag>10, replace 'other' with 'quest.interp. sst'
                     %else add 'quest. interp. SST'
                     hasother=zeros(size(junk,1),1);
                     for k=11:size(sub_flag,1)
                         ssearch=regexprep(char(sub_flag(k)),'(','\\('); ssearch=regexprep(ssearch,')','\\)'); ssearch=regexprep(ssearch,'+','\\+');
                         oneother=regexp(DataA(:,vsubfu),ssearch);%find  subflag
                         oneother=~cellfun('isempty',oneother);%0 and 1s
                         if k==17
                             qissti=oneother;
                         else
                             hasother=hasother | oneother;
                         end
                     end
                     
                     junksubf(~hasother & qissti)=cellstr(regexprep(junksubf(~hasother & qissti),sub_flag{10},sub_flag{17}));
                     junksubf(hasother & qissti)=cellstr(strcat(char(junksubf(hasother & qissti)), ';', sub_flag(17)));
                     
                     DataO(2:end,i+datao_shift)=junk(indf);
                     j=find(outvar==vsubf);
                     DataO(2:end,j+datao_shift)=junksubf(indf);
                 else
                     %vsubf taken care of just above
                     if  outvar(i)~=vsubf,DataO(2:end,i+datao_shift)=DataA(indf,outvar(i));end
                 end
             else % if OTin==0 all other columns unmodified
                 DataO(2:end,i+datao_shift)=DataA(indf,outvar(i));
             end
             DataO(1,i+datao_shift)=outtit(i);
         end
         datsave=(DataA(indf,vdutc));[ysave, msave, dsave, ~, ~, ~]=datevec(datsave,'mm/dd/yy');
%         ysave=year(datsave,'mm/dd/yy');[msave,junk]=month(datsave,'mm/dd/yy');dsave=day(datsave,'mm/dd/yy');
         DataO(2:end,2+datao_shift)= cellstr([num2str(dsave,'%02i') num2str(msave,'%02i') num2str(ysave,'%04i')]);
         
         clear indf1 indf2 indf3 indf4 indf indfatm junk patm tk datsave ysave msave dsave;
         DataO(1,:)=regexprep(DataO(1,:),'\%','\%\%');%replaces special characters % to export in .csv
         DataO(strcmp(DataO,'NaN')==1)={'-999'};%Replace  NaN by -999 (bad data)
         DataO=strtrim(DataO);

         %filename=[char(sp) filesep  char(sf) exts];
       
         [erreur]=cell2csv(handles,ffn,DataO);
         erreur=~erreur;
         
         if erreur==0
             mmm1=errordlg(sprintf('Error saving %s .',ffn),'file write error');  uiwait(mmm1);
         end
    
    clear DataO;%drawnow;
end


%--------------------------------    % save working data -------------------------------------------------------------------------------
    %exts=['_Working.mat'];
    %if sum(s)==2, msgbox(mess(3,:));  uiwait; end

    fversions=fversionn;    
    DataO=DataA;

    DataO(1,:)=regexprep(DataO(1,:),'\%','\%\%');%replaces special characters % to export in .xls

    lba=find(strcmp(DataO,'NaN')==1); DataO(lba)={'-999'}; 
    DataO=strtrim(DataO);
    clear lba;clear DataA;%drawnow;
 
%--------------------------------    %save mat working data     --------------------------------------------------------------------------
    matfn=regexprep(wfn,'_Working.csv','_Working.mat');
    h = waitbar(0,'\fontname{courier}Saving mat File...');
    CenterWindow(handles,h,handles.figure1);
    OTin=str2num(get(handles.AOvalue,'String'));

    hdltype=[];
    for i=1:10 %checks which type to check for outliers - atm, equ or both    % -or=out of range   -ol=outliers
        hdllohistr{i,1}=strtrim(handles.("lo" + num2str(i)).String);
        hdllohistr{i,2}=strtrim(handles.("hi" + num2str(i)).String);

        switch 10*get(handles.("atm" + num2str(i)),'Value') + get(handles.("equ" + num2str(i)),'Value')
            case 11
                hdltype{i}=equ_atm;
            case 10
                hdltype{i}=atm;
            case 1
                hdltype{i}=equ;
            case 0
                hdltype{i}=equ_atm;
        end
    end
    
    save(matfn , '-v7', 'DataO', 'nvar', 'Nonecol', 'vxaxiscok', 'posok', 'stdok','atmok','airiok','vfco2ok', 'dateok', 'xco2corrok',...
        'listf_file', 'datap', 'resultp', 'headp', 'hfile', 'sysinip', 'sysinif','stdv','usestdv','GSID','baroh','equpdiff', ...
        'hdllohistr', 'hdltype', 'fversions', 'vgroup', 'vship', 'pi_names','vcruiseid','OTin','vtiniok','vxaxis','dryxco2','expot','expod',...
        'wfn','ffn','save_opts');
    waitbar(1,h);
    close(h);
    %restore default config files
    headp=headp0;hfile=hfile0;sysinip=sysinip0;sysinif=sysinif0;

%--------------------------------    %end save mat working data     --------------------------------------------------------------------------

if workingok
%--------------------------------    %save csv working data     --------------------------------------------------------------------------
    
        [erreur]=cell2csv(handles,wfn,DataO);
        erreur=~erreur;

        if erreur==0
        mmm1=errordlg(sprintf('Error saving %s file.',wfn),'file write error');  uiwait(mmm1);
%         if error occurred, will try to save the file by splitting it. Each part is saved on a different 
%              sheet in the xls file. The sheets are named 'PART 1', 'PART 2', etc...
            options.Interpreter='tex';options.Resize='off';
            prompt = sprintf('%s\n%s\n%s\n%s','\fontname{Times New Roman}\fontsize{10}Data will be split','to solve memory problem.',...
                'Enter HOW MANY PARTS','you want to split the data into{\bf\color{red} (1-10)}:');
            dlg_title = 'Memory Problem While Saving Data';
            num_lines = 1.3;
            def = {'3'};mmmn=-1;
            while ~isempty(def) & (mmmn<1 | mmmn>10)
                def = inputdlg(prompt,dlg_title,num_lines,def,options);
                mmmn=str2num(char(def));
            end
            if ~isempty(def)
                partfmax=mmmn;
                p1=0;
                for partf=1:partfmax
                    p0=p1;p1=partf*round(size(DataO,1)/(partfmax));
                    if partf==partfmax, p1=size(DataO,1);end
                    [erreur]=cell2csv(handles,regexprep(wfn,'.csv', ['_' sprintf('Part %d',partf) '.csv']),DataO(1+p0:p1,:));
                    erreur=~erreur;
                end
            end
    end
    
    clear DataO;

end
%--------------------------------    %end save csv working data     --------------------------------------------------------------------------


  
activity(0,handles); 
   


% --------------------------------------------------------------------
function  fcmotion(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.frame6.
global motion ;

if motion~=0
    if strcmp(get(handles.motiontxt,'String'),'Cursor INACTIVE')==1
        %      zoom (handles.figure1,'on');
        motion=2;
        set(handles.motiontxt,'String','Cursor ACTIVE','TooltipString','Cursor ACTIVE')
        handles.axes1.Toolbar.Visible = 'off';
        handles.axes2.Toolbar.Visible = 'off';
       %to undo any tools user might have selected
        pan(handles.figure1,'on');
        pan(handles.figure1,'off');
        set(h,'BackgroundColor',[1    0    0]);
    else
        %      zoom (handles.figure1,'off');
        motion=1;
        set(handles.motiontxt,'String','Cursor INACTIVE','TooltipString','Cursor INACTIVE')
        handles.axes1.Toolbar.Visible = 'on';
        handles.axes2.Toolbar.Visible = 'on';
        set(h,'BackgroundColor',[0.9255    0.9137    0.8471]);
    end
else
    mmm=msgbox('No points plotted yet!','Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
end

%% Correct xCO2 for Standards

% --------------------------------------------------------------------
function  correctxco2_Callback(h, eventdata, handles, varargin) %#ok<DEFNU>
% Stub for Callback of the uicontrol handles.correctxco2.

global DataA DSize xco2corrok;
global vxaxis vtype vstdx  vxco2dry vlicorxcorr vflag vstdi;
global dlicorxcorr fxco2w fxco2a fxco2u vlicorxw vlicorxa;
global siw sia siu;
global stdv stdval dstd dotsize;

activity(1,handles); 
%drawnow;
DSize=size(DataA,1);
h = waitbar(0,'\fontname{technical}Initializing...');
CenterWindow(handles,h,handles.figure1);
dtype=DataA(:,vtype);
total=4;
A=str2num(char(DataA(2:end,vxaxis)));A(DSize)=0;
dyday=circshift(A,1) ;
waitbar(1/total,h);
A=(regexprep(DataA(2:end,vstdx),'[^0-9.+-]','0'));
A=str2num(char(A));A(DSize)=0; %keeps only numbers (for discrete pCO2)
dstdx=circshift(A,1) ;
waitbar(2/total,h);
A=str2num(char(DataA(2:end,vxco2dry)));A(DSize)=0;
dlicorx=circshift(A,1) ;
waitbar(3/total,h);
A=str2num(char(DataA(2:end,vflag)));A(DSize)=0;
dflag=circshift(A,1) ;
waitbar(4/total,h);
close(h);clear A;

if verLessThan('matlab','9.1')
    fxco2w=find(~cellfun(@isempty,strfind(DataA(:,vtype),'EQU')));
    fxco2a=find(~cellfun(@isempty,strfind(DataA(:,vtype),'ATM')));
    fxco2u=find(~cellfun(@isempty,strfind(DataA(:,vtype),'UNK')));
else
    fxco2w=find(contains(DataA(:,vtype),'EQU'));
    fxco2a=find(contains(DataA(:,vtype),'ATM'));
    fxco2u=find(contains(DataA(:,vtype),'UNK'));
end



if verLessThan('matlab','9.1')
    siw=find(~cellfun(@isempty,strfind(dtype(2:end),'EQU')) & dyday(2:end)~=-999);
    sia=find(~cellfun(@isempty,strfind(dtype(2:end),'ATM')) & dyday(2:end)~=-999);
    siu=find(~cellfun(@isempty,strfind(dtype(2:end),'UNK')) & dyday(2:end)~=-999);
else
    siw=find(contains(dtype(2:end),'EQU') & dyday(2:end)~=-999);
    sia=find(contains(dtype(2:end),'ATM') & dyday(2:end)~=-999);
    siu=find(contains(dtype(2:end),'UNK') & dyday(2:end)~=-999);
end
si=union(union(siw,sia),siu);
%stdig = index of STD in groups of at least 2. Single STD not taken into account.
%finds STDx which has another STDx before or after it, thus part of a group
%regexp(DataA(1:15,vtype),'^STD\d\d?(?![sz])') will find any STD#(#) except if z or s

%get average time between a standard and the standard right before it...
cond1=(~cellfun(@isempty,regexp(dtype(2:end),'^STD\d\d?(?![sz])'))); %is a std but no zero/span
cond1(DSize)=0; cond1=circshift(cond1,1) ;

%cond1after=circshift(cond1,-1) ; cond1after(end)=0; %std which has a std after
cond1before=circshift(cond1,1) ; cond1before(1)=0;cond1before=(cond1&cond1before);%std which has a std before
avgtimebefore=(dyday(find(cond1before))-dyday(find(cond1before)-1))*24*60;%all times between std and std before
avgtimebefore=mean(avgtimebefore(avgtimebefore<10));%average of all these times

%cond1=(~cellfun(@isempty,regexp(dtype(2:end),'^STD\d\d?(?![sz])'))); %is a std but no zero/span
cond2=(dyday(3:end)-dyday(2:end-1))*24*60<=avgtimebefore+1 ;% meas before is less than x min before (1st element=0)
cond2(DSize-1:DSize)=0;% when shifted once, cond2 is actually array of elements whose next meas is less than x min apart (last element = 0)
cond2=circshift(cond2,1);% array whose meas AFTER is < x min apart
cond3=circshift(cond2,1);% array whose meas BEFORE is < x min apart

stdig=find(cond1 & (cond2|cond3)); %indices of STDs having another STD within +/- x min.

if isempty(stdv) || strcmp(get(handles.stdo_f,'BackgroundColor'),'r')
    strmsg=sprintf('Problems with STD values\nCheck data!');
    mmm=msgbox(strmsg,'Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

% if get(handles.rdostdv,'Value')
%     stdv=get(handles.tblstdv,'Data');
%     tts=''; for i=1:length(stdv), tts= [tts 'STD ' num2str(i) sprintf(' : %0.2f ppm',stdv(i))];if i<length(stdv), tts=[tts sprintf('\n')]; end, end
%     set(handles.tblstdv,'TooltipString',tts);
% else %takes average of each STD value from the 'stdvalue' column 
%     mmm=-1;stdv=[];
%     maxstd=max(cellfun(@str2num,strrep(strrep(dtype(find(cond1)+1),'STD',''),'-DRAIN',''))); %get max STD number
%     cc='g';mess='Std Offset OK';
%     for i=1:maxstd
%         stdx=(~cellfun(@isempty,regexp(dtype(2:end),['^STD' num2str(i) '(?![sz])'])));
%         if std(dstdx(find(stdx)+1))>1e-2
%             strmsg=sprintf('Error in %s known values.\nStandard Deviation > 0.01\nCheck data!', ['STD' num2str(i)]);
%             mmm=msgbox(strmsg,'Error','error','modal');
%             CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
%             cc='r';mess=strmsg; 
%             set(handles.stdo_f,'BackgroundColor',cc,'TooltipString',mess);
%         else
%             if ~isempty(stdx)
%                 stdv=[stdv;round(mean(dstdx(find(stdx)+1)),2)];
%             else
%                 stdv=[stdv;-999];
%             end
%         end
%         if mmm~=-1,return;end
%     end
% end

%creates stdix= index of STDx (or STDx-DRAIN) part of a group
%also finds the size of each stdix array and store in stdsize array.
%stdi1=index of STD1 in groups= index of STD1 that will be used for interpolation.
stda='';
for x=1:length(stdv)
    %eval(['stdigg' num2str(x) '=find((strcmp(dtype(stdig),' char(39) 'STD' num2str(x)  char(39) ')==1 | strcmp(dtype(stdig),'  char(39) 'STD' num2str(x) '-DRAIN'  char(39) ')==1) & dyday(stdig)~=-999);']);
    eval(['stdigg' num2str(x) '=find(~cellfun(@isempty,regexp(dtype(stdig),''^STD' num2str(x) '(?![sz])'')) & dyday(stdig)~=-999);']);
    %st=find(~cellfun(@isempty,regexp(dtype(stdig),'^STD\d\d?(?![sz])')) & dyday(stdig)~=-999);
    eval(['stdi' num2str(x) '=stdig(stdigg' num2str(x) ');']);
    if x==1, stdsize=length(eval(['stdi' num2str(x)])); else stdsize=[stdsize;length(eval(['stdi' num2str(x)]))];end 
end

%if one std has been measured less than 50% of the time of the other standards, it is not used.
%creates 2 arrays: std_incl= included standards and std_excl
%creates stdi={stdi1;stdi2...}
stdprop=stdsize/max(stdsize);
stdistr='stdi={';
std_incl=[];std_excl=[];
for x=1:length(stdv)
    if ~isempty(eval(['stdi' num2str(x)]))
        if stdprop(x)<.5 
            strmsg=sprintf('Not enough data points for STD%d.\nSTD NOT INCLUDED!', x);
            std_excl=[std_excl x];
            mmm=msgbox(strmsg,'Error','error','modal');
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        elseif  stdv(x)==0 && ~get(handles.rdoUse0,'Value') %exclude zero if not checked
            std_excl=[std_excl x];
        else
            stda=[stda;'STD' num2str(x)];
            stdistr=[stdistr ' stdi' num2str(x)];
            std_incl=[std_incl x];
        end
    end
end
stdistr=[stdistr '};'];
eval(stdistr);

%Groups of standards. g1 = indices of 1st elements of groups          g2= indices of last elements of groups
g11=[true;stdig(2:end)~=stdig(1:end-1)+1]; % true adds the first one
g1=stdig(g11);
g22=[stdig(1:end-1)~=stdig(2:end)-1;true]; % true adds the last one
g2=stdig(g22);

%Calculate number of excluded STD in each group: array 'g1x'
g1x(1:length(g1),1)=0;
for i=1:length(std_excl)
    for j=1:length(eval(['stdi' num2str(std_excl(i))]))
        stdixj=eval(['stdi' num2str(std_excl(i)) '(j)']);% stdixj = index of excluded stds
        xj=find(stdixj-g1<6 & stdixj-g1>=0);% xj = index of group containing excluded stds
        if (g2(xj)- stdixj)<6, g1x(xj,1)=g1x(xj,1)+1; end %that's how many excluded std are in that group (xj)
    end
end

%(g2-g1)+1=size of group (which could have an excluded standard).
%Group size = [How many STD used] - Excl STD.
%If not, 1 or more STD are missing and index of 1st std in group saved in 'miss' array
%1st row of miss will contain indices of missing STD1
%2nd row, indices of missing STD2...etc.
[row,col]=find((5 - g1x)<size(stdi,2)); %#ok<USENS> %size(stdi,2)= number of stds that will be used
for xj=row' %for each group with less than all the standards, taking the excluded stds into account
    for k=1:size(stda,1)%check which std is missing and puts index in 'miss' for each std missing
        if sum(ismember(dtype(g1(xj):g2(xj)),stda(k,:)) | ismember(dtype(g1(xj):g2(xj)),[stda(k,:) '-DRAIN']))==0
            if exist('miss','var')
                miss(k,end+1)=g1(xj);
            else
                miss(k,1)=g1(xj);
            end
        end
    end
end
if exist('miss','var')
    if size(miss,1)<size(stdi,2), miss(end+1:size(stdi,2),:)=0; end
end

clear clear -regexp stdi\d+;
clear dtype;
%std1val=NaN;std2val=NaN;std3val=NaN;std4val=NaN;
% set(handles.std1v,'String',strtrim(get(handles.std1v,'String')));
% set(handles.std2v,'String',strtrim(get(handles.std2v,'String')));
% set(handles.std3v,'String',strtrim(get(handles.std3v,'String')));
% set(handles.std4v,'String',strtrim(get(handles.std4v,'String')));

clear dstdx;
clear headp hfile datap sysinip sysinif;

%% Std Interpolation 

%%%%%% interpolate standards over measurement series - LINEARLY BETWEEN EACH VALUES %%%%%%%%%%
fstdt=[];bstdt=[];oorstdt=[];std_num=[1:length(stdv)];
%       Checks for missing std value
%find longest arrays of std (istdmax) 

stdmax=max(stdsize);istdmax=std_num(stdsize==stdmax);
clear dstd;
dstd(1:DSize,1:size(stdi,2))=-999;dstdf(1:DSize,1:size(stdi,2))=0;
for i=1:size(stdi,2)
    if ~isempty(stdi{i})
        x=dyday(stdi{i});
        Y=dlicorx(stdi{i});
        nx=dyday(si+1);
        y=interp1(x, Y,nx,'linear','extrap');

        y(find(nx<x(1)))=Y(1); %sets first values to series start
        y(find(nx>x(end)))=Y(end); %sets last values to series end
        dstd(si+1,i)=y;

        dstdf(si+1,i)=2;
        %calculated Std values are flagged 4 when before and after a missing standard or a BAD measured standard
        %Calc. Std values flagged 4 will not be used to correct xCO2 values.But flagged 3 (questionable measured standard) will
        %be used as if good.
        for j=1:length(stdi{i})
            if dflag(stdi{i}(j))>3
                if j>1 & j<length(stdi{i})
                    dstdf(find(dstdf(stdi{i}(j-1)+1:stdi{i}(j+1)-1,i)==2)+stdi{i}(j-1),i)=4;
                elseif j==1
                    dstdf(find(dstdf(1:stdi{i}(j+1)-1,i)==2),i)=4;
                else
                    dstdf(find(dstdf(stdi{i}(j-1)+1:end,i)==2)+stdi{i}(j-1),i)=4;
                end
            end
        end
        if exist('miss','var')
            if size(miss,1)>=i
                for k=miss(i,:)
                    if k~=0
                        %if missing std is first or last, don't do anything (because I don't know what to do with it)
                        if k> min(stdi{i}) & k< max(stdi{i})
                            ibefore=stdi{i}(find(stdi{i}-k<0,1,'last'));
                            iafter=stdi{i}(find(stdi{i}-k>0,1,'first'));
                            dstdf(find(dstdf(ibefore+1:iafter-1,i)==2)+ibefore,i)=4;
                        end
                    end
                end
            end
        end
        bstd(i)={find(dstdf(:,i)==4)};
    end
end
%bstdt=unique(bstdt);
clear  stda ;
clear x Y nx y xx yy;
%%%%%%%%% Plots results of Standard Interpolation for each standard on a
%%%%%%%%% separate graph

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/4 scrsz(4)/6 scrsz(3)/2 scrsz(4)*2/3],'Name','Standard Interpolation Plots','NumberTitle','off')
drawnow;
CenterWindow(handles,gcf,handles.figure1);
for i=1:size(stdi,2)
    if ~isempty(stdi{i})
        sih(i)=subplot(size(stdi,2),1,i);
        plot(dyday(stdi{i}),dlicorx(stdi{i}),'ob','MarkerSize',dotsize+1);
        hold on;
        plot(dyday(si+1),dstd(si+1,i),'og','MarkerSize',dotsize-1);%Good Std : Green
        %             plot(dyday(fstd{i}),dstd(fstd{i},i),'or','MarkerSize',dotsize-1);% Flag 3 Std : Orange
        plot(dyday(bstd{i}),dstd(bstd{i},i),'or','MarkerSize',dotsize-1);%Bad Std : Red
        ylabel(['Std ',int2str(std_incl(i))]);
        %             legend(sprintf('Std %d',i),'Interpolation','Location','Best');
    end
end
linkaxes(sih,'x');
h = axes('Position',[0 0 1 1],'Visible','off');
set(gcf,'CurrentAxes',h)

text(.25,.97,['{\fontsize{14}\fontname{times new roman}'...
    '\color{blue}o Std Values   ' '\color{green}-- Interpolated Values   ' '\color{red}-- Bad (not used)}']);
uistack(h,'bottom');
drawnow;

    %Interpolated Standards on each side of a BAD or MISSING STD are replaced by NaN and will not be used in the 
%subsequent fits...
for i=1:size(stdi,2)
    if ~isempty(stdi{i})
        dstd(bstd{i},i)=NaN;
    end
end

clear dyday bstd;

% Construct arrays of interpolated standards at each datum for both ATM5 and EQU (xx and yy)
% h = waitbar(0,'\fontname{technical}Task ^1/_3: Calculating Standards at each datum...');total=2*length(si)+DSize;
%xx{1}=[];yy{1}=[];
stdval=[];
for i=1:size(stdi,2)
      stdval(i)=stdv(std_incl(i));
end



clear  stdi stdva;
%% Standard fits at each datum for air and water


%for each sample, find position and interpolate std-USE LINEAR REGRESSION
nostd=[];
clear dlicorxcorr;
for i=1:length(si)
    index=si(i)+1;
    if length(find(~isnan(dstd(index,:))))>=2
        p=polyfit(dstd(index,find(~isnan(dstd(index,:)))),stdval(find(~isnan(dstd(index,:)))),1); %#ok<FNDSB,FNDSB>
        if length(find(~isnan(dstd(index,:))))==2 & dflag(index)~=4
            if exist('bstdt','var'),   bstdt(end+1,1)=index; else bstdt=index;  end
        end
        dlicorxcorr(index,:)=polyval(p,dlicorx(index));
        %if measured xco2 (not corrected) is out of range (more than 100 ppm) from available stds at that point, then flagged
        if (dlicorx(index)<min(dstd(index,:)-100) | dlicorx(index)>max(dstd(index,:))+100) & dflag(index)~=4
            if exist('oorstdt','var'),   oorstdt(end+1,1)=index; else oorstdt=index;  end
            %sprintf('point#: %d - %1.2f',index,dlicorx(index));
        end

    else
        nostd=cat(1,nostd,index-1);
        dlicorxcorr(index,:)=-999;
    end
end

dlicorxcorr(length(dlicorxcorr)+1:DSize)=-999;
dlicorxcorr((dlicorx(si+1)==-999),:)=-999;
dlicorxcorr(dlicorxcorr==0)=-999;

xco2corrok=1;
% if rangeok==1,    set(handles.despiket,'String', 'Check Range - incl. xCO2');  end
clear xx yy;

for i=1:length(std_incl),  DataA(2:DSize,vstdi(std_incl(i)))=cellstr(num2str(dstd(2:DSize,i),'%1.3f'));end
DataA(2:DSize,vlicorxcorr)=cellstr(num2str(dlicorxcorr(2:DSize),'%1.3f'));
DataA(strcmp(DataA,'-999.000')==1)={'-999'};

DataA(fxco2w,vlicorxw)=DataA(fxco2w,vlicorxcorr);
DataA(fxco2a,vlicorxa)=DataA(fxco2a,vlicorxcorr);

%% Flags data corrected with missing standard

Set_subflag(bstdt,9,handles); %Flags points fitted with less than 3 stds
Set_subflag(oorstdt,1,handles); %Flags points with out of range xco2
Set_subflag(nostd,0,handles); %flags data corrected with less than 2 std to 4
clear bstdt oorstdt nostd;

set(handles.xco2aok,'visible','on');

activity(0,handles); 

% close (h);




function  Check_Outliers(h, eventdata, handles, varargin)
% Sub for Callback of the uicontrol handles.Check_Outliers.
global DataA hdllohistr hdltype;
global vdt vyday Nonecol vtiniok dateok vdutc vtutc
global vteq vpeq vwflo  vgflo  vpamb  vtype  vsal  vtin vtini ...
        vflag vpatm vpeqa vtcond vlicorcav ;
global out_of_range outliers

%save('beforeCheck');

% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';...
%         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';...
%         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';... %         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';... %         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'interpolated using questionable SST'};  variable missing';'interpolated using questionable SST'};

if  ~dateok
    mmm=msgbox(sprintf('Cruise Day / Year Day needs to be calculated first.'),'modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

if isempty(DataA), return; end 
activity(1,handles); 


equ_atm={'EQU';'ATM';'EQU-DRAIN';'ATM-DRAIN'};
equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};


%detect -999 values not flagged 4 which would result in flagged data but
%only in hdl columns
hdl=[vgflo vwflo vsal vtini vteq vdt vtcond vpeqa vpamb vpatm];
hdlt={'gas flow';'H2O flow'; 'Salinity';'in situ Temp';'EQU Temp';'Delta T';'Cond Temp';'EQU Pressure';'Ambient Pressure';'ATM Pressure'};
hdlollim=[0 0 1 1 1 0 0 0 0 0];%if not 0 - check for outliers. Value is max diff between 2 points
hdltype=[];out_of_range={};outliers={};out_of_range(1,:)=hdlt';outliers(1,:)=hdlt';
for i=1:10 %checks which type to check for outliers - atm, equ or both
    % -or=out of range   -ol=outliers
    handles.("or" + num2str(i)).String='-'; handles.("or" + num2str(i)).TooltipString='-';
    handles.("ol" + num2str(i)).String='-'; handles.("ol" + num2str(i)).TooltipString='-';
    hdllohistr{i,1}=strtrim(handles.("lo" + num2str(i)).String);
    hdllohistr{i,2}=strtrim(handles.("hi" + num2str(i)).String);
    handles.("lo" + num2str(i)).TooltipString=handles.("lo" + num2str(i)).String;
    handles.("hi" + num2str(i)).TooltipString=handles.("hi" + num2str(i)).String;
    
    switch 10*get(handles.("atm" + num2str(i)),'Value') + get(handles.("equ" + num2str(i)),'Value') 
        case 11
          hdltype{i}=equ_atm;
        case 10
          hdltype{i}=atm;
        case 1
          hdltype{i}=equ;
        case 0
          hdltype{i}=equ_atm;
          %set(handles.("equ" + num2str(i)),'Value',1) ;
          handles.("equ" + num2str(i)).Value=true;
          %set(handles.("atm" + num2str(i)),'Value',1) ;
          handles.("atm" + num2str(i)).Value=true;
          mmm=msgbox(['No Type (EQU or ATM) selected for range check of ' hdlt{i} '.\nBOTH were checked automatically'],"modal");
          CenterWindow(handles,mmm,handles.figure1); uiwait(mmm);
    end
end

dyday_limit=5; %limit in minutes to observe anomalous changes in temp, or salinity. If time between 2 measurements is longer, Delta is not anomalous anymore.

indx=find(hdl~=Nonecol);%If Data is assigned ==>
if ~isempty(indx)
    %check for -999 data points
    if ~isempty(find(strcmp(DataA(2:end,hdl(indx)),'-999')~=1 & strcmp(DataA(2:end,vflag),'4')~=1 ,1))
        indx2=[];
        for ix = indx
            if ~isempty(find(strcmp(DataA(2:end,hdl(ix)),'-999')==1 & strcmp(DataA(2:end,vflag),'4')~=1 & ismember(DataA(2:end,vtype),hdltype{ix}),1)) indx2=[indx2;ix]; end
        end
        if ~isempty(indx2)
            junk=sprintf('%s\n',hdlt{indx2});
            mmm=questdlg(sprintf('You have data set to -999 which haven''t been flagged 4 in:\n\n%s',junk),'-999 issue','STOP','CONTINUE','CONTINUE');
            if strcmp(mmm,'STOP'), return; end
        end
    end

    indspk1=find(strcmp(DataA(2:end,vflag),'4')~=1);

    if ~isempty(find(strcmp(DataA(2:end,hdl(indx)),'-999')~=1 & strcmp(DataA(2:end,vflag),'4')~=1,1))

        for i = indx

            %Out of range
            indspk2=[];indspk2=find(ismember(DataA(2:end,vtype),hdltype{i}));
            lo=str2num(get(handles.("lo" + num2str(i)),'string'));hi=str2num(get(handles.("hi" + num2str(i)),'string'));
            indspk3=[];indspk3=find(str2num(char(DataA(2:end,hdl(i))))<lo | str2num(char(DataA(2:end,hdl(i))))>hi ...
                & str2num(char(DataA(2:end,hdl(i))))~=-999);

            indspk=intersect(indspk1,intersect(indspk2,indspk3));
            %set(handles.("or" + num2str(i)),'String',num2str(numel(indspk)));
            handles.("or" + num2str(i)).String=num2str(numel(indspk));handles.("or" + num2str(i)).TooltipString=num2str(numel(indspk));

            orpc=numel(indspk)/numel(ismember(DataA(2:end,vtype),hdltype{i}))*100;%out of range percentage            
            if orpc > 10 %more than 10% outliers
                 handles.("or" + num2str(i)).ForegroundColor='red';
            else
                 handles.("or" + num2str(i)).ForegroundColor='black';
            end
            out_of_range(2:numel(indspk)+1,i)=cellstr(strcat(strvcat(DataA{indspk+1,vyday}),'-',strvcat(DataA{indspk+1,vdutc}),'-',strvcat(DataA{indspk+1,vtutc})));
            %Outliers
            if hdlollim(i) ~= 0

                if hdl(i)==vtini && ~vtiniok
                    mmm=msgbox(sprintf(['Problem with insitu T data.'...
                        '\nMake sure an Offset is displayed and re-calc fields.'...
                        '\nSST not checked']),'modal');
                    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                else
                    lo=hdlollim(i);

                    indspk2=[];indspk2=find(ismember(DataA(2:end,vtype),hdltype{i}) & strcmp(DataA(2:end,vflag),'4')~=1 & str2num(char(DataA(2:end,hdl(i))))~=-999);
                    VARi=str2num(char(DataA(indspk2+1,hdl(i))));VARyd=str2num(char(DataA(indspk2+1,vyday)));
                    dbefore=zeros(length(VARi),1);dafter=zeros(length(VARi),1);dd=zeros(length(VARi),1);dydayb=zeros(length(VARi),1);dydaya=zeros(length(VARi),1);
                    dbefore(2:end,1)=abs(VARi(2:end)-VARi(1:end-1));dbefore(1,1)=0; %diff of each point with data before
                    dafter(1:end-1,1)=abs(VARi(1:end-1)-VARi(2:end));dafter(end,1)=0;%diff of each point with data after
                    dd(2:end-1)=abs(dafter(2:end-1,1)-dbefore(2:end-1,1));dd(1,1)=0;dd(end,1)=0;
                    dydayb(2:end,1)=VARyd(2:end)-VARyd(1:end-1);dydayb(1,1)=0;
                    dydaya(1:end-1,1)=dydayb(2:end,1);dydaya(end,1)=0;
                    %outlier has diff before and after about the same...
                    indspk3=find(dbefore>lo & dafter>lo & dd<lo & dydayb<(dyday_limit/60/24) & dydaya<(dyday_limit/60/24));indspk3=indspk2(indspk3);

                    indspk=intersect(intersect(indspk1,indspk2),indspk3);
                    %set(handles.("ol" + num2str(i)),'string')=num2str(numel(indspk));
                    handles.("ol" + num2str(i)).String=num2str(numel(indspk));handles.("ol" + num2str(i)).TooltipString=num2str(numel(indspk));
                    olpc=numel(indspk)/numel(ismember(DataA(2:end,vtype),hdltype{i}))*100;%out of range percentage
                    if olpc > 10 %more than 10% outliers
                        handles.("ol" + num2str(i)).ForegroundColor='red';
                    else
                        handles.("ol" + num2str(i)).ForegroundColor='black';
                    end
                    outliers(2:numel(indspk)+1,i)=cellstr(strcat(strvcat(DataA{indspk+1,vyday}),'-',strvcat(DataA{indspk+1,vdutc}),'-',strvcat(DataA{indspk+1,vtutc})));
                end

            end
        end


    end
end



% --------------------------------------------------------------------
function Set_subflag(indspk,subflag,handles) %#ok<INUSD>
global DataA vsubf vsubfu sub_flag;
global  vflag;
%sets subflag and subflagu for subflag<=10
%for values >10 subflag should be set to 'other(see metadata)' and subflagu  'subflag user' to whatever
% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';...
%         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';...
%         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';... %         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';... %         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'interpolated using questionable SST'};  variable missing';'interpolated using questionable SST'};

if subflag>10, mmm=msgbox(sprintf('Program error in ''Set_subflag'' routine.\n\nContact Denis.Pierrot@noaa.gov.'),'Flagging Error');CenterWindow(handles,mmm,handles.figure1);uiwait(mmm);return;end


if length(indspk)>0
    for i=1:length(indspk)
        if subflag==0 %flagged 4
            DataA(indspk(i),vsubf)=cellstr('');
            DataA(indspk(i),vsubfu)=cellstr('');
            DataA(indspk(i),vflag)={'4'};
        else %flagged 3
            DataA(indspk(i),vflag)={'3'};
            if subflag==-1,continue;end
            %in 'anomalous DT', ()and + are interpreted as special characters so need to be replaced before being searched
            ssearch=regexprep(char(sub_flag(subflag)),'(','\\('); ssearch=regexprep(ssearch,')','\\)'); ssearch=regexprep(ssearch,'+','\\+');
            
            vsubbf=vsubf;
            for k=1:2%Once for vsubf, once for vsubfu - sets subflag + subflag user to same string
                if k==2, vsubbf=vsubfu; end
                if k==2 & subflag==10, continue;end %do not set subflag user to 'Other - See metadata'. Should be set using Set_subflaguser routine
                if iscellstr(DataA(indspk(i),vsubbf))% some subflag already present
                    ta=indspk(i);
                    for j=1:length(ta)
                        already=regexp(DataA(ta(j),vsubbf),ssearch);%this subflag already present or not
                        if isempty(already{1})%subflag not present
                            if length(char(DataA(indspk(i),vsubbf)))==0
                                DataA(indspk(i),vsubbf)=cellstr(sub_flag(subflag));
                            else
                                DataA(indspk(i),vsubbf)=cellstr(strcat(char(DataA(indspk(i),vsubbf)),';',sub_flag(subflag)));%adds this subflag to one present
                            end
                        end
                    end
                else
                    DataA(indspk(i),vsubbf)=cellstr(sub_flag(subflag));
                end
            end
            
        end
    end
end


% --------------------------------------------------------------------
function Set_subflaguser(indspk,subflag,handles) %#ok<INUSD>
global DataA  vsubfu sub_flag ;


if length(indspk)>0
    for i=1:length(indspk)
        if subflag>10
            if iscellstr(DataA(indspk(i),vsubfu))
                ta=indspk(i);
                for j=1:length(ta)
                    already=regexp(DataA(ta(j),vsubfu),char(sub_flag((subflag))));
                    if isempty(already{1})%subflag not present
                        if length(char(DataA(indspk(i),vsubfu)))==0
                            DataA(indspk(i),vsubfu)=cellstr(sub_flag(subflag));
                        else
                            DataA(indspk(i),vsubfu)=cellstr(strcat(char(DataA(indspk(i),vsubfu)),';',sub_flag(subflag)));
                        end
                    end
                end
            else
                DataA(indspk(i),vsubfu)=cellstr(sub_flag(subflag));
            end
        end
    end
end


%%

% --------------------------------------------------------------------
function Calc_fCO2_Callback(h, eventdata, handles, varargin) %#ok<DEFNU>
% Stub for Callback of the uicontrol handles.cfco2.

global DataA DSize;
global vtype vfco2w  vfco2a  vxco2icorr  vfco2i vdfco2 vlicorxcorr;
global vfco2ok  airiok;
% global dlicorxcorr dxco2icorr;
global  vartitles ;
global  vteq vsal vtin vtini vpeqa  vpatm  Nonecol;

activity(1,handles); 
%drawnow;
DSize=size(DataA,1);
h = waitbar(0,'\fontname{technical}Initializing...');
CenterWindow(handles,h,handles.figure1);
dtype=DataA(:,vtype);
total=7;
A=str2num(char(DataA(2:end,vpeqa)));A(DSize)=0;
dpeqa=circshift(A,1) ;
waitbar(1/total,h);
A=str2num(char(DataA(2:end,vsal)));A(DSize)=0;
dsal=circshift(A,1) ;
waitbar(2/total,h);
A=str2num(char(DataA(2:end,vteq)));A(DSize)=0;
dteq=circshift(A,1) ;
waitbar(3/total,h);
A=str2num(char(DataA(2:end,vtini)));A(DSize)=0;
dtini=circshift(A,1) ;
waitbar(4/total,h);
A=str2num(char(DataA(2:end,vpatm)));A(DSize)=0;
dpatm=circshift(A,1) ;
waitbar(5/total,h);
A=str2num(char(DataA(2:end,vxco2icorr)));A(DSize)=0;
dxco2icorr=circshift(A,1) ;
waitbar(6/total,h);
A=str2num(char(DataA(2:end,vlicorxcorr)));A(DSize)=0;
dlicorxcorr=circshift(A,1) ;
waitbar(7/total,h);

siw=find(strcmp(dtype(2:end),'EQU')==1 | strcmp(dtype(2:end),'EQU-DRAIN')==1);
sia=find(strcmp(dtype(2:end),'ATM')==1 | strcmp(dtype(2:end),'ATM-DRAIN')==1);

clear dtype;
%% Calculates fugacities %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Pressure Correction for height of Deck Barometer.
% from
% http://www.npl.co.uk/reference/faqs/how-do-i-calculate-pressure-height-corrections-(faq-pressure)
% ph = pi - dgh/u
% where   ph 		is the pressure at level h metres above the measuring instrument
% pi 		is the pressure at the height of the measuring instrument
% d 		is the density of the fluid in kg.m-3
% g = 9.8		is the local value of gravitational acceleration in m.s-2.  d=9.80 m.s-2 is the standard gravity (world average)
% h 		is the vertical distance (height) of the liquid surface above the level at which the pressure value is being determined.
% u (101325) is a factor which converts the height correction term from pascals to the pressure units being used
%     here, u=1 since pressure is in mbar and correction is in mbar as
%     well.
% dgh had units of kg.m-1.s-2 = Pascal = 1 N/m2 = 1 J/m3 (1 J= 1 N.m) so 1 J = 1 Pa.m3


% from http://en.wikipedia.org/wiki/Density_of_air
% d = P/(R'*T)  where R' = R/M  (M = molar mass of air = 28.97 g/mol), and P=atm5 pressure.
   % R = 8.314 J.mol-1.K-1 so d in kg.m-3
% T (air temperature) is approximated by the Equilibrator Temperature.
% Correction will have same unit as P unit in d equation.
% At 298K and 1013.25 mbar, correction=0.1161 mbar/m



%% fCO2 for SW
%fCO2 from DOE handbook
%For Water
%determines if all the data needed is assigned
hdl=[vpeqa vsal vteq vtin];
hdlt={'EQU Pressure'; 'Salinity';'EQU Temp';'in situ Temp'};
%hdlt=strtrim(char(hdlt));
[tf,indx]=ismember(Nonecol,hdl);%If Data not assigned ==>

if tf
    mmm=msgbox(sprintf('%s\n\nnot assigned to a column.\n\nfCO2(sw) not calculated.',hdlt{indx}),'Calculation Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
end
clear hdl hdlt indx;

height=str2num(get(handles.dbaroh,'string'));
if isempty(height), height=0;end
if height<0 | height>100, height=0; end
set(handles.dbaroh,'String',strtrim(num2str((height))));

% DataA(2:end,vfco2w)={'-999'};
t_k(1:DSize,:)=NaN; S(1:DSize,:)=NaN;pH2O(1:DSize,:)=NaN; 
pCO2_sw(1:DSize,:)=NaN; fCO2_equ(1:DSize,:)=NaN; dfco2w(1:DSize,:)=NaN;

if ~tf & ~isempty(siw)
    gsiw=find(dpeqa(siw+1)>0 & dlicorxcorr(siw+1)>0 & dsal(siw+1)>0 & dteq(siw+1)>-20 & dtini(siw+1)>-20);
    
    t_k(siw(gsiw)+1,:)=dteq(siw(gsiw)+1)+273.15;
    S(siw(gsiw)+1,:)=dsal(siw(gsiw)+1);
 
    %pH2O from Weiss & Price (1980). in atm5
    pH2O(siw(gsiw)+1,:)= exp(24.4543-67.4509.*(100./t_k(siw(gsiw)+1))-4.8489.*log(t_k(siw(gsiw)+1)./100) -0.000544.*S(siw(gsiw)+1));

    %fCO2 from DOE handbook
%     pCO2_insitu temp correction from Takahashi(93)
    pCO2_sw(siw(gsiw)+1,:)=dlicorxcorr(siw(gsiw)+1).*...
        (dpeqa(siw(gsiw)+1)./1013.25-pH2O(siw(gsiw)+1)).*...
        exp(0.0423.*(dtini(siw(gsiw)+1)-dteq(siw(gsiw)+1)));
    
%     xCO2_sw(siw(gsiw)+1,:)=pCO2_sw(siw(gsiw)+1)./((dpatm(siw(gsiw)+1)+(1.2*0.098*height))./1013.25);
    
    %fCO2=xco2*p*exp{[B+2(1-xco2)^2 delta]p/RT}  xco2(ppm) p(atm5) B&delta(m3.mol) R(82.05) = fCO2 in uatm
%new
    %     fCO2_equ(siw(gsiw)+1,:)=1e6.*pCO2_sw(siw(gsiw)+1).*1e-6.*exp((-1636.75+(12.0408.*t_k(siw(gsiw)+1))...
%         -(0.0327957.*t_k(siw(gsiw)+1).^2)+(0.0000316528.*t_k(siw(gsiw)+1).^3)+2.*(1-xCO2_sw(siw(gsiw)+1).*1e-6).^2.*...
%         (57.7-0.118.*t_k(siw(gsiw)+1)).*(((dpatm(siw(gsiw)+1)+(1.2*0.098*height))./1013.25)./(82.0575.*1013.25.*t_k(siw(gsiw)+1)))));
%old
%     fCO2_equ(siw(gsiw)+1,:)=1e6.*pCO2_sw(siw(gsiw)+1).*1e-6.*exp((-1636.75+(12.0408.*t_k(siw(gsiw)+1))...
%         -(0.0327957.*t_k(siw(gsiw)+1).^2)+(0.0000316528.*t_k(siw(gsiw)+1).^3)+2.*(1-xCO2_sw(siw(gsiw)+1)).*1e-6).^2.*...
%         (57.7-0.118.*t_k(siw(gsiw)+1)).*(((dpatm(siw(gsiw)+1)+(1.2*0.098*height))./1013.25)./(82.0575.*1013.25.*t_k(siw(gsiw)+1))));
    fCO2_equ(siw(gsiw)+1,:)=pCO2_sw(siw(gsiw)+1).*exp(((-1636.75+(12.0408.*t_k(siw(gsiw)+1))...
        -(0.0327957.*t_k(siw(gsiw)+1).^2)+(0.0000316528.*t_k(siw(gsiw)+1).^3))+2.*(1-dlicorxcorr(siw(gsiw)+1).*1e-6).^2.*...
        (57.7-0.118.*t_k(siw(gsiw)+1))).*((dpeqa(siw(gsiw)+1)./1013.25)./(82.0575.*t_k(siw(gsiw)+1))));
%    fCO2coeff(siw(gsiw)+1,:)=exp(((-1636.75+(12.0408.*t_k(siw(gsiw)+1))...
%        -(0.0327957.*t_k(siw(gsiw)+1).^2)+(0.0000316528.*t_k(siw(gsiw)+1).^3))+2.*(1-dlicorxcorr(siw(gsiw)+1).*1e-6).^2.*...
%        (57.7-0.118.*t_k(siw(gsiw)+1))).*((dpeqa(siw(gsiw)+1)./1013.25)./(82.0575.*t_k(siw(gsiw)+1))));
end
dfco2w=fCO2_equ;

dfco2w(find(isnan(dfco2w)))=-999;
% save -append pCO2_sw;
clear pCO2_sw fCO2_equ ;

%% fCO2 for Air Interpolated
fCO2_air_i(1:DSize,:)=NaN;
dfco2i(1:DSize,:)=NaN;
ddfco2(1:DSize,:)=NaN;

%determines if all the data needed is assigned
hdl=[vpeqa vsal vteq vtin vpatm];
hdlt={'EQU Pressure'; 'Salinity';'EQU Temp';'in situ Temp';'ATM Pressure'};
[tf,indx]=ismember(Nonecol,hdl);%If Data not assigned ==>

if tf
    mmm=msgbox(sprintf('%s\n\nnot assigned to a column.\n\nfCO2(air) [meas. and interp.] not calculated.',hdlt{indx}),'Calculation Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
end
clear hdl hdlt indx;


%fCO2=xco2*p*exp{[B+2(1-xco2)^2 delta]p/RT}
if ~tf & airiok==1 & ~isempty(siw) & ~isempty(sia)
    gsiw=find(dpeqa(siw+1)>0 & dxco2icorr(siw+1)>0 & dsal(siw+1)>0 & dteq(siw+1)>-20 & dtini(siw+1)>-20 & dpatm(siw+1)>0);
    t_k(siw(gsiw)+1,:)=dtini(siw(gsiw)+1)+273.15;
    sea_level_corr(siw(gsiw)+1,:)=(dpatm(siw(gsiw)+1)./t_k(siw(gsiw)+1)./(8.314/28.97e-3)).* 9.8 .* height;
    %pH2O from Weiss & Price (1980). in atm5
    pH2O(siw(gsiw)+1,:)= exp(24.4543-67.4509.*(100./t_k(siw(gsiw)+1))-4.8489.*log(t_k(siw(gsiw)+1)./100) -0.000544.*S(siw(gsiw)+1));
    fCO2_air_i(siw(gsiw)+1,:)=((dpatm(siw(gsiw)+1)+sea_level_corr(siw(gsiw)+1))./1013.25-pH2O(siw(gsiw)+1)).*dxco2icorr(siw(gsiw)+1)...
        .*exp(((-1636.75+(12.0408.*t_k(siw(gsiw)+1))-(0.0327957.*t_k(siw(gsiw)+1).^2)+(0.0000316528.*t_k(siw(gsiw)+1).^3))...
        +2.*(1-dxco2icorr(siw(gsiw)+1).*1e-6).^2.*(57.7-0.118.*t_k(siw(gsiw)+1)))...
        .*(((dpatm(siw(gsiw)+1)+sea_level_corr(siw(gsiw)+1))./1013.25)./(82.0575.*t_k(siw(gsiw)+1))));

end

dfco2i=fCO2_air_i;
dfco2i(find(isnan(dfco2i)))=-999;
if ~isempty(gsiw), ddfco2(siw(gsiw)+1,:)=dfco2w(siw(gsiw)+1)-dfco2i(siw(gsiw)+1); end
ddfco2(find(isnan(ddfco2)))=-999;
ddfco2(dfco2w==-999)=-999; %index gsiw is not the same for dfco2w and dfco2i. dfco2i uses dpatm so different condition (dpatm>0) for good point
ddfco2(dfco2i==-999)=-999;

% save  dfco2w dfco2i ddfco2;
clear  fCO2_air_i;

%% fCO2 for ATM5

fCO2_air(1:DSize,:)=NaN;
dfco2a(1:DSize,:)=NaN;

if ~tf & ~isempty(sia)
    gsia=find(dpatm(sia+1)>0 & dlicorxcorr(sia+1)>0 & dsal(sia+1)>0 & dtini(sia+1)>-20 );
    t_k(sia(gsia)+1,:)=dtini(sia(gsia)+1)+273.15;
    sea_level_corr(sia(gsia)+1,:)=(dpatm(sia(gsia)+1)./t_k(sia(gsia)+1)./(8.314/28.97e-3)).* 9.8 .* height;
    S(sia(gsia)+1,:)=dsal(sia(gsia)+1);
    pH2O(sia(gsia)+1,:)= exp(24.4543-67.4509.*(100./t_k(sia(gsia)+1))-4.8489.*log(t_k(sia(gsia)+1)./100) -0.000544.*S(sia(gsia)+1));
%     fCO2_air_i(siw(gsiw)+1,:)=((dpatm(siw(gsiw)+1)+(1.2*0.098*height))./1013.25).*(1-pH2O(siw(gsiw)+1)).*dxco2icorr(siw(gsiw)+1)...
%         .*exp(((-1636.75+(12.0408.*t_k(siw(gsiw)+1))-(0.0327957.*t_k(siw(gsiw)+1).^2)+(0.0000316528.*t_k(siw(gsiw)+1).^3))...
%         +2.*(1-dxco2icorr(siw(gsiw)+1).*1e-6).^2.*(57.7-0.118.*t_k(siw(gsiw)+1)))...
%         .*(((dpatm(siw(gsiw)+1)+(1.2*0.098*height))./1013.25)./(82.0575.*t_k(siw(gsiw)+1))));
    fCO2_air(sia(gsia)+1,:)=((dpatm(sia(gsia)+1)+sea_level_corr(sia(gsia)+1))./1013.25-pH2O(sia(gsia)+1)).*dlicorxcorr(sia(gsia)+1)...
        .*exp(((-1636.75+(12.0408.*t_k(sia(gsia)+1))-(0.0327957.*t_k(sia(gsia)+1).^2)+(0.0000316528.*t_k(sia(gsia)+1).^3))...
        +2.*(1-dlicorxcorr(sia(gsia)+1).*1e-6).^2.*(57.7-0.118.*t_k(sia(gsia)+1)))...
        .*(((dpatm(sia(gsia)+1)+sea_level_corr(sia(gsia)+1))./1013.25)./(82.0575.*t_k(sia(gsia)+1))));
end

dfco2a=fCO2_air;
dfco2a(find(isnan(dfco2a)))=-999;

clear  fCO2_air;
clear gsia;
clear t_k S pH2O ;
%dfco2a already in memory
% load  dfco2w dfco2i ddfco2;
% vfco2w  vfco2a  vxco2icorr  vfco2i vdfco2;

DataA(2:end,vxco2icorr)=cellstr(num2str(dxco2icorr(2:end),'%1.2f'));
DataA(2:end,vfco2w)=cellstr(num2str(dfco2w(2:end),'%1.2f'));
DataA(2:end,vfco2i)=cellstr(num2str(dfco2i(2:end),'%1.2f'));
DataA(2:end,vdfco2)=cellstr(num2str(ddfco2(2:end),'%1.2f'));
DataA(2:end,vfco2a)=cellstr(num2str(dfco2a(2:end),'%1.2f'));
DataA(strcmp(DataA,'-999.00')==1)={'-999'};

%sets Flag 4 for -999 data in columns 'varlist' .
%varlist=[vfco2w vfco2i]; %not vfco2a because not reported in final file anyway. Will keep xCO2a in final file.
varlist=[vfco2w]; %not vfco2i either because we don't want data flagged 4 when there is no good atm measurement but equ fCO2 is good. (changed july 2-15)
gvarlist=varlist(~ismember(varlist,Nonecol));
gindex={siw siw sia};j=0;
for i=gvarlist
    j=j+1;
    Set_subflag(intersect(find(strcmp(DataA(2:end,i),'-999')==1),gindex{j})+1,0,handles);
end

clear  siw sia gsiw ;
vfco2ok=1;
activity(0,handles); 
close (h);
%%


% --------------------------------------------------------------------
function Calc_YDay(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton24.
global DataA  listf Nonecol dateok nvar;
global vdutc  vtutc vyday vydayi vcrdayi;

%Also calculates speed
M=size(DataA,1);


activity(1,handles);

if ismember(Nonecol,[vdutc vtutc])
    mmm=msgbox('No Date or Time column selected','YDay Calculation Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
   return;
end

mmm='Unknown';junk='';
smax=max(100, size(DataA,1));
if sum(ismember(['T','Z',':'],char(DataA(2:smax,vdutc))))==3 | sum(ismember(['T','Z',':'],char(DataA(2:smax,vtutc))))==3
    mmm='ISO 8601 - Zulu';
    if sum(ismember(['T','Z',':'],char(DataA(2:smax,vdutc))))~=3 
        nvar(2)=vtutc;junk=sprintf('The assignment of ''Date UTC'' has been changed to %s',listf{vtutc});
    elseif sum(ismember(['T','Z',':'],char(DataA(2:smax,vtutc))))~=3
        nvar(3)=vdutc;junk=sprintf('The assignment of ''Time UTC'' has been changed to %s',listf{vdutc});
    end
    if ~isempty(junk)
        Update_Var(handles);
        mmm2=msgbox(junk,'YDay Calculation','warn','modal');
        CenterWindow(handles,mmm2,handles.figure1) ;  uiwait(mmm2);
    end 
elseif mean(str2num(char(DataA(2:end,vdutc))))>30000
    mmm='Excel'
else
    dsplit=split(DataA(2:end,vdutc),'/');
    if size(dsplit,2)>1
        maxd1=max(str2num(char(unique(dsplit(:,1)))));
        numd1=numel(str2num(char(unique(dsplit(:,1)))));
        if maxd1>12 && maxd1<=31
            %msgbox('date seems to be in European format','modal');
            mmm='Euro.';
        else
            maxd2=max(str2num(char(unique(dsplit(:,2)))));
            numd2=numel(str2num(char(unique(dsplit(:,2)))));
            if maxd2>12 && maxd2<=31
                %msgbox('date seems to be in US format','modal');
                mmm='US';
            else
                if maxd1<=12 && maxd2<=12 && numd1>numd2, mmm='Euro.';end
                if maxd1<=12 && maxd2<=12 && numd1<numd2, mmm='US';end
            end
        end
    end
end
%time format
if strcmp(mmm,'ISO 8601 - Zulu')<0
    tsplit=split(DataA(2:end,vtutc),':');
    if size(tsplit,2)==1
        if max(str2num(char(DataA(2:end,vtutc))))<=1  % excel time format
            DataA(2:end,vtutc)=cellstr(datetime(str2num(char(DataA(2:end,vtutc))),'convertfrom','excel','format','HH:mm:ss'));
        end
    end
end

%Date Format unknown
if strcmp(mmm,'Unknown')>0
    mmm=questdlg(sprintf(['Date Format coud not be determined.\n\nHere are some dates from the file:'...
        '\n\n      %s\n      %s\n      %s\n      %s\n      %s'...
        '\n\nSelect the Date Format:\n\n      US (Month/Day/Year)\n\n      European (Day/Month/Year)\n\n      Excel (Number~40,000)\n\n'],...
        DataA{2,vdutc},DataA{max([2, floor(M/4)]),vdutc},DataA{max([2, floor(M/2)]),vdutc},DataA{max([2, floor(3*M/4)]),vdutc},DataA{end,vdutc}),...
        'Date Format','US', 'Euro.', 'Excel', 'Euro.');
    if isempty(mmm),  return;  end
end
yd(1:size(DataA,1),1)={'-999'}; yd(1,1)={'YDay Calc'};
crd(1:size(DataA,1),1)={'-999'}; crd(1,1)={'Cruise Day'};
DataA(2:end,vdutc)=regexprep(DataA(2:end,vdutc),'(^./)','0$1'); %replaces single digit days or months to 2-digit
DataA(2:end,vdutc)=regexprep(DataA(2:end,vdutc),'/(./)','/0$1'); %replaces single digit days or months to 2-digit
DataA(2:end,vdutc)=regexprep(DataA(2:end,vdutc),'/\d\d(\d\d)','/$1'); %replaces 4- digit years to 2-digit
DataA(2:end,vtutc)=regexprep(DataA(2:end,vtutc),'(^.:)','0$1'); %replaces single digit hour to 2-digit
DataA(2:end,vtutc)=regexprep(DataA(2:end,vtutc),':(.:)',':0$1'); %replaces single digit minute  to 2-digit
DataA(2:end,vtutc)=regexprep(DataA(2:end,vtutc),':(.$)',':0$1'); %replaces  single digit second  to 2-digit
dutc=char(DataA(2:end,vdutc));
tutc=char(DataA(2:end,vtutc));

if strcmp(mmm,'ISO 8601 - Zulu')>0
    badzstri=find(cellfun(@isempty,regexp(DataA(2:end,vdutc),'\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\dZ')));
    if ~isempty(badzstri)
        %find good date near 1st bad
        if badzstri(1)>2, gzstri=badzstri(1)-1;else, gzstri=badzstri+1;end
        junk=sprintf('%d Bad Date/Time found (ex: %s).\n\n Look around %s\n\n Correct and re-import data.',numel(badzstri),char(DataA(badzstri(1)+1,vdutc)),char(DataA(gzstri+1,vdutc)))
        mmm2=msgbox(junk,'YDay Calculation','warn','modal');
        CenterWindow(handles,mmm2,handles.figure1) ;  uiwait(mmm2);
    end
    tsplit=split(DataA(2:end,vtutc),'T');
    dutc=char(tsplit(:,1));
    tsplit1=split(tsplit(:,2),'Z');
    tutc=char(tsplit1(:,1));
end

try
    dutcl=cellfun(@length, DataA(2:end,vdutc));%finds dates length
    %ind=find(ismember(dutcl,size(dutc,2)));% index of dates whose length is 8. any other is bad
    % index of dates whose length is shorter than 1/1/14 or greater than 01/01/2013T01:02:03Z is bad
    ind=find(dutcl>=6 & dutcl<=20);
    %ind=find(~ismember(dutc,{'-999'}) & ~ismember(dutc,{'//'}));
    if strcmp(mmm,'Euro.')
        dutc=datestr(datenum(dutc(ind,:),'dd/mm/yy'),2);
    elseif strcmp(mmm,'US')
        dutc=datestr(datenum(dutc(ind,:),'mm/dd/yy'),2);
    elseif strcmp(mmm,'ISO 8601 - Zulu')
        dutc=datestr(datenum(dutc(ind,:),'yyyy-mm-dd'),2);
    elseif strcmp(mmm,'Excel')
        ind=find(dutcl>=4 & dutcl<=6);
        dutc1=str2num(dutc(ind,:))+ datenum('30-Dec-1899');
        dutc1=datestr(dutc1(ind,:),'mm/dd/yy');
        mmm=questdlg(sprintf(['Are these dates correct in ''month/day/year'' format?:'...
            '\n\n      %s\n      %s\n      %s\n      %s\n      %s'],...
            dutc1(1,:),dutc1(floor(M/4),:),dutc1(floor(2*M/4),:),dutc1(floor(3*M/4),:),dutc1(end,:)),...
            'Date Format','Yes', 'No',  'Cancel','Yes');
        if strcmp(mmm,'Cancel') | isempty(mmm),  return;  end
        if strcmp(mmm,'No')
            dutc2=str2num(dutc(ind,:))+ datenum('1-Jan-1904');
            dutc2=datestr(dutc2(ind,:),'mm/dd/yy');
            mmm=questdlg(sprintf(['How about these (month/day/year)?:'...
                '\n\n      %s\n      %s\n      %s\n      %s\n      %s'],...
                dutc2(1,:),dutc2(floor(M/4),:),dutc2(floor(2*M/4),:),dutc2(floor(3*M/4),:),dutc2(end,:)),...
                'Date Format','Yes', 'No',  'Cancel','Yes');
            if strcmp(mmm,'Cancel') | isempty(mmm),  return;  end
            if strcmp(mmm,'No')
                mmm=msgbox('Can''t figure out date format! YDay Calculation aborted!','Date Format Error','modal');
                dateok=0;
                CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                return;
            else %excel 1904
                dutc=dutc2;
                clear dutc2;
            end      
        else%excel 1900
            dutc=dutc1;
            clear dutc1;
         end
        
         
    end
    
    sp(1:size(dutc,1),1)=char(32);%space

    y=datevec([dutc, sp, tutc]);
    for i=1:3, t1(1:size(y,1),i)=y(1,i);end % t1 first day of cruise
    for i=4:6, t1(1:size(y,1),i)=0;end % time 00:00
    t2(1:size(y,1),1)=y(:,1);
    for i=2:3, t2(1:size(y,1),i)=1;end %t2=Jan 1st, same year
    for i=4:6, t2(1:size(y,1),i)=0;end
    %cruise day: starts at 1 at MIDNIGHT OF THE FIRST DAY of the cruise.
    crd1=etime([y(:,1) y(:,2) y(:,3) y(:,4) y(:,5) y(:,6)],t1)/86400 +1;
    crd(ind+1,1)=cellstr(num2str(crd1,'%1.5f'));
    %Year day: starts at 1 at MIDNIGHT OF THE FIRST DAY of the cruise.
    yd1=etime([y(:,1) y(:,2) y(:,3) y(:,4) y(:,5) y(:,6)],t2)/86400+1;
    yd(ind+1,1)=cellstr(num2str(yd1,'%1.5f'));

    %insert yday & cruise day IN THAT ORDER
    insert_in_data(yd(:,1),'calc',handles);
    insert_in_data(crd(:,1),'calc',handles);
    nvar(vydayi)=length(listf)-2;nvar(vcrdayi)=length(listf)-1;
    
    dflt='Year Day';
    msg=sprintf('Which x-axis Time Scale do you want to use?');
    %if data over 2 different years, yday won't be continuously increasing
    % old   if ~isempty(find(yd1(min(ind)+1:end-1,1)- yd1(min(ind)+2:end,1)> 300 ,1,'first')) %crossing end of year
    if size(unique(y(:,1)),1)>1 % more than one year in data
        dflt='Cruise Day';
        msg=sprintf('Data spans at least 2 different years\nWhich x-axis Time Scale do you want to use?');
    end
    mmm=questdlg(msg,'YDay Calculation','Year Day', 'Cruise Day', dflt);
    
    if strcmp(mmm,'Year Day')==1, nvar(1)=nvar(vydayi); elseif strcmp(mmm,'Cruise Day')==1, nvar(1)=nvar(vcrdayi) ;end
    
    if (~isempty(find(yd1(min(ind)+1:end-1,1) > yd1(min(ind)+2:end,1) ,1,'first'))) & size(unique(y(:,1)),1)==1 %dates not increasing and only 1 year in data
        badt=find(yd1(min(ind)+1:end-1,1) > yd1(min(ind)+2:end,1));
        if length(badt)>1 ||  (y(ind(badt(1))+1,1)== y(ind(badt(1))+2,1))
            mmm=msgbox(sprintf(['A few (%d) Date/Time seem to be corrupted\nSee around Year Day %.3f'...
                '\nDate: %s %s'],length(badt),str2num(char(DataA(ind(badt(1))+2,nvar(vydayi)))),...
                char(DataA(ind(badt(1))+2,vdutc)) ,char(DataA(ind(badt(1))+2,vtutc))),'modal');
            dateok=0;
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
            return;
        end
    end
        
    DataA(ind+1,vdutc)=cellstr(num2str(dutc));

catch
    disp(lasterr);
    mmm=msgbox(sprintf('Date seems to be in wrong FORMAT\n(Here is the first good date value:\n%s   %s)\nSelect proper format or correct in data file.',char(DataA(ind(1)+1,vdutc)),char(DataA(ind(1)+1,vtutc))),'Year/Cruise Day Error','error','modal');
    dateok=0;
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);activity(0,handles);
    return;
end

clear dutc tutc crd1 yd1;

dateok=1;
set(handles.ydayok,'Visible','on');drawnow;
Update_List(handles);
Update_Var(handles);
Calc_Fields_Callback(handles);%also calc speed

activity(0,handles);
%--------------------------------------------------------------------------


% --------------------------------------------------------------------
function rdoCalc_Fields(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.rdostdv.
Calc_Fields_Callback(handles);

% --------------------------------------------------------------------

% --- Executes on button press in Config_Sys_Load.
function Config_Sys_Load(hObject, eventdata, handles)
% hObject    handle to LoadPara (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global datap resultp  headp hfile vgroup vship vcruiseid sysinip sysinif where_xml def_folder expot pi_names save_opts;
global clear_savew off_on stdv dateok hdllohi hdllohistr hdltype;

equ_atm={'EQU';'ATM';'EQU-DRAIN';'ATM-DRAIN'};
equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};
hdltype={equ_atm equ equ equ equ equ_atm equ equ equ_atm equ_atm};
for i=1:10
    handles.("atm" + num2str(i)).Value=0;handles.("equ" + num2str(i)).Value=0;
    if ismember('ATM',hdltype{i}), handles.("atm" + num2str(i)).Value=true;end
    if ismember('EQU',hdltype{i}), handles.("equ" + num2str(i)).Value=true;end
    handles.("or"+num2str(i)).String='-'; handles.("or"+num2str(i)).TooltipString='-';
    handles.("ol"+num2str(i)).String='-'; handles.("ol"+num2str(i)).TooltipString='-';
end

if dateok==1
    mmm=questdlg('This will reset the data processing.','System Config Load','Reset','Cancel','Cancel');
    if strcmp(mmm,'Cancel')==1, return; end
end

if isobject(hObject), activity(1,handles); end

save_opts=[1,1,0,0,0];%default

if exist([sysinip filesep sysinif],'file')~=2,  sysinip=[def_folder , filesep , 'System Configurations'];end
if exist([sysinip filesep sysinif],'file')~=2,  sysinif='Initial System Config.ini';sysinip=def_folder;   end
if exist([sysinip filesep sysinif],'file')~=2,  if sum(ismember('WIN',computer))==3, sysinip='c:'; elseif sum(ismember('MACI',computer))==4, sysinip=''; end,    end

if isobject(hObject) | exist([sysinip filesep sysinif],'file')~=2 %#ok<*OR2>
%     ccd=pwd;cd(sysinip);
    mmm=msgbox(sprintf('%s\n\n%s\n\n','Select a SYSTEM CONFIGURATION file','in the next popup window!'),'System Config','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    [fname,fpath] = uigetfile('*.ini','Select System Config File',[sysinip filesep sysinif]);
    if (fpath==0),return;end
    sysinif=fname;
    if strcmp(fpath(end),filesep)>0,fpath=fpath(1:end-1);end
    sysinip=fpath;
end

%code below also used when loading .mat file (has to be modified there too)
if exist([sysinip filesep sysinif],'file')==2
    fid=fopen([sysinip filesep sysinif]);

    while ~feof(fid)
        lines=deblank(fgetl(fid));
        switch lines
            case 'Expocode'
                expot=deblank(fgetl(fid));set(handles.expot,'string',expot,'TooltipString',expot);%#ok<NASGU> %Expocode
            case 'Group'
                vgroup=deblank(fgetl(fid));%#ok<NASGU> %Group name
            case 'Ship'
                vship=deblank(fgetl(fid));%#ok<NASGU> %Ship name
            case 'Cruise ID'
                vcruiseid=deblank(fgetl(fid));%#ok<NASGU> %Cruise ID
            case 'Barometer Height (m)'
                lines=fgetl(fid);dbaroh=sscanf(lines,'%f');  dbaroh=num2str(dbaroh); set(handles.dbaroh,'string',dbaroh,'TooltipString',dbaroh);
            case 'EQU Pressure Differential(1-Yes 0-No)'
                lines=fgetl(fid);  set(handles.rdoequpdiff,'value',sscanf(lines,'%f'));
            case {'Standard Values (1,2,3,..)','Standard Values (1,2,3,4)'}
                lines=fgetl(fid);%#ok<NASGU> %Std values
                count=0;
                if lines~=-1
                    [stdv,count]=sscanf(lines,'%f');
                    set(handles.tblstdv,'Data',stdv,'ColumnFormat',({'bank'}));
                    tts=''; for i=1:length(stdv), tts= [tts 'STD ' num2str(i) sprintf(' : %0.2f ppm',stdv(i))];if i<length(stdv), tts=[tts sprintf('\n')]; end, end
                    set(handles.tblstdv,'TooltipString',tts);
                    stdl={};
                    for i=1:length(stdv)  stdl(i,:)={true,['STD ' num2str(i,0)]};end
                    tts=''; for i=1:length(stdv), tts= [tts 'STD ' num2str(i)];if i<length(stdv), tts=[tts sprintf('\n')]; end, end
                    set(handles.tblstdl,'Data',stdl,'TooltipString',tts);
                end
            case 'Use these Values(1-Yes 0-No)'
                lines=fgetl(fid);    set(handles.rdostdv,'value',sscanf(lines,'%f'));   
            case 'Use Zero Std in fit(1-Yes 0-No)'
                lines=fgetl(fid);    handles.rdoUse0.Value=sscanf(lines,'%f');   
            case 'min-max for gas flow, water flow, sst, sss, condsr, dt' %for version <1.30
                lines=fgetl(fid);
                [junka,count]=sscanf(lines,'%f');junka=num2str(junka);
                j=1;
                for i=[1,2,4,3,7,6]
                    handles.("lo" + num2str(i)).String=junka(j,:);handles.("lo" + num2str(i)).TooltipString=junka(j,:);
                    handles.("hi" + num2str(i)).String=junka(j+1,:);handles.("hi" + num2str(i)).TooltipString=junka(j+1,:);
                    j=j+2;
                end
                % junk='ctrla=[';
                % for i=[1,2,4,3,7,6], junk=[junk 'handles.lo' num2str(i) ' handles.hi' num2str(i) ' '];end
                % junk=[junk ']'];
                % eval(junk);%sets ctrla to array of handles lo and hi for parameters limits
                % ctrla=[handles.lo(1),handles.hi(1),handles.lo(1),handles.hi(1),handles.sstmin,handles.sstmax,handles.sssmin,handles.sssmax,...
                %     handles.condsrmin,handles.condsrmax,handles.dtmin,handles.dtmax];
                % ranges(1,:)={junka(1,:),junka(2,:),'',''}; ranges(2,:)={junka(3,:),junka(4,:),'',''};  ranges(4,:)={junka(5,:),junka(6,:),'',''};
                % ranges(3,:)={junka(7,:),junka(8,:),'',''};  ranges(7,:)={junka(9,:),junka(10,:),'',''};  ranges(6,:)={junka(11,:),junka(12,:),'',''};
                %%for i=1:length(ctrla),set(ctrla(i),'string',junka(i,:),'Tooltipstring',junka(i,:));end
                % for i=1:20,set(hdllohi(i),'string',junka(i,:),'Tooltipstring',junka(i,:));end
            case 'p min-max and t min-max for equ, licor, deck box'%for version <1.30
                lines=fgetl(fid);
                [junka,count]=sscanf(lines,'%f');
                junka([7:8,11:12])=[];%delete licor and deck box t limits - not needed anymore for v1.40
                junka=num2str(junka);
                j=1;
                for i=[8,5,9,10]
                    handles.("lo" + num2str(i)).String=junka(j,:);handles.("lo" + num2str(i)).TooltipString=junka(j,:);
                    handles.("hi" + num2str(i)).String=junka(j+1,:);handles.("hi" + num2str(i)).TooltipString=junka(j+1,:);
                    j=j+2;
                end
                % junk='ctrla=[';
                % for i=[8,5,9,10], junk=[junk 'handles.lo' num2str(i) ' handles.hi' num2str(i) ' '];end
                % junk=[junk ']'];
                % eval(junk);%sets ctrla to array of handles lo and hi for parameters limits
                % ctrla=[handles.peqmin,handles.peqmax,handles.teqmin,handles.teqmax,handles.plimin,handles.plimax,handles.tlimin,handles.tlimax,...
                %     handles.pdeckmin,handles.pdeckmax,handles.tdeckmin,handles.tdeckmax,];
                % % [junka,count]=sscanf(lines,'%f');
                % % junka([7:8,11:12])=[];%delete licor and deck box t limits - not needed anymore for v1.40
                % % junka=num2str(junka);
                % ranges(8,:)={junka(1,:),junka(2,:),'',''}; ranges(5,:)={junka(3,:),junka(4,:),'',''};  ranges(9,:)={junka(5,:),junka(6,:),'',''};
                % ranges(10,:)={junka(9,:),junka(10,:),'',''};%uses Deckbox P limits for ATM P limits
                %%for i=1:length(ctrla),set(ctrla(i),'string',junka(i,:),'Tooltipstring',junka(i,:));end
            case 'Abnormal values for sst and sss'%for version <1.30 Hard Coded in v1.40+
                lines=fgetl(fid);
            case 'min-max limits for gas flow, water flow, SSS, SST, EQU T, Delta T, Cond T, EQU P, Ambient P, ATM P'%for version >1.40 
                lines=fgetl(fid);
                [junka,~]=sscanf(lines,'%f');junka=num2str(junka);
                for i=1:10
                    handles.("lo" + num2str(i)).String=junka(2*i-1,:);handles.("lo" + num2str(i)).TooltipString=junka(2*i-1,:);
                    handles.("hi" + num2str(i)).String=junka(2*i,:);handles.("hi" + num2str(i)).TooltipString=junka(2*i,:);
                end   
                % for i=1:size(junka,1)/2,ranges(i,:)={junka(2*i-1,:),junka(2*i,:),'',''};end
                % for i=1:20,set(hdllohi(i),'string',junka(i,:),'Tooltipstring',junka(i,:));end
            case 'Range Check Flags' %not needed anymore for v1.40
                lines=fgetl(fid);%Range Check Flags
            case 'Insitu Temperature Offset(minutes)'
                lines=fgetl(fid);
                if lines~=-1, set(handles.AOvalue,'String',strtrim(lines),'Tooltipstring',strtrim(lines)); end
            case 'Correct xCO2 for water(1-Yes 0-No)'
                lines=fgetl(fid);
                if lines~=-1
                    set(handles.rdodryxco2,'value',sscanf(lines,'%f'));
                    set([handles.xco2d_f,handles.xco2d_ft],'Visible',off_on(get(handles.rdodryxco2,'value')+1,:));
                end
            case 'PI Names for this data'
                lines=fgetl(fid);
                if lines~=-1, pi_names=lines; end
            case 'Save Options (1-Yes 0-No)(Expocode as col, Group as col, Ship as col, Save flags 4, Exclude ATMs)'
                lines=fgetl(fid);
                if lines~=-1, [save_opts,count]=sscanf(lines,'%d');  end
                if length(save_opts)<5, while length(save_opts)<5 , save_opts=[save_opts ;0]; end, end
                clear_savew=1;
            otherwise
                mess=sprintf(['This System Configuration File does not have the expected format.'...
                    '\n\n Unexpected String (%s) found.'...
                    '\n\n The results are not guaranteed.'...
                    '\n\n Check all system parameters.'],lines);
                mmm=msgbox(mess);
                CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);

        end
    end
    fclose(fid);
    % tts=''; for i=1:size(ranges,1), tts= [tts char(ranges(i,1))  char(9) char(ranges(i,2)) char(9)  char(ranges(i,3)) char(9)  char(ranges(i,4)) sprintf('\n')]; end
    % set(handles.frmrange,Tooltipstring',tts);
   
    if ~exist('expot','var')
        mess=sprintf(['This Configuration File was created with an earlier version of the program.'...
            '\n\n The Ship''s Expocode will be set to ''XXXX''.'...
            '\n\n Enter the correct Expocode and save the Configuration File.']);
        mmm=msgbox(mess);
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        expot='XXXX';set(handles.expot,'string',expot,'TooltipString',expot);%#ok<NASGU> %Expocode
    end
    if ~exist('vgroup','var') | ~exist('vship','var') | ~exist('vcruiseid','var') | ~exist('dbaroh','var') | ...
            ~exist('stdv','var') |  ~exist('off_on','var') | ~exist('pi_names','var') | ...
            ~exist('save_opts','var') | ~exist('clear_savew','var') 
        mess=sprintf(['This Configuration File was created with an earlier version of the program.'...
            '\n\n The results are not guaranteed.'...
            '\n\n Check all system parameters.']);
        mmm=msgbox(mess);
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    end
    
    set(handles.paraframe,'title',sysinif);            


end

if exist([sysinip filesep sysinif],'file')==2
    Save_Ini(datap,resultp,headp,hfile,sysinip,sysinif,where_xml);
end

if isobject(hObject), Reset_Reduction(-2, 0, handles);end
if isobject(hObject),activity(0,handles); end

clear hdllohi hdla hdle hdltype;

%------------------------------------------------------------------------
% function Save_Ini(datap, resultp, headp, hfile, sysinip, sysinif, where_xml)
function Save_Ini(varargin)
global  def_folder where_ini;

hd=strvcat('Default Data File Path','Default Result File Path','Default Header File Path','Default Header File',...
           'Default System Configuration File Path','Default System Configuration File','Default xml Data Path'); %#ok<VCAT>

fid=fopen(where_ini,'r');
line='';
for i=1:14
        line=strvcat(line,fgetl(fid)); %#ok<VCAT>
end
fclose(fid);
fid=fopen(where_ini,'w');
if fid==-1
    [stat,struc] =fileattrib(where_ini);
    if struc.UserWrite~=1
        mess=fprintf('You do not have writing privileges in the following folder:\n%s\n%s.ini file not saved.',def_folder,mfilename);
        mmm=msgbox(mess,'Saving ini file');
        return;
    end
end
for i=[1 2]% Prints Header, then path under it
    fprintf(fid,'%s\n',strtrim(hd(i,:)));
    if exist(char(varargin(1,i)),'dir')~=7% if new path is not a path, prints old one
        fprintf(fid,'%s\n',strtrim(line(2*i,:)));
    else
        fprintf(fid,'%s\n',strtrim(char(varargin(1,i))));
    end
end
for i=[3 5]%Prints header for path, then path under it, then header for file, then file name 
    fprintf(fid,'%s\n',strtrim(hd(i,:)));
    % if new 'file' does not exists, prints old file
    if exist([char(varargin(1,i)),filesep,char(varargin(1,i+1))],'file')~=2
        fprintf(fid,'%s\n',strtrim(line(2*i,:)));
        fprintf(fid,'%s\n',strtrim(hd(i+1,:)));
        fprintf(fid,'%s\n',strtrim(line(2*(i+1),:)));
    else
        fprintf(fid,'%s\n',strtrim(char(varargin(1,i))));
        fprintf(fid,'%s\n',strtrim(hd(i+1,:)));
        fprintf(fid,'%s\n',strtrim(char(varargin(1,i+1))));
    end
end
for i=[7]% Prints Header, then path under it
    fprintf(fid,'%s\n',strtrim(hd(i,:)));
    if exist(char(varargin(1,i)),'dir')~=7% if new path is not a path, prints old one
        fprintf(fid,'%s\n',strtrim(line(2*i,:)));
    else
        fprintf(fid,'%s\n',strtrim(char(varargin(1,i))));
    end
end

fclose(fid);


% --- Executes on button press in SavePara.
function Config_Sys_Save(hObject, eventdata, handles)

global datap  resultp headp hfile vgroup vship vcruiseid sysinip sysinif where_xml expot pi_names save_opts;
global stdv usestdv use0 GSID baroh equpdiff dryxco2;
% hObject    handle to Config_Sys_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if isobject(hObject), button was clicked so activity updated and user asked LOCATION to save file
%if hv == -1, activity not updated but user asked LOCATION to save file
%if hv == -2, activity not updated and user NOT asked LOCATION to save file

if isobject(hObject),hv=-999; else hv=hObject; end
if isobject(hObject), activity(1,handles);end 
%drawnow;

if hv ~= -2 | exist(sysinip,'file')~=7
%     ccd=pwd;cd(sysinip);
    mmm=msgbox(sprintf('%s\n\n%s\n\n','Select a location to save the SYSTEM CONFIGURATION file','in the next popup window!'),'System Config','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    [fname,fpath] = uiputfile('*.ini','Save System Config File',[sysinip filesep sysinif]);
    if (fpath==0),return;end
    sysinif=fname;sysinip=fpath;
    %addpath (sysinip);%cd(ccd);
end

junk='hdllohi=[';
for i=1:10, junk=[junk 'handles.lo' num2str(i) ' handles.hi' num2str(i) ' '];end
junk=[junk '];'];
eval(junk);

stdv=get(handles.tblstdv,'Data');
tts=''; for i=1:length(stdv), tts= [tts 'STD ' num2str(i) sprintf(' : %0.2f ppm',stdv(i))];if i<length(stdv), tts=[tts sprintf('\n')]; end, end
set(handles.tblstdv,'TooltipString',tts);
equpdiff=get(handles.rdoequpdiff,'value'); 
usestdv=get(handles.rdostdv,'value'); 
use0=handles.rdoUse0.Value;
dryxco2=get(handles.rdodryxco2,'value'); 
baroh=get(handles.dbaroh,'string');
if isempty(str2num(baroh)), baroh=0;end
if str2num(baroh)<0 | str2num(baroh)>100, baroh='0'; end
set(handles.dbaroh,'String',strtrim((baroh)),'TooltipString',strtrim((baroh)));
expot=get(handles.expot,'string');
set(handles.expot,'String',strtrim((expot)),'TooltipString',strtrim((expot)));

GSID=strvcat(vgroup,vship,vcruiseid); %#ok<VCAT>

if (exist(sysinip,'dir')==7)
    fid=fopen([sysinip filesep sysinif],'w');
    if fid==-1
        [stat,struc] =fileattrib([sysinip filesep sysinif]);
        if struc.UserWrite~=1
            mess=fprintf('You do not have writing privileges in the following folder:\n%s\n%s file not saved.',sysinip,sysinif);
            mmm=msgbox(mess,'Saving System Config file');
            return;
        end
    end
    fprintf(fid,'%s\n','Expocode');
    fprintf(fid,'%s\n',expot);
    fprintf(fid,'%s\n','Group');
    fprintf(fid,'%s\n',vgroup);
    fprintf(fid,'%s\n','Ship');
    fprintf(fid,'%s\n',vship);
    fprintf(fid,'%s\n','Cruise ID');
    fprintf(fid,'%s\n',vcruiseid);
    fprintf(fid,'%s\n','Barometer Height (m)');
    fprintf(fid,'%s\n',baroh);
    fprintf(fid,'%s\n','EQU Pressure Differential(1-Yes 0-No)');
    fprintf(fid,'%d\n',equpdiff);
%     stdv=[sprintf('%8s',get(handles.std1v,'String')); sprintf('%8s',get(handles.std2v,'String')); sprintf('%8s',get(handles.std3v,'String')); sprintf('%8s',get(handles.std4v,'String'))];
%     fprintf(fid,'%s\n','Standard Values (1,2,3,4)');
%     fprintf(fid,'%s\t%s\t%s\t%s\n',stdv(1,:),stdv(2,:),stdv(3,:),stdv(4,:));
    stdv=get(handles.tblstdv,'Data');
    fprintf(fid,'%s\n','Standard Values (1,2,3,..)');
    for i=1:length(stdv), fprintf(fid,'%0.2f',stdv(i));if i<length(stdv), fprintf(fid,'\t'); else fprintf(fid,'\n'); end, end
    fprintf(fid,'%s\n','Use these Values(1-Yes 0-No)');
    fprintf(fid,'%d\n',usestdv);
    fprintf(fid,'%s\n','Use Zero Std in fit(1-Yes 0-No)');
    fprintf(fid,'%d\n',use0);
    fprintf(fid,'%s\n','min-max limits for gas flow, water flow, SSS, SST, EQU T, Delta T, Cond T, EQU P, Ambient P, ATM P');%for version >1.40 
    % ranges=get(handles.tblrange,'Data');
    % for i=1:size(ranges,1)-1,fprintf(fid,' %s\t %s\t',strtrim(char(ranges(i,1))),strtrim(char(ranges(i,2))));end
    % fprintf(fid,' %s\t %s\n',strtrim(char(ranges(size(ranges,1),1))),strtrim(char(ranges(size(ranges,1),2))));
    for i=1:19, fprintf(fid,'%s\t',strtrim(get(hdllohi(i),'string')));end
    fprintf(fid,'%s\n',strtrim(get(hdllohi(20),'string')));
    fprintf(fid,'%s\n','Insitu Temperature Offset(minutes)');
    fprintf(fid,'%s',strtrim(get(handles.AOvalue,'string')));set(handles.AOvalue,'TooltipString',strtrim(get(handles.AOvalue,'string')));
    fprintf(fid,'\n%s\n','Correct xCO2 for water(1-Yes 0-No)');
    fprintf(fid,'%d\n',dryxco2);    
    fprintf(fid,'%s\n','PI Names for this data');
    fprintf(fid,'%s\n',pi_names);    
    fprintf(fid,'%s\n','Save Options (1-Yes 0-No)(Expocode as col, Group as col, Ship as col, Save flags 4, Exclude ATMs)');
    for i=1:length(save_opts)-1,fprintf(fid,'%d\t',(save_opts(i)));end
    fprintf(fid,'%d\n',(save_opts(end)));

    fclose(fid);

    Save_Ini(datap,resultp,headp,hfile,sysinip,sysinif,where_xml);
    set(handles.paraframe,'title',sysinif);

end

if isobject(hObject),activity(0,handles);end 




% --------------------------------------------------------------------
function assignvar(h,varargin)

global DataA nvar listf titles vartitles;

fh=guidata(h);
nvar(get(fh.lb,'value'))=find(strcmp(strtrim(listf),strtrim(get(h,'label')))==1);
Update_Var(fh);

for i= 1:size(vartitles,1)-3
        ntitles{i}=strtrim(char(DataA(1,nvar(i))));
end
titles=strtrim(char(ntitles{1,:}));

set(fh.SaveConfig,'Visible','on');


%-------------------------------------------------------------------------
function yd = dayofyear(varargin)
global crd;
%DAYOFYEAR Ordinal number of day in a year.
%
%   DAYOFYEAR(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the ordinal
%   day number in the given year plus a fractional part depending on the
%   time of day.
%   modified from:
%   Author:      Peter J. Acklam      E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam


[year, month, day, hour, minute, second] = deal(varargin{:});

days_in_prev_months = [0; 31; 59; 90; 120; 151; 181; 212; 243; 273; 304; 334];

% Day in given month.
yd = days_in_prev_months(month) ...               % days in prev. months
    + ( isleapyear(year) & ( month > 2 ) ) ...   % leap day
    + day ...                                    % day in month
    + ( second + 60*minute + 3600*hour )/86400;  % part of day

yd=cellstr(num2str(yd,'%1.5f'));
 t1(1:length(year),1)=year(1);t1(1:length(year),2)=month(1);t1(1:length(year),3)=day(1);
% t1(1:length(year),4)=hour(1);t1(1:length(year),5)=minute(1);t1(1:length(year),6)=second(1);
t1(1:length(year),4)=0;t1(1:length(year),5)=0;t1(1:length(year),6)=0;
%cruise day: starts at 0 at MIDNIGHT OF THE FIRST DAY of the cruise.
crd=etime([year month day hour minute second],t1)/86400;
 t2(1:length(year),1)=year;t2(1:length(year),2)=(1);t2(1:length(year),3)=(1);
t2(1:length(year),4)=0;t2(1:length(year),5)=0;t2(1:length(year),6)=0;
yd2=etime([year month day hour minute second],t2)/86400;
%-------------------------------------------------------------------------
function ydf = dayofyear_fractional(varargin)
global lpy_interval;
%If data covers different years, yday is replaced by fractional year
%starting at 0 for the first year

%DAYOFYEAR Ordinal number of day in a year.
%
%   DAYOFYEAR(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the ordinal
%   day number in the given year plus a fractional part depending on the
%   time of day.
%   modified from:
%   Author:      Peter J. Acklam      E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam


[year, month, day, hour, minute, second] = deal(varargin{:});

days_in_prev_months = [0; 31; 59; 90; 120; 151; 181; 212; 243; 273; 304; 334];
lpy=isleapyear(year);
% Day in given month.
ydf = days_in_prev_months(month) ...               % days in prev. months
    + ( lpy & ( month > 2 ) ) ...   % leap day
    + day ...                                    % day in month
    + ( second + 60*minute + 3600*hour )/86400;  % part of day
ydf=(year-min(year))+ydf./(365+lpy);

ydf=cellstr(num2str(ydf,'%1.10f'));
%determine leap year intervals i0=indices for start of leap year  i1= indices for end of leap year
i0=find(lpy(1:end-1)-lpy(2:end)==-1);i1=find(lpy(1:end-1)-lpy(2:end)==1);

%lpy_interval gives indices when year is leap year.
%Used to back calculate YDay from fractional year day on graph
lpy_interval=[];
if sum(lpy)~=0
    if ismember(min(union(i0,i1)),i1),i0=cat(1,1,i0);end %adds index 1 to i0
    if ismember(max(union(i0,i1)),i0),i1=cat(1,i1,length(lpy));end %adds index 1 to i0
    if sum(lpy)==length(lpy)
        lpy_interval=[ydf(1):ydf(end)];
    else
        if length(i0)~=length(i1), mmm=msgbox('leap year problem','modal');CenterWindow(handles,mmm,handles.figure1); uiwait(mmm);end  
        for i=1:length(i0)
            lpy_interval=union(lpy_interval,[i0(i):i1(i)]);
        end    
    end
end    

%-------------------------------------------------------------------------
function t = isleapyear(year)
%ISLEAPYEAR True for leap years.
%
%   ISLEAPYEAR(YEAR) returns 1's for the elements of YEAR that are leap
%   years and 0's for those that are not.  If YEAR is omitted, the current
%   year is used.  Gregorian calendar is assumed.
%
%   A year is a leap year if the following returns true
%
%       ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400)
%
%   A year is not a leap year if the following returns true
%
%      rem(year, 4) | ( ~rem(year, 100) & rem(year, 400) )

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 12:51:45 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

error(nargchk(0, 1, nargin));

if nargin == 0               % If no input argument...
    clk = clock;              % ...get current date and time...
    year = clk(1);            % ...and extract year.
end

t = ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400);


% % --- Executes on button press in showlb.
% function showlb_Callback(hObject, eventdata, handles)
% % hObject    handle to showlb (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% dasst=get(handles.dasst,'String');
% if findstr('Show',dasst)>0
%     set(handles.lb,'visible','on');
%     set(handles.lbinstr,'visible','on');
%     set(handles.lbh,'visible','on');
%     set(handles.dasst,'String','Hide Data Assignment');
%     set(hObject,'BackgroundColor','g');
% else
%     set(handles.lb,'visible','off');
%     set(handles.lbinstr,'visible','off');
%     set(handles.lbh,'visible','off');
%     set(handles.dasst,'String','Show Data Assignment');
%     set(hObject,'BackgroundColor','r');
% end
% %drawnow;
% 

% --- Executes on button press in air_interpolation.
function Air_Interp(hObject, eventdata, handles)
% hObject    handle to air_interpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DataA vxaxis  vtype  vlicorxcorr  vxco2icorr Nonecol vflag xco2corrok airiok vsubf sub_flag dotsize;

%air interpolation

if  vxaxis== Nonecol | vtype== Nonecol 
    if vxaxis== Nonecol,txts='x-axis Time';end
    if vtype== Nonecol,txts='Type';end
    mmm=msgbox(sprintf('\"%s\" is not assigned.\nSelect column to assign.',txts),'Data Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
if xco2corrok~=1
    mmm=msgbox(sprintf('%s\n%s','xCO2 not corrected yet...','Cannot interpolate'),'Data Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
    
airiok=0;
DataA(2:end,vxco2icorr)={'-999'};
%interpolated values will be calculated for all EQU measurements
if verLessThan('matlab','9.1')
    fxco2w=find(~cellfun(@isempty,strfind(DataA(2:end,vtype),'EQU')) & str2num(char(DataA(2:end,vxaxis)))>0)+1;
else
    fxco2w=find(contains(DataA(2:end,vtype),'EQU') & str2num(char(DataA(2:end,vxaxis)))>0)+1;
end

%fxco2w=find((strcmp(DataA(2:end,vtype),'EQU')==1 | strcmp(DataA(2:end,vtype),'EQU-DRAIN')==1)& str2num(char(DataA(2:end,vxaxis)))>0)+1;%interpolated values will be calculated for all EQU measurements
s=regexp(DataA(:,vsubf),sub_flag(8));%Questionable Air Value
for i=1:length(s)
    f8(i)=~isempty(s{i});
end
%%%%%%%%%% not really needed below %%%%%%%%%%%%%%%%%%%%%
inda00=find(f8(:)'~=1); 
inda01=find(strcmp(DataA(:,vlicorxcorr),'-999')~=1 & strcmp(DataA(:,vxaxis),'-999')~=1);
inda0=intersect(inda00,inda01);
if verLessThan('matlab','9.1')
    inda1=find(~cellfun(@isempty,strfind(DataA(:,vtype),'ATM')) & strcmp(DataA(:,vflag),'3')==1);
    inda2=find(~cellfun(@isempty,strfind(DataA(:,vtype),'ATM')) & strcmp(DataA(:,vflag),'2')==1);
    fxco2aAF=find(~cellfun(@isempty,strfind(DataA(:,vtype),'ATM')) & strcmp(DataA(:,vxaxis),'-999')==0 & strcmp(DataA(:,vlicorxcorr),'-999')==0);% Air with all flags
else
    inda1=find(contains(DataA(:,vtype),'ATM') & strcmp(DataA(:,vflag),'3')==1);
    inda2=find(contains(DataA(:,vtype),'ATM') & strcmp(DataA(:,vflag),'2')==1);
    fxco2aAF=find(contains(DataA(:,vtype),'ATM') & strcmp(DataA(:,vxaxis),'-999')==0 & strcmp(DataA(:,vlicorxcorr),'-999')==0);% Air with all flags
end

inda3=intersect(inda0,inda1);% ATM NOT with quest. air value
fxco2a=union(inda3,inda2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(fxco2w)|isempty(fxco2a)|isempty(fxco2aAF)
    mmm=questdlg(sprintf('Either no EQU or no ATM to interpolate from or to..\n\nIf you choose to continue, interpolated ATM fugacities\n\n and fugacity differences will be set to ''-999'''),'Air Interpolation','Continue','Stop','Stop');
    %uiwait(mmm);
    if strcmp(mmm,'Continue')==1
       airiok=1;
       set(handles.cfco2,'Visible','on');
    end 
    return;
end
j=1;group=0;
for i=1:length(fxco2aAF)
    if i==1
        i0=1;group=1;i=i+1; %#ok<FXSET>
    elseif (fxco2aAF(i-1)~=fxco2aAF(i)-1) & i<=length(fxco2aAF) %1st of the series
       i0=i;group=1;
       if i<length(fxco2aAF)
           i=i+1; %#ok<FXSET>
       end
    end
    while (fxco2aAF(i-1)==fxco2aAF(i)-1) & group==1 & i<length(fxco2aAF)%while atm successive
        i=i+1; %#ok<FXSET>
    end
    if i==length(fxco2aAF)
       if ((fxco2aAF(i-1)==fxco2aAF(i)-1) | (i0==i)) & group==1 %i0==i means last data point is good and by itself
           i=i+1; %#ok<FXSET>
       end           
    end
    i=i-1; %#ok<FXSET>
    if group==1
        i1=i;
        gairi0=find(f8(fxco2aAF(i0):fxco2aAF(i1))'~=1);%all NOT subflags "Quest. Air"
        gairi00=find(strcmp(DataA(fxco2aAF(i0):fxco2aAF(i1),vlicorxcorr),'-999')~=1);
        gairi01=find(strcmp(DataA(fxco2aAF(i0):fxco2aAF(i1),vflag),'3')==1);
        gairi02=find(strcmp(DataA(fxco2aAF(i0):fxco2aAF(i1),vflag),'2')==1);
        gairi=union(gairi02,intersect(intersect(gairi0,gairi00),gairi01));
        if ~isempty(gairi)
            gairv=str2num(char(DataA(fxco2aAF(i0)+gairi-1,vlicorxcorr)));
            gdayv=str2num(char(DataA(fxco2aAF(i0)+gairi-1,vxaxis)));

            xco2aAVG(j)=mean(gairv);
            ydayAVG(j)=mean(gdayv);
            j=j+1;
        end
        group=0;
    end
end
xout=str2num(char(DataA(fxco2w,vxaxis)));

xin0=str2num(char(DataA(fxco2aAF,vxaxis)));
yin0=str2num(char(DataA(fxco2aAF,vlicorxcorr)));
%interpolates air values , not taking avg of 5 on
%each side.
xin=str2num(char(DataA(fxco2a,vxaxis)));
yin=str2num(char(DataA(fxco2a,vlicorxcorr)));

%interpolates air values , averaging xco2 and YDay
xin2=ydayAVG;
yin2=xco2aAVG;

if length(xin2)<2 | length(yin2)<2 | length(xout)==0
    mmm=msgbox (sprintf('Problem with Data.\nInterpolation uses only ''Flag 2''\nand ''Flag 3'' which  are not a ''Questionable Air Value''.\nCheck Data or Flags'),'Interpolation Error!','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end
% yout=interp1(xin,yin,xout,'linear');
% yout(isnan(yout))=-999;
% subplot(1,3,1);
% plot(xin0,yin0,'ok',xin,yin,':k',xout,yout,'Xg');

yout2=interp1(xin2,yin2,xout,'linear');
yout2(isnan(yout2))=-999;
yout2(find(xout<xin2(1)))=yin2(1);%sets first values to series start
yout2(find(xout>xin2(end)))=yin2(end);%sets last values to series end
% subplot(1,3,2);
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/4 scrsz(4)/6 scrsz(3)/2 scrsz(4)*2/3],'Name','ATM Interpolation','NumberTitle','off');drawnow;
CenterWindow(handles,gcf,handles.figure1);
plot(xin0,yin0,'ok',xin2,yin2,'-or',xout,yout2,'xb','MarkerSize',dotsize);
 legend('Air Values','Average', 'Interpolated','Location','Best');
% subplot(1,3,3);
% plot(xin0,yin0,'ok',xin,yin,':k',xin2,yin2,'-or',xout,yout,'Xg',xout,yout2,'+b');

for j=1:length(yout2)
    DataA(fxco2w(j),vxco2icorr)={num2str(yout2(j))};
end


airiok=1;

set(handles.cfco2,'Visible','on');

% --- Executes on button press in ydayok.posok.stdok.xco2aok
function Status(hObject, eventdata, handles)
global vxaxiscok posok stdok atmok;

offon={'off' 'on'};rangetext={'Check Range - Not xCO2' 'Check Range - incl. xCO2'};
switch  hObject
    case handles.ydayok
        [fstring s1 cc]=Status_string(handles,handles.ydayokt);
        set(handles.ydayokt,'String',sprintf(fstring,s1));
        set(handles.posok,'visible',offon{1+strcmp(cc,'g')});
        vxaxiscok=strcmp(cc,'g');
    case handles.posok
        [fstring s1 cc]=Status_string(handles,handles.posokt);
        set(handles.posokt,'String',sprintf(fstring,s1));
        set(handles.latlong_interpolation,'visible',offon{1+strcmp(cc,'g')});
        set(handles.stdok,'visible',offon{1+strcmp(cc,'g')});
        posok=strcmp(cc,'g');
%     case handles.rangeok
%         [fstring s1 cc]=Status_string(handles,handles.rangeokt);
%         set(handles.rangeokt,'String',sprintf(fstring,s1));
%         rangeok=strcmp(cc,'g');
%         set(handles.despiket,'String',sprintf('%s',rangetext{1+xco2corrok}));
    case handles.stdok
        [fstring s1 cc]=Status_string(handles,handles.stdokt);
        set(handles.stdokt,'String',sprintf(fstring,s1));
        set(handles.correctxco2,'visible',offon{1+strcmp(cc,'g')});
        stdok=strcmp(cc,'g');
%     case handles.xco2rangeok
%         [fstring s1 cc]=Status_string(handles,handles.xco2rangeokt);
%         set(handles.xco2rangeokt,'String',sprintf(fstring,s1));
%         cc1=get(handles.xco2aok,'BackgroundColor');
%         set(handles.air_interpolation,'visible',offon{max(1,strcmp(cc,'g')+isequal(cc1,[0 1 0]))});
    case handles.xco2aok
        [fstring s1 cc]=Status_string(handles,handles.xco2aokt);
        set(handles.xco2aokt,'String',sprintf(fstring,s1));
%         cc1=get(handles.xco2rangeok,'BackgroundColor');
        set(handles.air_interpolation,'visible',offon{1+strcmp(cc,'g')});
        atmok=strcmp(cc,'g');
end
set(hObject,'BackgroundColor',cc);


function [fstring s1 cc]=Status_string(handles,hObject)

if hObject==handles.ydayok
    okstr='OK';
else
    okstr='Flagged';
end
objs=get(hObject,'String');
s1 = objs(1:(findstr('-',objs)-1));
%s1=sscanf(objs,'%s - %*s %*s');
s2=sscanf(objs(findstr('-',objs):end),'- %s %*s');
s2ok=strcmp(s2,okstr);
if s2ok
    fstring=['%s- Not ' okstr];
    cc='r';
else
    fstring=['%s- ' okstr];
    cc='g';
end


% --- Executes on button press in lyRY.
function SetYaxes(hObject, eventdata, handles)
% hObject    handle to lyRY or  ryLY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch hObject
    case handles.lyRY
        gy=handles.axes2;
        by=handles.axes1;
        set(by,'YLim',get(gy,'YLim'));
    case handles.ryLY
        by=handles.axes2;
        gy=handles.axes1;
        set(by,'YLim',get(gy,'YLim'));
    case handles.blyry%both y axes limits set to highest
        ry=get(handles.axes2,'YLim');
        ly=get(handles.axes1,'YLim');
        ll=min(ry(1),ly(1));hl=max(ry(2),ly(2));
        set(handles.axes2,'YLim',[ll hl]);set(handles.axes1,'YLim',[ll hl]);
    case handles.blxrx%both x axes limits set to limit of active X axis
        rx=get(handles.axes2,'XLim');
        lx=get(handles.axes1,'XLim');
        ll=min(rx(1),lx(1));hl=max(rx(2),lx(2));
%         set(handles.axes2,'XLim',[ll hl]);set(handles.axes1,'XLim',[ll hl]);
        aa=get(handles.popGselect,'Value');
        nl=get(handles.("axes" + num2str(aa)),'XLim');set(handles.("axes" + num2str(3-aa)),'XLim',nl);
        clear aa;
end



% --- Executes on button press in latlong_interpolation.
function latlong_interpolation_Callback(hObject, eventdata, handles)
% hObject    handle to latlong_interpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DataA Nonecol vxaxis vlat vlong  vflag;

if  vxaxis== Nonecol | vlat== Nonecol | vlong== Nonecol
    if vxaxis== Nonecol,txts='x-axis Time';end
    if vlat== Nonecol,txts='Latitude';end
    if vlong== Nonecol,txts='Longitude';end
    mmm=msgbox(sprintf('\"%s\" is not assigned.\nShow Data Assignment \nand select column to assign.',txts),'Data Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

baddata=0;gooddata=0;

%interpolates lat
fdataok=find(strcmp(strtrim(DataA(2:end,vlat)),'-999')==0 & strcmp(strtrim(DataA(2:end,vlong)),'-999')==0 & strcmp(DataA(2:end,vflag),'4')==0)+1;
fdataoknot=find(strcmp(DataA(2:end,vflag),'4')==1)+1;
if isempty(fdataoknot)
    baddata=baddata+1;
elseif  isempty(fdataok)
    gooddata=gooddata+1;
else

fdatamnn=find((strcmp(strtrim(DataA(2:end,vlat)),'-999')==1 | strcmp(strtrim(DataA(2:end,vlong)),'-999')==1) & strcmp(DataA(2:end,vflag),'4')~=1)+1;
if ~isempty(fdatamnn)
    mmm=msgbox(sprintf('Some (%d) records had missing position data !\n\nThey will not be interpolated unless you flag them 4.\n\nPlot missing values by selecting the ''-999'' button.',length(fdatamnn)));
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
end
    
    xin=str2num(char(DataA(fdataok,vxaxis)));
%LATITUDE
    yin=str2num(char(DataA(fdataok,vlat)));
    xout=str2num(char(DataA(fdataoknot,vxaxis)));
    if length(xin)==0 | length(yin)==0 | length(xout)==0
        mmm=msgbox ('Problem with Lat. Data.\n Check Data or Flags','Interpolation Error!','error','modal');
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        return;
    end
    youtlat=interp1(xin,yin,xout,'linear');


%LONGITUDE
    yin=str2num(char(DataA(fdataok,vlong)));
    if length(xin)==0 | length(yin)==0 | length(xout)==0
        mmm=msgbox ('Problem with Long. Data.\n Check Data or Flags','Interpolation Error!','error','modal');
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        return;
    end
    youtlon=interp1(xin,yin,xout,'linear');

    if ~isempty(find(isnan(youtlon), 1)) | ~isempty(find(isnan(youtlat), 1))
        howmanylon='ALL';howmanylat='ALL';
        if length(find(isnan(youtlon)))<length(fdataoknot),howmanylon=sprintf('%d',length(find(isnan(youtlon))));end
        if length(find(isnan(youtlat)))<length(fdataoknot),howmanylat=sprintf('%d',length(find(isnan(youtlat))));end
        mmm=msgbox(sprintf(['%s latitudes and %s Longitudes \n\nWere Outside the Range of Flag 2 Data\n\nAnd Were Not Extrapolated.'...
            '\n\nPoints Should Be Bracketed by Flag 2 Data to Be Interpolated.\n\nThese Points Will Remain Flagged 4.'],howmanylat,howmanylon),'Error');
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        if strcmp(howmanylat,'ALL') | strcmp(howmanylon,'ALL'), return;end
    end
    
    xin=str2num(char(DataA(fdataok,vlong)));
    yin=str2num(char(DataA(fdataok,vlat)));

    if get(handles.rdodateline,'value')
        xin(xin<-3)=xin(xin<-3)+360; 
        youtlon(youtlon<-3)=youtlon(youtlon<-3)+360;       
    end
    
    figi=figure;
    CenterWindow(handles,figi,handles.figure1);
    plot(xin,yin,'-k',youtlon,youtlat,'or');
    legend('Good Position', 'Interpolated Pos.','Location','Best');
            
    set(gcf,'Toolbar','figure')
    axl=0.15;axb=0.2;axw=0.75;axh=0.75;
    gcfp=get(gcf,'position'); set(gca,'position',[axl axb axw axh]);
    fl=gcfp(1);fb=gcfp(2);fw=gcfp(3);fh=gcfp(4);
    butw=fw/8;buth=0.1*fh;clear gcbf;
    bt(1)=uicontrol('Style','radiobutton','Position',[(fw-butw)/2 2 butw buth],'String','Accept...',...
        'Callback','uiresume','Value',0);
    bt(2)=uicontrol('Style','radiobutton','Position',[(fw+butw)/2 2 butw buth],'String','Reject...',...
        'Callback','uiresume','Value',0);
    uiwait(figi);
    mmm=find(cell2mat(get(bt,'Value')));
    close(figi);
    %mmm=questdlg(sprintf('Accept Interpolation?'),'Interpolation Results','Yes','No','No');
    %close ;
    if mmm==1
        DataA(fdataoknot(~isnan(youtlat)),vlat)=cellstr(num2str(youtlat(~isnan(youtlat))));
        DataA(fdataoknot(~isnan(youtlon)),vlong)=cellstr(num2str(youtlon(~isnan(youtlon))));
%         for j=1:length(fdataoknot)
%             DataA(fdataoknot(j),vlat)={num2str(youtlat(j))};
%             DataA(fdataoknot(j),vlong)={num2str(youtlon(j))};
%         end
        Calc_Speed(handles);
        Set_subflag(fdataoknot(~isnan(youtlat)),10,handles);
        Set_subflaguser(fdataoknot(~isnan(youtlat)),16,handles);
        Check_for_999(handles);
    end
end
%interpolates long
% fdataok=find(strcmp(strtrim(DataA(2:end,vlong)),'-999')==0 & strcmp(strtrim(DataA(2:end,vlong)),'0')==0 & strcmp(DataA(2:end,vflag),'4')==0)+1;
% fdataoknot=find(strcmp(strtrim(DataA(2:end,vlong)),'-999')==1 | strcmp(strtrim(DataA(2:end,vlong)),'0')==1)+1;
% if isempty(fdataoknot)
%     baddata=baddata+2;
% elseif  isempty(fdataok)
%     gooddata=gooddata+2;
% else
% 
%     xin=str2num(char(DataA(fdataok,vxaxis)));
%     yin=str2num(char(DataA(fdataok,vlong)));
%     xout=str2num(char(DataA(fdataoknot,vxaxis)));
%     if length(xin)==0 | length(yin)==0 | length(xout)==0
%         msgbox ('Problem with Long. Data.\n Check Data or Flags','Interpolation Error!','error','modal');
%         uiwait;
%         return;
%     end
%     yout=interp1(xin,yin,xout,'linear');
%     yout(isnan(yout))=-999;
%     subplot (2,1,2);
%     plot(xin,yin,'-k',xout,yout,'or');
%     legend('Good Longitude', 'Interpolated Long.','Location','Best');
%     for j=1:length(yout)
%         DataA(fdataoknot(j),vlong)={num2str(yout(j))};
%     end
% end
% mstring='';
% switch baddata
%     case 1
%         mstring='Lat.';
%     case 2
%         mstring='Long.' ;
%     case 3
%         mstring='Lat. or Long.';
% end
if baddata~=0
    mmm=msgbox('Couldn''t find bad position data','interpolation Error!','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

if gooddata~=0
    mmm=msgbox('Couldn''t find good position data','interpolation Error!','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

 
% set(hObject,'BackgroundColor','g');



% --- Executes on button press in SaveConfig.
function Config_Data_Save(hObject, eventdata, handles)
% hObject    handle to SaveConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DataA vartitles datap resultp hfile headp nvar sysinip sysinif where_xml;
global vydayi vcrdayi Nonecol;

if isobject(hObject),hv=-999; else hv=hObject; end

if hv~=-2 | exist(headp,'file')~=7    
    mmm=msgbox(sprintf('%s\n\n%s\n\n','Select a location to save the DATA CONFIGURATION file','in the next popup window!'),'Data Config','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    [fname,fpath] = uiputfile('*.csv','Save Data Config File',[headp filesep hfile]);
    if (fpath==0),return;end
    headp=fpath;hfile=fname;%addpath(headp);%cd(ccd);
    Save_Ini(datap,resultp,headp,hfile,sysinip,sysinif,where_xml);
end

fid=fopen([headp filesep hfile],'w');
if fid==-1
    [stat,struc] =fileattrib([headp filesep hfile]);
    if struc.UserWrite~=1
        mess=fprintf('You do not have writing privileges in the following folder:\n%s\n%s file not saved.',headp,hfile);
        mmm=msgbox(mess,'Saving Data Config file');
        return;
    end
end

%does not save QC,SubFlag or SubFlag_User column number. Assume these are never loaded from data file
%hence the '- 3'
mess='';
for i= 1:size(vartitles,1)-3
        mess=[mess sprintf('%s, ',strtrim(vartitles(i,:)))];
end
fprintf(fid,'%s\n',mess(1:end-2));

mess='';%nox=[vydayi;vcrdayi];%x-axis cannot be assigned to YDay calc or cruise day calc
for i= 1:size(vartitles,1)-3
    mess=[mess sprintf('%s, ',strtrim(char(DataA(1,nvar(i)))))];
    ntitles{i}=strtrim(char(DataA(1,nvar(i))));
end
fprintf(fid,'%s\n',mess(1:end-2));
titles=strtrim(char(ntitles{1,:}));

fclose(fid);

set(handles.loadhct,'String',hfile);
clear ntitles;



% --- Executes on button press in Stat_Update.
function Stat_Update_Callback(hObject, eventdata, handles)
% hObject    handle to Stat_Update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DataA vflag vsubf sub_flag ;

if isempty(DataA), return; end 
set(handles.Stat_Update,'Visible','off');pause(0.5);set(handles.Stat_Update,'Visible','on');
sub_flagt=sub_flag;
%replaces special characters ()+
sub_flagt=regexprep(sub_flagt,'\(','\\\(');
sub_flagt=regexprep(sub_flagt,'\+','\\\+');
sub_flagt=regexprep(sub_flagt,'\)','\\\)');
set(handles.sfl2,'string',num2str(length(find(strcmp(DataA(:,vflag),'2')==1))),'Tooltipstring',num2str(length(find(strcmp(DataA(:,vflag),'2')==1))));
set(handles.sfl3,'string',num2str(length(find(strcmp(DataA(:,vflag),'3')==1))),'Tooltipstring',num2str(length(find(strcmp(DataA(:,vflag),'3')==1))));
set(handles.sfl4,'string',num2str(length(find(strcmp(DataA(:,vflag),'4')==1))),'Tooltipstring',num2str(length(find(strcmp(DataA(:,vflag),'4')==1))));

for k=1:10
    if k==9
        k=k+1; %#ok<FXSET> %Skips Interpolated Std
    end
    s=regexp(DataA(:,vsubf),sub_flagt(k));
    j=sum((~cellfun(@isempty,s)));
    set(handles.("sf" + num2str(k)),'string',num2str(j),'TooltipString',num2str(j));
end


% --- Despike flag buttons.
function despikefl(hObject, eventdata, handles)
nflc={'r'; 'b'};
nfl=get(hObject , 'value') + (get(hObject , 'value')==3) - (get(hObject , 'value')==4);
nfls=sprintf('Flag %d',nfl);
set(hObject,'value',nfl,'String',nfls,'TooltipString',nfls,'ForegroundColor',char(nflc(1+(get(hObject , 'value')==4))));

%% OFFSET Tin

function Offset(hObject, eventdata, handles)
global DataA vtin vtini vteq vxaxis;
global t dotsize;

gs=get(handles.popGselect,'String');
gsn=get(handles.popGselect,'Value');
if hObject == handles.OSet && handles.OManual.Value==1
    if ~ismember(DataA(1,vtin),gs) | ~ismember(DataA(1,vteq),gs)
        mmm=msgbox (sprintf('Both "InSitu Temperature" and "EQU Temperature"\nneed to be on the graph'),'Offset Determination Error!','error','modal');
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        return;
    end
end

switch hObject
    case handles.OSet
        if strcmp(handles.OSettxt.String,'Set Offset')>0
            if handles.OManual.Value==1
                handles.OSettxt.String='Stop';handles.OSet.BackgroundColor='red';
                handles.OManual.Enable='off';handles.OAuto.Enable='off';
                handles.OApply.Visible='off';handles.OCheck.Visible='off';
                %if vtin is on graph but not active, make active
                if strcmp(DataA(1,vtin),gs(gsn))==0
                    set(handles.popGselect,'Value',3-gsn);
                    lstGselect_Callback(handles.popGselect,1,handles,1);
                    mmm=msgbox (sprintf('"%s" made active on the graph',char(gs(3-gsn))),'Offset Determination Error!','correction','modal');
                    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                end
                set(handles.popGselect,'Enable','off');

                % set offset to value in box
                x2l=get(handles.axes2,'XLim');
                tOTin=str2num(get(handles.Ovalue,'String'));
                set(handles.axes1,'XLim',[(x2l(1)-tOTin/60/24) (x2l(2)-tOTin/60/24)]);
                % start panning
                pan(handles.figure1,'on');    set(handles.blxrx,'Visible','off');
                %         t = timer('TimerFcn',{@OTin_Calc,handles});
                t = timer('TimerFcn',{@OTin_Calc,handles});
                set(t,'ExecutionMode','fixedRate','BusyMode','drop','Period',0.1);
                start(t);
            else %Automatic
                set(handles.Ovalue,'String','-','TooltipString','-');
                tOTin = fminsearch(@Auto_SST_Offset,0);
                if tOTin<0
                    set(handles.Ovalue,'String',num2str(tOTin,'%.2f'),'TooltipString',num2str(tOTin,'%.2f')); 
                    mmm=msgbox (sprintf('Offset Calculated is negative ( %.2f ).\nBetter check manually.',tOTin),'Offset Auto Determination Error!','correction','modal');
                    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                elseif tOTin>5
                    set(handles.Ovalue,'String',num2str(tOTin,'%.2f'),'TooltipString',num2str(tOTin,'%.2f')); 
                    mmm=msgbox (sprintf('Offset Calculated is greater than 5 ( %.2f ).\nBetter check manually.',tOTin),'Offset Auto Determination Error!','correction','modal');
                    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                else %minimum found
                    set(handles.AOvalue,'String',num2str(tOTin,'%.2f'),'TooltipString',num2str(tOTin,'%.2f')); 
                    set(handles.Ovalue,'String',num2str(tOTin,'%.2f'),'TooltipString',num2str(tOTin,'%.2f')); 
                    Calc_Fields_Callback(handles);%offset applied in Calc_Fields routine
                    %Config_Sys_Save(-2, -1, handles);
                    set(handles.popGselect,'Enable','on');
                end
            end
        else %Stop 
            handles.OSettxt.String='Set Offset';handles.OSet.BackgroundColor=[0.502 0.502 1.000];
            handles.OManual.Enable='on';handles.OAuto.Enable='on';
            handles.OApply.Visible='on';handles.OCheck.Visible='on';
            pan(handles.figure1,'off'); set(handles.blxrx,'Visible','on');
            if isvalid(t),  stop(t);delete(t);clear t;  end
            set(handles.popGselect,'Enable','on');

        end
    case handles.OApply
        pan(handles.figure1,'off'); set(handles.blxrx,'Visible','on');
        if handles.OManual.Value==1, if isvalid(t),  stop(t);delete(t);clear t;  end, end
        gs=get(handles.popGselect,'String');
        gsn=get(handles.popGselect,'Value');
        tOTin=str2num(get(handles.Ovalue,'String'));% when Tin original is  plotted, new offset is full offset b/c determined from original data
        set(handles.AOvalue,'String',num2str(tOTin,'%10.2f'),'TooltipString',num2str(tOTin,'%10.2f')); set(handles.Ovalue,'String','0.00','TooltipString','0.00');
        Calc_Fields_Callback(handles);%offset applied in Calc_Fields routine
        Config_Sys_Save(-2, -1, handles);
        set(handles.popGselect,'Enable','on');
    case handles.OCheck
        pan(handles.figure1,'off');set(handles.blxrx,'Visible','on');
        set(handles.popGselect,'Enable','on');
        if handles.OManual.Value==1, if isvalid(t),  stop(t);delete(t);clear t;  end, end
        tOTin=str2num(get(handles.Ovalue,'String'));
        if tOTin==0
            mmm=msgbox (sprintf('Offset Calculated is 0.\nNo interpolation done.'),'Offset Check','none','modal');
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        else
            xin=str2num(char(DataA(2:end,vxaxis)));
            yin=str2num(char(DataA(2:end,vtin)));
            yeq=str2num(char(DataA(2:end,vteq)));
            xout=xin-tOTin/60/24;
            yout=interp1(xin,yin,xout,'linear','extrap');
            plot(xin,yin,'-ok',xin,yout,'or',xin,yeq,'ob','MarkerSize',dotsize);
            legend('InSitu T','Interpolated','EQU T','Location','Best');
            axes(gca);
        end 
end


function OTin_Calc(obj, event,handles)

%Displays Time Offset between the 2 x-axes
%Used to display SST time offset

x1l=get(handles.axes1,'XLim');x2l=get(handles.axes2,'XLim');
tOTin=(x2l(1)-x1l(1))*24*60;
set(handles.Ovalue,'String',num2str(tOTin,'%10.2f'),'TooltipString',num2str(tOTin,'%10.2f'));


function infer_temp_Callback(hObject, eventdata, handles)
global DataA ind1 ind2 Currentplot InferTO VarO;
global listf;
global mom dotsize monitors;
global vtype vxaxis Nonecol;
global  vlat vlong vpeq vteq vwflo  ...
    vgflo  vlicorcav  vpamb vsal  vtin vlicorx  ...
    vlicorw vpatm vtcpu vtdbox vtcond vtini;

%Either determines offset between pairs of plotted data
%Or Applies previously determined offset to selected data
%Offset can be entered manually.

activity(1,handles); 
equ_atm={'EQU';'ATM';'EQU-DRAIN';'ATM-DRAIN'};
equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};

if hObject==handles.gd1eo | hObject==handles.gd2eo% Enter Offset
    
    mom=1; %mom = manual offset multiplier...
    mmm=inputdlg('Enter Offset Value');   
    if (~isempty(mmm) & isnumeric(str2num(char(mmm))))
        InferTO=str2num(char(mmm));
        mess=sprintf('Apply Offset(~ %g) on Selected Data',InferTO);
        set(handles.gd1ao,'Label',mess);
        set(handles.gd2ao,'Label',mess);
    end
    return;
end

gs=get(handles.popGselect,'String');
%allowed=[vtin vteq vpamb vpatm];
suggested=[vlat vlong vpeq vteq vwflo  ...
    vgflo  vlicorcav  vpamb  vsal  vtin vlicorx  ...
    vlicorw vpatm vtcpu vtdbox vtcond vtini];

    
gsok(1)=strmatch(gs(1,:),listf);
gsok(2)=strmatch(gs(2,:),listf);

% if sum(ismember(gsok,allowed))~=2
%     mess=sprintf('Pair of variables on graph:\n\n\t(%s - %s)\n\n',char(gs(1,:)),char(gs(2,:)));
%     mess=[mess sprintf('Acceptable pairs of variables for this function are:\n\n     (SST - EQU t) = (%s - %s)\n\n      or\n\n      (Licor P - ATM5 P) = (%s - %s)',...
%         char(listf(allowed(1))),char(listf(allowed(2))),char(listf(allowed(3))),char(listf(allowed(4))))];
%     mmm=msgbox(mess,'Infer Inactive Data','modal');
% %    mmm=msgbox(sprintf('Acceptable pairs of variables for this function are:\n\t\tEQU t/SST  or  ATM5 P/LICOR P'),'Infer Inactive Data','modal');
%     uiwait(mmm);
%     return;
% end

if sum(ismember(gsok,suggested))~=2
    if ~ismember(gsok(1),suggested), mess1=gs(1,:);else mess1=gs(2,:);end
    mess=sprintf('The variable: "%s" on the graph\n\n%s\n\n%s\n\n%s\n\n%s\n',char(mess1),'is not on the list of variables','suggested for this function.',...
        'Are you sure you want to continue?', 'List of suggested variables:');
    for i=1:size(suggested,2), if suggested(i)~=Nonecol, mess=sprintf('%s\n   - %s',mess, char(listf(suggested(i))));end;end
    mmm=questdlg(mess,'Infer Data','Continue','Abort','Continue');
%    mmm=msgbox(sprintf('Acceptable pairs of variables for this function are:\n\t\tEQU t/SST  or  ATM5 P/LICOR P'),'Infer Inactive Data','modal');
%     uiwait(mmm);
    if strcmp(mmm,'Abort')==1, return;end
end

rh=findobj(handles.axes1,'Type','Rectangle');
axY=handles.LeftY; otherY=handles.RightY; ind=ind1; omo=1; %omo = offset multiplier...offset is (2)-(1)...when applied to (1)...(1)+1*offset=(1)+[(2)-(1)]=(2)
if isempty(rh)
    rh=findobj(handles.axes2,'Type','Rectangle');
    if ~isempty(rh)
        axY=handles.RightY; otherY=handles.LeftY; ind=ind2; omo=-1;%when offset applied to (2): (2) - [(2)-(1)]= (1)
    else
        mmm=msgbox(sprintf('No Area found on either axis...\nAn area over which to do the inference\nneeds to be selected.'),'Infer Inactive Data','modal');
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        clear axY otherY ind omo ;
        return;
    end
end



if hObject==handles.gd1do | hObject==handles.gd2do% Determine Offset
    ptagActLall=zeros(length(ind1{1}),1);ptagActRall=zeros(length(ind2{1}),1); 
    for r=1:length(rh)% for each rectangle
        ptagActL=[];ptagActR=[];
        rv=get(rh(r),'Position');
        xv=[rv(1) rv(1) rv(1)+rv(3) rv(1)+rv(3)];
        yv=[rv(2) rv(2)+rv(4) rv(2)+rv(4) rv(2)];
        maxly=max(str2num(char(DataA(ind{1},handles.LeftY.Value))))+1;
        minly=min(str2num(char(DataA(ind{1},handles.LeftY.Value))))-1;
        ymv=[minly maxly maxly minly];%will select any flag2(ind{1}) point between x values since y is min and max
        % bad is flag3(ind{2}) or flag4(ind{3}) and has x value between first and last good data point used for interpolation
        % average of good right values - good left values  b/c ind1{1} is only flags 2
        mom=0; %mom = manual offset multiplier...
        ptagActL = inpolygon(str2num(char(DataA(ind1{1},handles.CurrentX.Value))),str2num(char(DataA(ind1{1},handles.LeftY.Value))),xv,ymv);
        ptagActR = inpolygon(str2num(char(DataA(ind2{1},handles.CurrentX.Value))),str2num(char(DataA(ind2{1},handles.RightY.Value))),xv,ymv);
        ptagActLall=ptagActLall | ptagActL;  %ptagAct(L or R) are logical arrays
        ptagActRall=ptagActRall | ptagActR;
    end
    ptagActA= unique(union(ind1{1}(ptagActLall),ind2{1}(ptagActRall)));
    diffA = find(strcmp(DataA(ptagActA,handles.RightY.Value),'-999')==0 & strcmp(DataA(ptagActA,handles.LeftY.Value),'-999')==0);
    %take difference between same indices for both arrays (right and left)
    DeltaRmL=mean(str2num(char(DataA(ptagActA(diffA),handles.RightY.Value)))-str2num(char(DataA(ptagActA(diffA),handles.LeftY.Value))));
    StdevRmL=std(str2num(char(DataA(ptagActA(diffA),handles.RightY.Value)))-str2num(char(DataA(ptagActA(diffA),handles.LeftY.Value))));
    if isnan(DeltaRmL)
        mmm=msgbox(sprintf('No Offset can be calculated...\nWith the data selected.\nTry again.'),'Infer Inactive Data','modal');
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        return;
    end
    prec=int2str(2*(DeltaRmL>0.5)+3*(DeltaRmL<=0.5));
    frmt=['%1.' prec 'f'];
    mess=sprintf(['Store Offset?\n ' frmt ' (+- ' frmt ')      n = %i'],DeltaRmL,StdevRmL,numel(diffA));
    mmm=questdlg(mess,'Infer Data','Yes','No','No');
    if strcmp(mmm,'Yes')
        InferTO=DeltaRmL;
        mess=sprintf(['Apply Offset(~ ' frmt ') on Selected Data'],InferTO);
        set(handles.gd1ao,'Label',mess);
        mess=sprintf(['Apply Offset(~ ' frmt ') on Selected Data'],-InferTO);
        set(handles.gd2ao,'Label',mess);
        VarO=[handles.LeftY.Value handles.RightY.Value];%variables [(1), (2)]used to determine offset (2)-(1) = (right - left)
    end
    
else % Apply Offset
    if ~exist('InferTO','var')
        mmm=msgbox(sprintf('No Offset has been determined yet!\nTry again.'),'Infer Inactive Data','modal');
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        return;
    end
    
    %if (VarO(1)~=get(axY,'Value') | VarO(2)~= get(otherY,'Value')) %NOT same plot as when offset determined
    if ~isempty(VarO)  % if VarO empty, offset not determined from graph (entered manually)
        if sum(ismember([get(axY,'Value') get(otherY,'Value')],VarO))~=2 %NOT same variables plotted as when offset determined
            mess=sprintf('%s.\n\n%s\n\n    Active variable: %s\n\n    Inactive variable: %s)','WARNING: Variables on graph was modified since Offset was determined:',...
                'Offset determined with:',char(listf(VarO(1))),char(listf(VarO(2))));
            mmm=msgbox(mess,'Infer Inactive Data','modal');
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        end
    end
    
    type='';types='';
    switch get(otherY,'Value')
        case {vtin,vtini}
            flag=3;subflag=2;type=equ;types='EQU';
        case vteq
            flag=3;subflag=3;type=equ;types='EQU';
        case vpamb
            flag=3;subflag=6;type=equ_atm;types='EQU and ATM';
        case vpatm
            flag=3;subflag=6;type=atm;types='ATM';
        case vsal
            mess='Choose how to flag the generated data!';
            mmm=questdlg(mess,'Infer Salinity','No Change','Flag 3','No Change');
            flag=2; if strcmp(mmm,'Flag 3'), flag=3;end
            subflag=5;type=equ;types='EQU';
        otherwise
            flag=3;subflag=10;type=equ_atm;types='EQU and ATM';
    end
    ptagActall=zeros(length(DataA)-1,1);
    for r=1:length(rh)% for each rectangle
        ptagAct=[];
        rv=get(rh(r),'Position');
        xv=[rv(1) rv(1) rv(1)+rv(3) rv(1)+rv(3)];
        yv=[rv(2) rv(2)+rv(4) rv(2)+rv(4) rv(2)];
        maxly=max(str2num(char(DataA(ind{1},handles.LeftY.Value))))+1;
        minly=min(str2num(char(DataA(ind{1},handles.LeftY.Value))))-1;
        ymv=[minly maxly maxly minly];%will select any flag2(ind{1}) point between x values since y is min and max
        % bad is flag3(ind{2}) or flag4(ind{3}) and has x value between first and last good data point used for interpolation
        ptagAct = inpolygon(str2num(char(DataA(2:end,handles.CurrentX.Value))),str2num(char(DataA(2:end,get(axY,'Value')))),xv,yv);
        ptagActall=ptagActall | ptagAct;
    end
    ptagActall=cat(1,false,ptagActall);%adds a 0 to first position
    %all T data selected
    xx=str2num(char(DataA(ptagActall,handles.CurrentX.Value)));
    ys=str2num(char(DataA(ptagActall,get(axY,'Value'))));
    
    if isempty(ys)
        mmm=msgbox(sprintf('No data to apply the offset to !\n\nTry again.'),'Infer Inactive Data','modal');
        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        return;
    end
    
    if isempty(VarO), omo=1;
    elseif (VarO(1)==get(axY,'Value')),omo=1;elseif (VarO(2)==get(axY,'Value')),omo=-1;
    end
    
    %all other T data in X selected
    yns=str2num(char(DataA(ptagActall,get(otherY,'Value'))));
    %ADJUSTED other T data from FLAG2 selected
    if (mom==1), om=1; else om=omo;end
    yas=ys + (om * InferTO);
    
    figi=figure;
    CenterWindow(handles,figi,handles.figure1);
    plot(xx,ys,'d','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',dotsize+1);
    gcfp=get(gcf,'position');set(gcf,'position',[gcfp(1) gcfp(2) 0.5*monitors(1,3) 0.5*monitors(1,4)]);CenterWindow(handles,gcf,handles.figure1);

    hold on;
    plot(xx,yns,'o','MarkerEdgeColor','b','MarkerSize',dotsize+1);
    
    %Plots adj. data in green color
    c2=[0.5 .8 0.5];%dark green
    plot(xx,yas,'o','MarkerFaceColor',c2,'MarkerEdgeColor',c2,'MarkerSize',dotsize-1);
    set(gca,'Color',get(handles.figure1,'Color'));
    if flag==3
        lg={'Reference Data';'Data Replaced';sprintf('New Data (will be flagged %d)',flag)};
    else
        lg={'Reference Data';'Data Replaced';sprintf('New Data (flag not changed)')};
    end
    legend(lg,'Location','Best');
    axes(gca);
    
    set(gcf,'Toolbar','figure')
    axl=0.15;axb=0.2;axw=0.75;axh=0.75;
    gcfp=get(gcf,'position'); set(gca,'position',[axl axb axw axh]);
    fl=gcfp(1);fb=gcfp(2);fw=gcfp(3);fh=gcfp(4);
    butw=fw/3;buth=0.1*fh;clear gcbf;
    if flag==3
        bt(1)=uicontrol('Style','radiobutton','Position',[fw/4 2 butw buth],...
            'String',sprintf('Accept New Data? (Only %s will be flagged %d)',types,flag),...
            'Callback','uiresume','Value',0);
    else
        bt(1)=uicontrol('Style','radiobutton','Position',[fw/4 2 butw buth],...
            'String','Accept New Data? (not flagged)',...
            'Callback','uiresume','Value',0);
    end
    bt(2)=uicontrol('Style','radiobutton','Position',[(fw+butw)/2 2 butw buth],...
        'String','Reject New Data?', 'Callback','uiresume','Value',0);
    uiwait(figi);
    mmm=find(cell2mat(get(bt,'Value')));
    close(figi);

    
%     uicontrol('Position',[(fw-butw)/2 2 butw buth],'String','Continue...','Callback','uiresume');
%     uiwait;
    
%     if flag==3
%         mmm=questdlg(sprintf('Accept New Data?\n\nOnly %s will be flagged %d',types,flag),'Inference Results','Yes','No','No');
%     else
%         mmm=questdlg(sprintf('Accept New Data?'),'Inference Results','Yes','No','No');
%     end
%     close ;
    
    if mmm==1
        
        if get(otherY,'Value')==vtini
            %  tin is shifted and interpolated. offset needs to be applied to original Tin
            %  so that shifted interpolated data has correct offset. That's not possible (values diverge
            %  quickly). Next best thing is to set original Tin values to the new values. The interpolated values
            %  will be close (though not with the exact offset one wanted, especially when gradient is strong)
            %
            
            OTin=str2num(get(handles.AOvalue,'String'));
            ydaymin=min(str2num(char(DataA(ptagActall,vxaxis))))+OTin/60/24;
            %ydaymax=max(str2num(char(DataA(ptagAct,vxaxis))))+OTin/60/24;
            iTinmin=(str2num(char(DataA(2:end,vxaxis)))<ydaymin);
            iTinmin=cat(1,true,iTinmin);   iTin1=find(iTinmin,1,'last');
            %iTinmax=(str2num(char(DataA(2:end,vxaxis)))>ydaymax);
            %iTinmax=cat(1,false,iTinmax);   iTin2=find(iTinmax,1,'first') -1;
            clear iTinmin  ydaymin ;
            
            DataA(iTin1:iTin1+length(yas)-1,vtin)=cellstr(num2str(yas));
            clear  iTin1;
        else
            DataA(ptagActall,get(otherY,'Value'))=cellstr(num2str(yas));
        end
        
        %Set Flag and subflag
        if flag==3
            indl= ptagActall & ismember(DataA(:,vtype),type);%selected indices filtered for right type
            Set_subflag(find(indl),subflag,handles);
        end
        Calc_Fields_Callback(handles);
        feval(Currentplot,1,1,handles,1);
        Check_for_999(handles);
    end
end


 
function merge_butt_Callback(hObject, eventdata, handles)
global IMM DataA replaced dateok;
%     set(handles.figure1,'visible','off');
if ~dateok
    mmm=msgbox(sprintf('Year Day not calculated yet.\n\nYou need it to merge data.'),'Merge Error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end


activity(1,handles);
    
    merge;
    uiwait;
    col0=size(DataA,2);replaced=0;

    if iscell(IMM)
        for i=1:size(IMM,2)
            insert_in_data(IMM(:,i),'file',handles) %will be inserted as 'file' columns.Header will be checked for duplicates
%             samei=find(strcmpi(IMM(1,i),DataA(1,:))==1, 1);
%             if ~isempty(samei)
%                 ss=''; if length(samei)==1, ss='s'; end
%                 mmm=questdlg(sprintf('%s%s%s already exist%s.\nREPLACE or INSERT?',char(39),char(DataA(1,samei)),char(39),ss),'Data Import...','Replace','Insert','Cancel','Cancel');
%                 if strcmp(mmm,'Insert')
%                     IMM(1,i)=cellstr(mmm2);
%                     insert_in_data(IMM(:,i),'file',handles) %will be inserted as 'file' columns
% %                     %change name if header
% %                     options.Resize='on';options.WindowStyle='normal';
% %                     mmm2=inputdlg('Enter NEW Name for Column','Column Header Change',1.3,{[char(IMM(1,i)) '_2']},options);
% %                     if ~isempty(find(strcmpi(mmm2,DataA(1,:))==1, 1));% 2 columns with same header not allowed.
% %                         mmm=msgbox(sprintf('Column %s%s%s not imported',char(39),char(IMM(1,i)),char(39)),'Import Error','modal');
% %                         uiwait(mmm);
% %                     else   
% %                         IMM(1,i)=cellstr(mmm2);
% %                         insert_in_data(IMM(:,i),'file',handles) %will be inserted as 'file' columns
% %                     end
%                 elseif strcmp(mmm,'Replace')
%                     DataA(:,samei)=IMM(:,i);j=j+1;
%                 end
%             else
%                 insert_in_data(IMM(:,i),'file',handles); %will be inserted as 'file' columns
%             end
        end
    end
    col1=size(DataA,2);
    if (col1-col0)+replaced>0
        msgbox(sprintf('%d column(s) were imported!',col1-col0+replaced),'Import Result','modal');        
    else
        msgbox(sprintf('NO Column Imported!'),'Import Result','modal');
    end
    DataA=strtrim(DataA);DataA(strcmp(DataA,'NaN'))={'-999'};
    clear replaced;
activity(0,handles);




function insert_in_data(IA,aswhat,handles)
global DataA listf listf_file listf_calc Nonecol nvar replaced;
%Adds columns of data in DataA ONE AT A TIME
%either at the end of the 'data' columns (columns imported from data file)
%or at the end of the 'calculated data' columns (after 'data' cols and before 'none' which is last)

adds=size(IA,2);%adds = column in array to add
if adds>1, return; end
ss='';mmm2='';

if strcmp(aswhat,'file')
%Checks if header exists already
samei=find(strcmpi(IA(1,1),DataA(1,:))==1);
while (size(samei,2)>1)% fix multiple instances of headers
    mmm=questdlg(sprintf(['This header (%s) is already present more than once.' '\n\nThis could create issues with the program.' ...
                                        '\n\nDo you want the program to AUTOCORRECT this' '\n\nOr correct it yourself and RE-IMPORT the data?'],IA{1,1}),...
                                        'Duplicate Header','Autocorrect','Re-import','Autocorrect');
    if strcmpi(mmm,'Re-import'),return;end
    if max(samei)>length(listf_file) %one duplicate is a calculated column - should be last one (not renumbered)
        for k=1:size(samei,2)-1, DataA(1,samei(k)) = {[DataA{1,samei(k)} '_' strtrim(num2str(size(samei,2)+1))]}; end   % col renumbered starting at number of occurrences+1
    else
        for k=2:size(samei,2), DataA(1,samei(k)) = {[DataA{1,samei(k)} '_' strtrim(num2str(size(samei,2)))]}; end   % col renumbered starting at number of occurrences+1
    end
    samei=find(strcmpi(IA(1,1),DataA(1,:))==1);
end
if (size(samei,2)==1)
    mmm=questdlg(sprintf('%s%s%s already exist%s.\nREPLACE or INSERT?',char(39),char(DataA(1,samei)),char(39),ss),...
        'Data Import...','Replace','Insert','Cancel','Cancel');
    
    if strcmp(mmm,'Cancel'),return;end
    if strcmp(mmm,'Replace'),DataA(:,samei)=IA(:,1);replaced=replaced+1;return;end    
end

while ~isempty(samei) & strcmpi(mmm2,'Cancel')~=1 %will loop until name is unique or user cancelled
    options.Resize='on';options.WindowStyle='normal';
    k=2;while (k<=10 & ~isempty(samei)),samei=find(strcmpi([IA{1,1} '_' num2str(k)],DataA(1,:))==1); if ~isempty(samei),k=k+1;end,end
    if k>10 k=0;end
    mmm2 = inputdlg('Enter NEW Header',sprintf('Existing Header (%s)',IA{1,1}),1.3,{[IA{1,1} '_' num2str(k)]},options);
    samei = find(strcmpi(mmm2,DataA(1,:))==1);
    IA(1,1) = mmm2;
end
if strcmpi(mmm2,'Cancel'),return;end

% if max(samei)>length(listf_file) %one duplicate is a calculated column - should be last one (not renumbered)
%     IA(1,1) = {[IA{1,1} '_' strtrim(num2str(size(samei,2)+1))]};    % col renumbered starting at number of occurrences+1
%     ss='It could conflict with a column created by the program.';
% else
%     IA(1,1) = {[IA{1,1} '_' strtrim(num2str(size(samei,2)+1))]};     % col renumbered starting at number of occurrences+1
% end
% mmm=msgbox(sprintf('This header already exists:\n\n%s\n%s\n\nIt has been modified\n\nto differentiate it.',dupcoltxt,ss),'Column Insert...','Warn','modal');
% uiwait(mmm);

end

switch  aswhat
    case 'file'% add data as a file column
        flim=length(listf_file);%limit between file column and calc columns
        
        DataA(:,flim+1+adds:end+adds)=DataA(:,flim+1:end);
        DataA(:,flim+adds)=IA(:,1);Nonecol=Nonecol+adds;
        nvar(find(nvar>flim))=nvar(find(nvar>flim))+adds;
        listf_file=DataA(1,1:flim+adds);%list of file columns.
        listf_file=strtrim(listf_file);
        
    otherwise% add data as a field column
        flim=length(listf)-1;%limit between calc column and none column
        if strfind(char(DataA(1,flim)),'Cruise Day')  %YDay and cruise day have already been calc
%                                   Replaces previous 'YDay Calc' and 'Cruise Day'
            if strfind(char(IA(1,1)),'YDay Calc')
               DataA(:,flim-1)=IA(:,1);%replace YDay
            else
               DataA(:,flim)=IA(:,1);%replace Cruise Day
            end
        else
            DataA(:,flim+1+adds:end+adds)=DataA(:,flim+1:end);
            DataA(:,flim+1)=IA(:,1);Nonecol=Nonecol+adds;
            nvar(find(nvar>flim))=nvar(find(nvar>flim))+adds;
        end    
              
end

listf_calc=DataA(1,length(listf_file)+1:end-1);%list of calc fields...after file columns and before 'none'
listf_calc=strtrim(listf_calc);
listf=DataA(1,:);
listf=strtrim(listf);

Update_List(handles);
Update_Var(handles);


function delete_col_butt_Callback(hObject, eventdata, handles)
%deletes a column from DataA

global DataA listf listf_file listf_calc Nonecol nvar;

[s,v] = listdlg('PromptString','Select Columns to delete:',...
                'SelectionMode','multiple',...
                'ListString',listf_file);
if ~v, return; end %exits if Cancel

for i=1:length(s)
    nvar(find(nvar==s(i)))=Nonecol;
    nvar(find(nvar>s(i)))=nvar(find(nvar>s(i)))-1;
    Nonecol=Nonecol-1;
    DataA(:,s(i):end-1)=DataA(:,s(i)+1:end); DataA(:,end)=[];%shifts Data left and deletes last column 
    s(i+1:end)=s(i+1:end)-1;
end
%DataB=DataA(:,1:end-length(s));
%clear DataA;DataA=DataB;clear DataB;

flim=length(listf_file)-length(s);
listf_file=DataA(1,1:flim);%list of file columns.
listf_file=strtrim(listf_file);
listf_calc=DataA(1,length(listf_file)+1:end-1);%list of calc fields...after file columns and before 'none'
listf_calc=strtrim(listf_calc);
listf=DataA(1,:);
listf=strtrim(listf);

Update_List(handles);
Update_Var(handles);
Calc_Fields_Callback(handles);

function activity(n,handles)

%n=0 means no activity. Activity_stat button set to green
%n=1 means program is active. activity_stat button text set to 'Busy' in bold red

if n==0
    set(handles.activity_stat,'BackgroundColor',get(handles.btnplot,'BackgroundColor'));
    set(handles.activity_stat,'String','Idle');

else
    set(handles.activity_stat,'BackgroundColor',get(handles.figure1,'Color'));
    set(handles.activity_stat,'String','Busy');

end


function QuestAirVal_Button_Callback(hObject, eventdata, handles, varargin)
      %button to flag air values
     Change_flag(handles.rf138, eventdata, handles, varargin);
     
function [x,y]=gpos(h_figure,h_axes)
%GPOS Get current position of cusor and return its coordinates in axes with handle h_axes
% h_axes - handle of specified axes
% [x,y]  - cursor coordinates in axes h_aexs
%
% -------------------------------------------------------------------------
% Note:
%  1. This function should be called in the figure callback WindowButtonMotionFcn.
%  2. It works like GINPUT provided by Matlab,but it traces the position
%       of cursor without click and is designed for 2-D axes.
%  3. It can also work even the units of figure and axes are inconsistent,
%       or the direction of axes is reversed.
% -------------------------------------------------------------------------

% Written by Kang Zhao,DLUT,Dalian,CHINA. 2003-11-19
% E-mail:kangzhao@student.dlut.edu.cn

% h_figure=gcf;
%h_figure=handles.figure1;
units_figure = get(h_figure,'units');
units_axes   = get(h_axes,'units');

if_units_consistent = 1;

if ~strcmp(units_figure,units_axes)
    if_units_consistent=0;
    set(h_axes,'units',units_figure); % To be sure that units of figure and axes are consistent
end

% Position of origin in figure [left bottom]
pos_axes_unitfig    = get(h_axes,'position');
width_axes_unitfig  = pos_axes_unitfig(3);
height_axes_unitfig = pos_axes_unitfig(4);

xDir_axes=get(h_axes,'XDir');
yDir_axes=get(h_axes,'YDir');

% Cursor position in figure
pos_cursor_unitfig = get( h_figure, 'currentpoint'); % [left bottom]

if strcmp(xDir_axes,'normal')
    left_origin_unitfig = pos_axes_unitfig(1);
    x_cursor2origin_unitfig = pos_cursor_unitfig(1) - left_origin_unitfig;
else
    left_origin = pos_axes_unitfig(1) + width_axes_unitfig;
    x_cursor2origin_unitfig = -( pos_cursor_unitfig(1) - left_origin_unitfig );
end

if strcmp(yDir_axes,'normal')
    bottom_origin_unitfig     = pos_axes_unitfig(2);
    y_cursor2origin_unitfig = pos_cursor_unitfig(2) - bottom_origin_unitfig;
else
    bottom_origin_unitfig = pos_axes_unitfig(2) + height_axes_unitfig;
    y_cursor2origin_unitfig = -( pos_cursor_unitfig(2) - bottom_origin_unitfig );
end

xlim_axes=get(h_axes,'XLim');
width_axes_unitaxes=xlim_axes(2)-xlim_axes(1);

ylim_axes=get(h_axes,'YLim');
height_axes_unitaxes=ylim_axes(2)-ylim_axes(1);

x = xlim_axes(1) + x_cursor2origin_unitfig / width_axes_unitfig * width_axes_unitaxes;
y = ylim_axes(1) + y_cursor2origin_unitfig / height_axes_unitfig * height_axes_unitaxes;

% Recover units of axes,if original units of figure and axes are not consistent.
if ~if_units_consistent
    set(h_axes,'units',units_axes); 
end

function err=cell2csv(handles,filename,cellArray,delimiter,mode)
global expot expod vgroup vship pi_names;
% Writes a character cell array content into a *.csv file.
%
% CELL2CSV(handles,filename,cellArray,delimiter,mode)
%
% filename = Name of the file to save. [ i.e. 'text.csv' ]
% cellarray = Name of the Char Cell Array containing the data
% delimiter = default:','    '\t' for tab
% mode = specifies the mode of opening the file. See fopen() for a detailed
% list (default is overwrite i.e. 'w')
%
% by Denis Pierrot, 2009

if nargin<4
    delimiter = ',';
end
if nargin<5
    mode = 'w';
end

err=0;
fid = fopen(filename,mode);
if fid<0
    err=1;
    return ;
end
if ispc
    nbfn=regexp(filename,'\\');
else
    nbfn=regexp(filename,filesep);
end
nfn=filename(nbfn(length(nbfn))+1:end);%Separates file name from path

%generate format string    
fstr='';
for i=1:size(cellArray,2)-1, fstr=[fstr '%s' delimiter];end
fstr=[fstr '%s\n'];

%progress bar
mess1=' ';
if ~isempty(strfind(filename,'_Final')), mess1=' final ';elseif ~isempty(strfind(filename,'_Working')),mess1=' full ';end
mess=['\fontname{courier}Saving' mess1 'csv Data...'];
h = waitbar(0,mess);
CenterWindow(handles,h,handles.figure1);

%text for SOCAT dashboard
if ~isempty(strfind(filename,'_Final'))
    fprintf(fid,'Expocode: %s\n',strcat(expot,expod));
    fprintf(fid,'Ship: %s\n',vship);
    fprintf(fid,'Group: %s\n',vgroup);
    fprintf(fid,'Investigators: %s\n',pi_names);
end

cellArray=cellArray';j=0;jm=size(cellArray,2);
for  i=1:size(cellArray,2)
    j=j+1;
    if j>500
        j=0;waitbar(i/jm,h,[mess '(' num2str(i) '/' num2str(jm) ')']);
    end
    fprintf(fid,fstr,cellArray{:,i});
end
close(h);
fclose(fid);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over versiont.
function versiont_ButtonDownFcn(hObject, eventdata, handles)
global headp;
% hObject    handle to versiont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
% Show EXE path:
    if isdeployed % Stand-alone mode.
        if ispc
            [status, result] = system('set PATH');
            projDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
            mess=sprintf('Application is deployed.\n\nPath of .exe file is:\n\n%s', projDir);
        elseif ismac
            projDir=ctfroot;
            mess=sprintf('Application is deployed.\n\nPath of .app file is:\n\n%s', projDir);
        end
    else % Running from MATLAB.
        if exist(mfilename,'file')==2
            projDir = strrep([mfilename('fullpath') '.m'],[filesep mfilename '.m'],'');
        else
            projDir =pwd;
        end   
        mess=sprintf('Running from MATLAB.\n\nPath of .m file is:\n\n%s', projDir);

    end

    if isdeployed, junk=ctfroot;else junk=matlabroot;end
    mmm=msgbox(sprintf('%s\n\nMatlab root:\n\n%s', mess,junk),'','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    
%     mmm = questdlg('Do you want to copy some files or folder to the Project Directory?','Move files','Move Files','Move Folder','No','No');
%     if strcmp(mmm,'Move Files')>0
%         [fname,fpath] = uigetfile([headp filesep '*.*'],'Select Files to Copy','MultiSelect','on');
%         if (fpath==0), return; end
%     elseif strcmp(mmm,'Move Folder')>0
%         fname=uigetdir([headp filesep],'Select Folder to Move');
%         if (fname==0), return; end
%     else
%         return;
%     end
%     
% 
%     junk=[];
%     if iscell(fname) % more than 1 file selected
%         for i=1:length(fname)
%             [status, junk]= copyfile([fpath fname{i}],[projDir filesep]);
%         end
%     else
%         if exist(fname,'dir')==7
%             [status, junk]= copyfile(fname,[projDir filesep]);
%         else
%             [status, junk]= copyfile([fpath fname],[projDir filesep]);
%         end 
%     end
%     if ~isempty(junk)
%         msgbox(junk,'','modal');
%         uiwait;
%     end


    

function Make_Data_Assignment_Box_Reappear(hObject, eventdata, handles)
%executes when clicking the 'Variable - File Column' text box
%Used to make the assignment box reappear when mistakes (like command-clicking
%the box) makes it disappear for the GUI

set(handles.lb,'Value',1.0);


function xml_create(hObject, eventdata, handles)
global datap resultp headp hfile sysinip sysinif where_xml;

nf=dir(where_xml);    
%check presence of necessary file
if ~ismember('xml.tsv.txt',{nf.name})
       mmm=questdlg(sprintf('could not locate the xml file: n\n%s\n in\n\n%s\nDo you want to locate it?','xml.tsv.txt',where_xml),'xml file error','Yes', 'Abort','Yes');
       if strcmp(mmm,'Yes')>0
           [fname,where_xml] = uigetfile('*.txt','Select xml Data  File',[where_xml filesep 'xml.tsv.txt']);
           if isempty(where_xml), return; end
           if strcmp(where_xml(end),filesep)>0,where_xml=where_xml(1:end-1);end
           Save_Ini(datap,resultp,headp,hfile,sysinip,sysinif,where_xml);
       else
           return;
       end
end


Edit_xml;
uiwait;


% --------------------------------------------------------------------
function Transform_Data_Callback(hObject, eventdata, handles)
global DataA listf_file listf_calc listf;
global vlicorx vlicorw Nonecol;
% hObject    handle to td1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

activity(1,handles); %drawnow;

gs=get(handles.popGselect,'String');gv=get(handles.popGselect,'Value');
if strcmp(gs(1,:),'Graph 1')==1 & strcmp(gs(2,:),'Graph 2')==1
    mmm=msgbox('No data plotted','Transformation','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

options.Resize='on';
mmm=inputdlg({'New Column Header:','Equation y=log10(3.5*x^2) + exp(x^(1/2))'},'Data Transformation',1,{[gs{gv,:}, ' mod'],'y='},options);
if isempty(mmm), return;end
hv=[handles.LeftY,handles.RightY];
op={'*';'/';'\^'};
try
    if strcmpi(mmm{2},'LI7000') %will perform re-calc of xCO2 to correct for firmware 2.0 bug
        mm='Use xCO2';
        if (Nonecol==vlicorx | Nonecol==vlicorw),mm=msgbox(sprintf('Either xCO2 or xH20 not assigned.\n\nPlease correct.'),'modal');uiwait(mm);return;end
        if get(hv(gv),'Value')~= vlicorx,  mm=questdlg(sprintf('Only xCO2 directly output from LICOR\n\ncan be re-calculated.'),'LI-7000 Bug', 'Use xCO2','Cancel','Use xCO2');end
        if ~strcmp(mm,'Use xCO2'), return;end

        
        x=str2num(char(DataA(2:end,vlicorx)));x(x==-999)=NaN;
        xw=str2num(char(DataA(2:end,vlicorw)));xw(xw==-999)=NaN;
        x(x<0)=0;xw(xw<0)=0; %equations not defined for negative values
        
        y=x;%initialization
        y(x<700)=x(x<700)-(3.0277e-07).*x(x<700).^2.0482.*xw(x<700).^0.99154;
        y(x>=700)=x(x>=700)-(3.6781e-06).*x(x>=700).^1.6707.*xw(x>=700).^0.97429;
        
        y(2:end+1)=y(1:end);y=cellstr(num2str(y));
        y{1}=mmm{1};
        y=strtrim(y);y(strcmp(y,'NaN'))={'-999'};
        insert_in_data(y,'file',handles);
        mm=msgbox(sprintf('Now you can re-assign the variable ''LICOR xCO2'' to ''%s''\n\nAnd re-do all the calculations.',mmm{1}),'modal');uiwait(mm);
        clear x xw;
        
    else     
        x=str2num(char(DataA(2:end,get(hv(gv),'Value'))));x(x==-999)=NaN;
        for i=1:length(op),mmm{2}=regexprep(mmm{2},op{i},['.' op{i}]); end  % change to "array operations"
        mmm{2}=regexprep(mmm{2},'y','y(:)'); mmm{2}=regexprep(mmm{2},'x','x(:)');
        y=x;%initializes y
        eval([mmm{2},';']);%evaluate equation    

        y(2:end+1)=y(1:end);y=cellstr(num2str(y));y{1}=mmm{1};  %insert header  
        y=strtrim(y);y(strcmp(y,'NaN'))={'-999'};
        insert_in_data(y,'file',handles);
end
    
   
    
catch
        disp(lasterr);
        if strcmp(lasterr,'Attempt to reference field of non-structure array.')~=1
            mmm=msgbox(sprintf('Error has occurred\n(%s)',lasterr),'Error','error','modal');
            CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
        end
        drawnow;

end



activity(0,handles);


% --- Executes on button press in Print_Graph.
function Print_Graph_Callback(hObject, eventdata, handles)
% hObject    handle to Print_Graph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global resultp expot expod hcb;

% need to create a new figure, recreate the axes on it, redo the plot and
% print the whole thing.
p3d=get(handles.p3d,'Value');
% fig = get(0, 'CurrentFigure');
% export_fig test0.jpg;

f2=figure;movegui(f2,'center');
if ~p3d, nax2=copyobj(handles.axes2,f2);  nax1=copyobj(handles.axes1,f2); else nax1=copyobj([handles.axes1,hcb],f2);colormap(nax1(1),jet);end
set(nax1(1),'outerposition',[0 0 1 1]);pos1=get(nax1(1),'Position');
if ~p3d, set(nax2,'position',pos1);end

t=inputdlg({'File Name (no ext)','Input Title for figure'},'Graph Title',1,{['pCO2_Sys_plot'],strcat(expot,expod)});
if size(t,1)>=1, title(nax1(1),t{2});end
if exist([resultp filesep ],'dir')~=7,  resultp=uigetdir([resultp filesep],'Select Folder to print');end

%set(f2,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
print (f2,'-djpeg','-r300', [resultp filesep t{1} '.jpg']);
close (f2);


% --------------------------------------------------------------------
function Replace_Data_Callback(hObject, eventdata, handles)
% hObject    handle to rd2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DataA ind1 ind2 hp hp2 Currentplot dotsize;
global   vteq vtin vtini vsal vpeqa vpatm vpamb vtype listf_file sub_flag;

% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';...
%         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';...
%         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';... %         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';... %         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'interpolated using questionable SST'};  variable missing';'interpolated using questionable SST'};


activity(1,handles); %drawnow;

equ_atm={'EQU';'ATM';'EQU-DRAIN';'ATM-DRAIN'};
equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};

mmm=0;
gs=get(handles.popGselect,'String');
if strcmp(gs(1,:),'Graph 1')==1 & strcmp(gs(2,:),'Graph 2')==1
    mmm=msgbox('No data plotted','Replacement','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

ax_fl=get(hObject,'Tag');  %Tag is 'rd'+axe #
axstr=ax_fl(3);

if strcmp(axstr,'1')==1
    ax=handles.axes1;   yax=handles.LeftY;     ind=ind1;     chp=hp;
else
    ax=handles.axes2;   yax=handles.RightY;    ind=ind2;     chp=hp2;
end
xlimits = get(ax,'XLim');    ylimits = get(ax,'YLim');
xp = (xlimits(2)-xlimits(1))/100;    yp = (ylimits(2)-ylimits(1))/100;

%     if get(yax,'Value')~=vteq & get(yax,'Value')~=vtin & get(yax,'Value')~=vlicorcav & get(yax,'Value')~=vsal
if get(yax,'Value')>length(listf_file)
    mmm=msgbox(sprintf('Can Only Replace Original Data.\n\n%s is Calculated.\n\nAction Cancelled',char(DataA(1,get(yax,'Value')))),'Replacement Error','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm); activity(0,handles);
    return;
end

rh=findobj(ax,'Type','Rectangle');
if isempty(rh)
    mmm=msgbox(sprintf('No area detected,\n\nSelect Data first.'),'Data Replacement','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm); activity(0,handles);
    return;
end

ptatoiall=[];
for r=1:length(rh)% for each rectangle
    ptatoi=[];
    rv=get(rh(r),'Position');
    xv=[rv(1) rv(1) rv(1)+rv(3) rv(1)+rv(3)];
    yv=[rv(2) rv(2)+rv(4) rv(2)+rv(4) rv(2)];
    %index of points to replace
    ptatoi = inpolygon(str2num(char(DataA(2:end,handles.CurrentX.Value))),str2num(char(DataA(2:end,get(yax,'Value')))),xv,yv);
    ptatoi=find(ptatoi)+1;
    indall=[];
    for i=1:size(ind,1), indall=union(indall,ind{i});end % take only plotted points
    ptatoi=intersect(ptatoi,indall);clear indall;
    ptatoiall=union(ptatoiall,ptatoi);
end

if mmm==0 & ~isempty(ptatoiall)  %if area selected, replaces data
        
    %--------------------------------------------------------------------------
    % Creates interface
    %--------------------------------------------------------------------------
    [scrsz]=get(0,'ScreenSize');   wd=400;ht=200;spc=50;ctrln=13;
    
    fig1 = figure('position', [(scrsz(3)-wd)/2 (scrsz(4)-ht)/2 wd ht],'menubar','none','numbertitle','off');
    bg=uibuttongroup('Position',[0 0.25 1 0.75]);
    uicontrol (bg,'Style','text','position', [(wd-200)/2 0.75*ht-30 200 20],'String','Replace Data With...','FontWeight','bold','backgroundcolor',get(fig1,'color'));
    b(1) = uicontrol (bg,'Style','radiobutton','position', [30 0.75*ht-85 100 50],'String','Constant Value','backgroundcolor',get(fig1,'color'));
    b(2) = uicontrol (bg,'Style','radiobutton','position', [30 0.75*ht-135 100 50],'String','Column Data','backgroundcolor',get(fig1,'color'));
    set(b(2),'Value',1);
    r(1) = uicontrol (bg,'Style','edit','position', [200 0.75*ht-85 100 20],'String','Enter Value here','backgroundcolor',get(fig1,'color'));
    r(2) = uicontrol (bg,'Style','popupmenu','position', [200 0.75*ht-135 150 20],'String',listf_file,'backgroundcolor',get(fig1,'color'));
    ccb=uicontrol ('Style','togglebutton','position', [(wd-100) 20 40 20],'String','Cancel','Callback','uiresume(gcbf)');
    okb=uicontrol ('Style','togglebutton','position', [(wd-50) 20 30 20],'String','OK','Callback','uiresume(gcbf)');
    align([b(1) r(1)],'None','Middle');  align([b(2) r(2)],'None','Middle');   align([r(1) r(2)],'Center','None');
    uicontrol(okb);
    uiwait(fig1);
    
    set(fig1,'visible','off');
    if ~get(okb,'value'), return;end
    
    xout=str2num(char(DataA(ptatoiall,handles.CurrentX.Value)));
    yorig=str2num(char(DataA(ptatoiall,get(yax,'Value'))));
    % use good point less those to interp
    xin=str2num(char(DataA(setdiff(ind{1},ptatoiall),handles.CurrentX.Value)));
    yin=str2num(char(DataA(setdiff(ind{1},ptatoiall),get(yax,'Value'))));
    
    if get(b(1),'Value') %constant value
        if isempty(str2num(get(r(1),'String'))),mmm=msgbox(sprintf('Value entered is not recognized as a number.\n\nTry again.'),'Error');uiwait(mmm);return;end
        yout=zeros(size(xout,1),1)+str2num(get(r(1),'String')); %create array filled with constant
    else %file column
        yout=str2num(char(DataA(ptatoiall,get(r(2),'Value'))));
    end
    close(fig1);
    
    if length(yout)==0 | length(xout)==0,mmm=msgbox(sprintf('No Data Replaced.\n\nTry again.'),'Error');uiwait(mmm);return;end
    
    wd=500;ht=500;
    fig3 = figure('position', [(scrsz(3)-wd)/2 (scrsz(4)-ht)/2 wd ht],'numbertitle','off');
    plot(xin,yin,'ok',xout,yout,'or',xout,yorig,'ob','MarkerSize',dotsize);
    set(gca,'XLim',xlimits);%set(gca,'YLim',ylimits);
    legend('Good Data','New Data','Replaced Data','Location','Best');
    axes(gca);
    
    set(gcf,'Toolbar','figure')
    axl=0.15;axb=0.2;axw=0.75;axh=0.75;
    gcfp=get(gcf,'position'); set(gca,'position',[axl axb axw axh]);
    fl=gcfp(1);fb=gcfp(2);fw=gcfp(3);fh=gcfp(4);
    butw=fw/8;buth=0.1*fh;clear gcbf;
    uicontrol('Position',[(fw-butw)/2 2 butw buth],'String','Continue...',...
        'Callback','uiresume');
    uiwait(fig3);
    close(fig3);
    
    mmm=questdlg(sprintf('Accept Replacement?'),'Replacement Results','Yes','No','No');
    
    if strcmp(mmm,'Yes')==1
    
        type='';types='';
        switch get(yax,'Value')
            case {vtin,vtini}
                type=equ;types='EQU';
            case vteq
                type=equ;types='EQU';
            case vpamb
                type=equ_atm;types='EQU and ATM';
            case vpatm
                type=atm;types='ATM';
            case vsal
                type=equ;types='EQU';
            otherwise
                type=equ_atm;types='EQU and ATM';
        end

        DataA(ptatoiall,get(yax,'Value'))=cellstr(num2str(yout));
        tf=ismember([vtin vtini vteq vsal vpeqa vpatm vpamb],get(yax,'Value'));
        if sum(tf)>0
            mmm=questdlg(sprintf('FLAG DATA?\n\nAUTO = program does it.\nMANUAL = User chooses'),'Replacement Flag','AUTO (3)', 'MANUAL', 'NO FLAG','AUTO (3)');
            if strcmp(mmm,'AUTO (3)')==1
                tf=ismember([vtin vtini vteq vsal],get(yax,'Value'));
                if sum(tf)>0
                    if sum(tf)>1
                        mmm=questdlg(sprintf(['These Data Are Assigned to at Least 2 of the Following Variables:'...
                            '\n\n\"in situ Temp\"\n\"EQU Temp\"\n\"Salinity\"\n\nHow Do You Want to Flag Them?'...
                            '\n\n2 - %s\n3 - %s\n5 - %s'],sub_flag{2,:},sub_flag{3,:},sub_flag{5,:}),'Interpolation Flag','2', '3', '5', '2');
                        flv=int8(str2num(mmm));
                    else
                        flv=2*(ismember(get(yax,'Value'),[vtin vtini])) + 3*(get(yax,'Value')==vteq)+5*(get(yax,'Value')==vsal);
                    end
                    
                    toflagi=intersect(find(ismember(DataA(:,vtype),type)),ptatoiall(~isnan(yout)));
                    Set_subflag(toflagi,flv,handles);
                end
                
                tf=ismember([vpeqa vpatm vpamb],get(yax,'Value'));
                if sum(tf)>0
                    flv=6; %Questionable/Interpolated P
                    toflagi=intersect(find(ismember(DataA(:,vtype),type)),ptatoiall(~isnan(yout)));
                    Set_subflag(toflagi,flv,handles);
                end
            elseif strcmp(mmm,'MANUAL')==1
                %--------------------------------------------------------------------------
                % Creates interface
                %--------------------------------------------------------------------------
                [scrsz]=get(0,'ScreenSize');   wd=400;ht=200;spc=50;
                
                fig2 = figure('position', [(scrsz(3)-wd)/2 (scrsz(4)-ht)/2 wd ht],'menubar','none','numbertitle','off');
                uicontrol ('Style','text','position', [(wd-200)/2 ht-30 200 20],'String','Select Subflag to Erase','FontWeight','bold','backgroundcolor',get(fig2,'color'));
                r(1) = uicontrol ('Style','popupmenu','position', [(wd-250)/2 ht-60 250 20],'String',sub_flag,'backgroundcolor',get(fig2,'color'));
                ccb=uicontrol ('Style','togglebutton','position', [(wd-100) 20 40 20],'String','Cancel','Callback','uiresume(gcbf)');
                okb=uicontrol ('Style','togglebutton','position', [(wd-50) 20 30 20],'String','OK','Callback','uiresume(gcbf)');
                uicontrol(okb);
                uiwait(fig2);
                
                %set(fig2,'visible','off');
                if ~get(okb,'value'), return;end
                flv=min(get(r(1),'Value'),10);
                toflagi=intersect(find(ismember(DataA(:,vtype),type)),ptatoiall(~isnan(yout)));
                if flv>10,  Set_subflag(toflagi,10,handles); end
                Set_subflag(toflagi,flv,handles);
                close(fig2);
            end
        end
        Calc_Fields_Callback(handles,ptatoiall(~isnan(yout)));
        Stat_Update_Callback(1, 1, handles);
        feval(Currentplot,1,1,handles,1);
    end
end

if isempty(ptatoiall)
    mmm=msgbox (sprintf('No Data to replace.\n\nCheck Data or Flags'),'Replacement Error!','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

Check_for_999(handles);
activity(0,handles);


%------------------------------------------------------------------------    
function Erase_Subflag(hObject, eventdata, handles)
global DataA ind1 ind2 hp hp2 Currentplot;
global sub_flag;
global vflag vsubf vsubfu ;

% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';...
%         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';...
%         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';... %         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';... %         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'interpolated using questionable SST'};  variable missing';'interpolated using questionable SST'};

activity(1,handles); %drawnow;

mmm=0;
gs=get(handles.popGselect,'String');
if strcmp(gs(1,:),'Graph 1')==1 & strcmp(gs(2,:),'Graph 2')==1
    mmm=msgbox('No data plotted','Flag Replacement','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    return;
end

ax_fl=get(hObject,'Tag');  %Tag is 'esfl'+axe #
axstr=ax_fl(5);

if strcmp(axstr,'1')==1
    ax=handles.axes1;   yax=handles.LeftY;     ind=ind1;     chp=hp;
else
    ax=handles.axes2;   yax=handles.RightY;    ind=ind2;     chp=hp2;
end
xlimits = get(ax,'XLim');    ylimits = get(ax,'YLim');
xp = (xlimits(2)-xlimits(1))/100;    yp = (ylimits(2)-ylimits(1))/100;

rh=findobj(ax,'Type','Rectangle');
if isempty(rh)
    mmm=msgbox(sprintf('No area detected,\n\nSelect Data first.'),'Flag Replacement','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm); activity(0,handles);
    return;
end

ptatoiall=[];
for r=1:length(rh)% for each rectangle
    ptatoi=[];
    rv=get(rh(r),'Position');
    xv=[rv(1) rv(1) rv(1)+rv(3) rv(1)+rv(3)];
    yv=[rv(2) rv(2)+rv(4) rv(2)+rv(4) rv(2)];
    %index of points to replace
    ptatoi = inpolygon(str2num(char(DataA(2:end,handles.CurrentX.Value))),str2num(char(DataA(2:end,get(yax,'Value')))),xv,yv);
    ptatoi=find(ptatoi)+1;
    indall=[];
    for i=1:size(ind,1), indall=union(indall,ind{i});end % take only plotted points
    ptatoi=intersect(ptatoi,indall);clear indall;
    ptatoiall=union(ptatoiall,ptatoi);
end


if mmm==0 & ~isempty(ptatoiall)  %if area selected, replaces data
        
    %--------------------------------------------------------------------------
    % Creates interface
    %--------------------------------------------------------------------------
    [scrsz]=get(0,'ScreenSize');   wd=400;ht=200;spc=50;
    
    fig2 = figure('position', [(scrsz(3)-wd)/2 (scrsz(4)-ht)/2 wd ht],'menubar','none','numbertitle','off');
    uicontrol ('Style','text','position', [(wd-200)/2 ht-30 200 20],'String','Select Subflag to Erase','FontWeight','bold','backgroundcolor',get(fig2,'color'));
    r(1) = uicontrol ('Style','popupmenu','position', [(wd-250)/2 ht-60 250 20],'String',sub_flag,'backgroundcolor',get(fig2,'color'));
    ccb=uicontrol ('Style','togglebutton','position', [(wd-100) 20 40 20],'String','Cancel','Callback','uiresume(gcbf)');
    okb=uicontrol ('Style','togglebutton','position', [(wd-50) 20 30 20],'String','OK','Callback','uiresume(gcbf)');
    uicontrol(okb);
    uiwait(fig2);
    
    %set(fig2,'visible','off');
    if ~get(okb,'value'), return;end
    sufl2er2=sub_flag{get(r(1),'Value')};
    sufl2er=sub_flag{min(get(r(1),'Value'),10)};
    %replaces special characters ()+
    sufl2er=regexprep(sufl2er,'\(','\\\(');   sufl2er=regexprep(sufl2er,'\+','\\\+');   sufl2er=regexprep(sufl2er,'\)','\\\)');
    sufl2er2=regexprep(sufl2er2,'\(','\\\(');   sufl2er2=regexprep(sufl2er2,'\+','\\\+');   sufl2er2=regexprep(sufl2er2,'\)','\\\)');
    
    close(fig2);

    chfl=[];nfl=[];
    ta=ptatoiall(str2num(char(DataA(ptatoiall,vflag)))==3);%all indices of selection that are flag 3
    already2=regexp(DataA(ta,vsubfu),sufl2er2);
    already=regexp(DataA(ta,vsubfu),sufl2er);
    if sum((~cellfun(@isempty,already2)))==0, mmm=msgbox(sprintf('Subflag Not Found in Selected Data!\n\nTry again.'));uiwait(mmm);return;end
    for j=1:length(ta)
        if ~isempty(already{j})%subflag usr  present
            nfl=[nfl ta(j)];
            %deletes sub_flag
            if already{j}>1
                DataA(ta(j),vsubf)=regexprep(DataA(ta(j),vsubf),[';' sufl2er],''); 
            elseif ismember(';',char(DataA(ta(j),vsubf)))
                DataA(ta(j),vsubf)=regexprep(DataA(ta(j),vsubf),[sufl2er ';'],'');
            else
                DataA(ta(j),vsubf)=regexprep(DataA(ta(j),vsubf),[sufl2er],'');
            end
            %deletes sub_flag user
            if already2{j}>1
                DataA(ta(j),vsubfu)=regexprep(DataA(ta(j),vsubfu),[';' sufl2er2],''); 
            elseif ismember(';',char(DataA(ta(j),vsubfu)))
                DataA(ta(j),vsubfu)=regexprep(DataA(ta(j),vsubfu),[sufl2er2 ';'],'');
            else
                DataA(ta(j),vsubfu)=regexprep(DataA(ta(j),vsubfu),[sufl2er2],'');
            end
            if length(char(DataA(ta(j),vsubfu)))==0, DataA(ta(j),vflag)={'2'};chfl=[chfl ta(j)];end
            
        end
    end
    Stat_Update_Callback(1, 1, handles);
    feval(Currentplot,1,1,handles,1);

    mmm=msgbox(sprintf('%d Subflags Found.\n\n%d Flags 3 Changed Back to 2.',length(nfl),length(chfl)),'Erase Flags');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm); activity(0,handles);
end


% --- Executes on button press in rdodryxco2.
function rdodryxco2_Callback(hObject, eventdata, handles)
% hObject    handle to rdodryxco2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdodryxco2
global off_on;

set([handles.xco2d_f,handles.xco2d_ft],'Visible',off_on(get(hObject,'Value')+1,:));
%set(handles.xco2d_ft,'Visible',off_on(get(hObject,'Value')+1,:));


function CenterWindow(handles,wh,ph)
global monitors;
%wh is handle of window to center
%ph is handle of parent or main window
pd=10^30;
m=1;w=getpixelposition(wh);p=getpixelposition(ph);
if size(monitors,1)>1
    for mi=1:size(monitors,1)
        if (p(1)+0.5*p(3)- monitors(mi,1))>0 && (p(1)+0.5*p(3) - monitors(mi,1))<pd
            m=mi ;pd=(p(1) - monitors(mi,1));
        end
    end
end
setpixelposition(wh,[monitors(m,1)+monitors(m,3)/2-w(3)/2 monitors(m,4)/2-w(4)/2 w(3) w(4)]);


% --- Executes on slider movement.
function map_res_Callback(hObject, eventdata, handles)
global map_folder;
% hObject    handle to map_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

map_res_s=['lo ';'med';'hi ']; map_s=['l';'h';'f']; 
handles.txt_map_res.String=map_res_s(hObject.Value,:); drawnow;%Slider has value of 1, 2 or 3
if isempty(dir([map_folder filesep 'gshhs_' map_s(hObject.Value,:) '.b']))
    mmm=msgbox(sprintf('Corresponding Map file (%s) NOT FOUND!\n Select another resolution!',['gshhs_' map_s(hObject.Value,:) '.b']));
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
end


% --------------------------------------------------------------------
function tblstd_modify(hObject, eventdata, handles)
% hObject    handle to remstd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stdv DataA nvar Nonecol vstdx vtype listf listf_calc listf_file;

try

    %save conditions in case of fail
    nvar0=nvar; Nonecol0=Nonecol; stdv0=stdv; DataA0=DataA; listf0=listf;
    listf_calc0=listf_calc;listf_file0=listf_file;
     
    %if ~isempty(stdv)
        stdvl=get(handles.tblstdv,'Data');
        if strcmp(hObject.Tag,'addstd'),junk='add';nm=length(stdvl)+1;
        elseif strcmp(hObject.Tag,'remstd'),junk='remove';nm=length(stdvl);end
        
        if strcmp(hObject.Tag,'addstd') | strcmp(hObject.Tag,'remstd')
            %string function doesn't work before 2016b
            if verLessThan('matlab','9.1'),str='';for i=1:nm, str=[str; num2str(i)];end
            else, str = string(1:nm); end
        
            [s,v] = listdlg('PromptString',['Select Std to ' junk ':'],'SelectionMode','single', 'ListString',str);
            if v==0, return;end  % no selection
        end
        
        if strcmp(hObject.Tag,'addstd'), stdvl(s+1:end+1)=stdvl(s:end); stdvl(s)=-999;
        elseif strcmp(hObject.Tag,'remstd')
            if s<length(stdvl),stdvl(s:end-1)=stdvl(s+1:end);end
            stdvl=stdvl(1:end-1);
        end
        %Update DataA columns
        if ~isempty(DataA) %data has been loaded already
            switch hObject.Tag
                case 'addstd'
                    for i=length(stdvl)-1:-1:s
                        %temptit={['Std' num2str(i) ' values']};
                        DataA(1,:)=regexprep(DataA(1,:),['Std' num2str(i) ' values'],['Std' num2str(i+1) ' values']);
                    end
                    n=find(strcmpi(['Std' num2str(s+1) ' values'],DataA(1,:)));
                    if isempty(n),n=find(strcmpi(['Std' num2str(s-1) ' values'],DataA(1,:)))+1;end  %adding last std
                    DataA(:,n+1:end+1)=DataA(:,n:end);
                    DataA{1,n}=['Std' num2str(s) ' values'];
                    DataA(2:end,n)={'-999'};
                    nnvar=find(nvar==n);%nnvar= nvar BEFORE which new nvar will be added
                    nvar(nnvar+1:end+1)=nvar(nnvar:end);
                    nvar(nvar>=n)=nvar(nvar>=n)+1;Nonecol=Nonecol+1;nvar(nnvar)=n;
                    mess=sprintf('One or more Std(s) have been added.\n\nStd Offsets should be re-calculated.\n\nCorrect and Re-Calc Fields');
                    cc='r';
                    set(handles.stdo_f,'BackgroundColor',cc,'TooltipString',mess);
                    mmm=msgbox(sprintf('Don''t forget to re-calculate the ''fields''.\nafter entering the STD value'),'STD addition','modal');
                    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                    
                case 'remstd'
                    n=find(strcmpi(['Std' num2str(s) ' values'],DataA(1,:)));
                    if isempty(n)
                        mmm=msgbox(sprintf('STD %d not found',s),'Std removal Error','modal');
                        CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                        return;
                    else
                        DataA(:,n:end-1)=DataA(:,n+1:end);
                        DataA=DataA(:,1:end-1);
                        nnvar=find(nvar==n);%nnvar=which nvar refers to col to be deleted
                        nvar(nnvar:end-1)=nvar(nnvar+1:end);nvar=nvar(1:end-1);
                        nvar(nvar>n)=nvar(nvar>n)-1;Nonecol=Nonecol-1;
                    end
                    
                case 'updatestd'
                    
                    if verLessThan('matlab','9.1'),stdr=find(~cellfun(@isempty,strfind(DataA(:,vtype),'STD')));
                    else, stdr=find(contains(DataA(:,vtype),'STD'));end
                    stds=unique(DataA(stdr,vtype));stdn=sort(unique(regexprep(stds,'\D','')));
                    stdvl=[]; stdm=max(str2num(char(stdn)));
                    for s=1:stdm
                        junk=find(strcmp(DataA(:,vtype),['STD' stdn{s}]),1);
                        if ~isempty(junk),stdvl=[stdvl;str2num(DataA{junk,vstdx})];
                        else
                            junk=find(strcmp(DataA(:,vtype),['STD' stdn{s} 'z']),1);
                            if ~isempty(junk),stdvl=[stdvl;0];
                            else
                                stdvl=[stdvl;-999];end
                        end
                        
                        n=find(strcmpi(['Std' num2str(s) ' values'],DataA(1,:)),1);
                        if isempty(n) %add stdx values column 
                            for i=stdm-1:-1:s %for all but last std
                                DataA(1,:)=regexprep(DataA(1,:),['Std' num2str(i) ' values'],['Std' num2str(i+1) ' values']);
                            end
                            k=1;
                            while isempty(n) && (s+k)<=stdm
                                n=find(strcmpi(['Std' num2str(s+k) ' values'],DataA(1,:)));
                                if isempty(n),k=k+1;end
                            end
                            if isempty(n),k=1;end
                            while isempty(n) && (s-k)>0%adding last std
                                n=find(strcmpi(['Std' num2str(s-k) ' values'],DataA(1,:)))+1;
                                if isempty(n),k=k+1;end
                            end
                            
                            if isempty(n)
                                mmm=msgbox(sprintf('NO ''STDx values'' found\n\nCheck your data.'),'Std update Error','modal');
                                CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                                return;
                            end
                            
                            DataA(:,n+1:end+1)=DataA(:,n:end);
                            DataA{1,n}=['Std' num2str(s) ' values'];
                            DataA(2:end,n)={'-999'};
                            nnvar=find(nvar==n);%nnvar= nvar BEFORE which new nvar will be added
                            nvar(nnvar+1:end+1)=nvar(nnvar:end);
                            nvar(nvar>=n)=nvar(nvar>=n)+1;Nonecol=Nonecol+1;nvar(nnvar)=n;
                        end
                    end
                    %search for any 'stdx values' column other than those related to STDx type
                    stdvc=find(~cellfun(@isempty,regexp(DataA(1,:),'Std[0-9]+ values')));
                    stdsc=unique(DataA(1,stdvc));stdnc=sort(unique(regexprep(stdsc,'\D','')));
                    for dc=stdnc(~ismember(stdnc,stdn)) %stdnc not found in stdn - delete these columns
                        n=find(strcmpi(['Std' dc{:} ' values'],DataA(1,:)));
                        DataA(:,n:end-1)=DataA(:,n+1:end);
                        DataA=DataA(:,1:end-1);
                        nnvar=find(nvar==n);%nnvar=which nvar refers to col to be deleted
                        nvar(nnvar:end-1)=nvar(nnvar+1:end);nvar=nvar(1:end-1);
                        nvar(nvar>n)=nvar(nvar>n)-1;Nonecol=Nonecol-1;
                    end
                    
            end
            
            stdv=stdvl;
            listf_calc=DataA(1,length(listf_file)+1:end-1);%list of calc fields...after file columns and before 'none'
            listf_calc=strtrim(listf_calc);
            listf=DataA(1,:);
            listf=strtrim(listf);
            
            Update_List(handles);
            Update_Var(handles);
           
        end
        
              
        %Update display
        set(handles.tblstdv,'Data',stdvl,'ColumnFormat',({'bank'}));
        tts=''; for i=1:length(stdvl), tts= [tts 'STD ' num2str(i) sprintf(' : %0.2f ppm',stdvl(i))];if i<length(stdvl), tts=[tts sprintf('\n')]; end, end
        set(handles.tblstdv,'TooltipString',tts);
        stdl={};
        for i=1:length(stdvl)  stdl(i,:)={true,['STD ' num2str(i,0)]};end
        tts=''; for i=1:length(stdvl), tts= [tts 'STD ' num2str(i)];if i<length(stdvl), tts=[tts sprintf('\n')]; end, end
        set(handles.tblstdl,'Data',stdl,'TooltipString',tts);
        stdv=stdvl;
        
        if strcmp(hObject.Tag,'updatestd'),  Calc_Fields_Callback(handles);end
        
        clear s v i k n nnvar mess cc stdr stds stdn stdm stdvc stdsc stdnc dc;
    %end
    
catch
    
    mmm=msgbox(lasterr,'Error');

    CenterWindow(handles,mmm,handles.figure1)
    uiwait(mmm);
    %restore conditions in case of fail
    if exist('Data0','var')      
        nvar=nvar0; Nonecol=Nonecol0; stdv=stdv0; DataA=DataA0;
        listf=listf0;listf_calc=listf_calc0;listf_file=listf_file0;
    end
    
    clear nvar0 Nonecol0 stdv0 DataA0 listf0 listf_calc0 listf_file0 ;
    
end


function tblunk_modify(hObject, eventdata, handles)
% hObject    handle to remstd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DataA  vtype ;

try


        unkl=get(handles.tblunkl,'Data');
        unknl=size(unkl,1);
        
        if strcmp(hObject.Tag,'addunk'), unknl=unknl+1;
        elseif strcmp(hObject.Tag,'remunk')
            if unknl>0,unknl=unknl-1;end
        end
        
        
        if strcmp(hObject.Tag,'updateunk')
            
            if verLessThan('matlab','9.1'),unkr=find(~cellfun(@isempty,strfind(DataA(:,vtype),'UNK')));
            else, unkr=find(contains(DataA(:,vtype),'UNK'));end
            unknl=0;
            if ~isempty(unkr)
                unks=unique(DataA(unkr,vtype));unkn=sort(unique(regexprep(unks,'\D','')));
                unknl=max(str2num(char(unkn)));
            end
        end
        
        %Update display
        unksl={};
        for i=1:unknl  unksl(i,:)={true,['UNK ' num2str(i,0)]};end
        tts=''; for i=1:unknl, tts= [tts 'UNK ' num2str(i)];if i<unknl, tts=[tts sprintf('\n')]; end, end
        set(handles.tblunkl,'Data',unksl,'TooltipString',tts);
                
        clear s v i unkr unks unkn str unkl unknl;
   
    
catch
    
    mmm=msgbox(lasterr,'Error');

    CenterWindow(handles,mmm,handles.figure1)
    uiwait(mmm);
    
end






% --------------------------------------------------------------------
function Display_DataA(hObject, eventdata, handles)
% hObject    handle to Data File frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DataA 

if strcmp(handles.figure1.SelectionType, 'alt') %right-click
    if ~exist('DataA','var') || isempty(DataA)
        mmm=msgbox('''DataA'' array doesn''t exist yet or is empty.','display Error','error','modal');
        CenterWindow(handles,mmm,handles.figure1)
        uiwait(mmm);
        return;
    end
    if ~exist('f','var'),f=figure('Name','DataA Array');f.Position=f.Position*diag([1 1 0.5 0.5]);CenterWindow(handles,f,handles.figure1);end
    t=uitable('Parent',f,'Position',[20 20 f.Position(3)-40 f.Position(4)-40],'Data',DataA);
    t.Units='normalized';
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text330(# out of range).
function Display_OutOfRange(hObject, eventdata, handles)
% hObject    handle to text330 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  out_of_range ;

if strcmp(handles.figure1.SelectionType, 'alt') %right-click
    if ~exist('out_of_range','var') || isempty(out_of_range)
        mmm=msgbox('''out_of_range'' array doesn''t exist yet or is empty.','display Error','error','modal');
        CenterWindow(handles,mmm,handles.figure1)
        uiwait(mmm);
        return;
    end
    if ~exist('f','var'),f=figure('Name','out_of_range Array');f.Position=f.Position*diag([1 1 0.5 0.5]);CenterWindow(handles,f,handles.figure1);end
    t=uitable('Parent',f,'Position',[20 20 f.Position(3)-40 f.Position(4)-40],'Data',out_of_range);
    t.Units='normalized';
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text330(# outliers).
function Display_Outliers(hObject, eventdata, handles)
% hObject    handle to text331 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  outliers ;

if strcmp(handles.figure1.SelectionType, 'alt') %right-click
    if ~exist('outliers','var') || isempty(outliers)
        mmm=msgbox('''outliers'' array doesn''t exist yet or is empty.','display Error','error','modal');
        CenterWindow(handles,mmm,handles.figure1)
        uiwait(mmm);
        return;
    end
    if ~exist('f','var'),f=figure('Name','outliers Array');f.Position=f.Position*diag([1 1 0.5 0.5]);CenterWindow(handles,f,handles.figure1);end
    t=uitable('Parent',f,'Position',[20 20 f.Position(3)-40 f.Position(4)-40],'Data',outliers);
    t.Units='normalized';
end



function err=Check_nSTD_from_file()
global DataA stdv vtype vstdx Nonecol;
        
       
        if ~isempty(find(ismember([vtype vstdx],Nonecol)))           
            mess=sprintf('Either ''Sample Type'' or ''Std Value'' column not set.\n\nCannot determine the number of STD used in the file.');
            mmm=msgbox(mess);err='1 column not set';
            return ;
        end
        err='';
        cond1=(~cellfun(@isempty,regexp(DataA(2:end,vtype),'^STD'))); %is a std 
        maxstd=max(cellfun(@str2num,regexprep(DataA(find(cond1)+1,vtype),'\D','')));
       % maxstd=max(cellfun(@str2num,strrep(strrep(DataA(find(cond1)+1,vtype),'STD',''),'-DRAIN',''))); %get max STD number
        
        if maxstd ~= length(stdv)
                     mess=sprintf(['File seems to have an Inconsistent number of STDs.\n\nMake sure the correct system configuration.ini file is loaded!'...
                         '\n\nEdit table of STDs (Add or Remove) by right-clicking the table.'...
                         '\n\nFields such as Std Offset or SST interpolated are being automatically re-calculated']);
                     mess2=sprintf(['\n\nYour options are:' '\n\n\t 1 - Let the program fix it automatically (Program fix).' '\n\n\t 2 - Do it yourself (Manual fix)']);
                     err=questdlg([mess mess2],'Data Load warning','1 - Program fix','2 - Manual fix','2 - Manual fix');
        end


        
function create_ini_file(fdn,fn)
global def_folder vartitlesc;    
        
 switch fn
     case [mfilename, '.ini']
        fid=fopen([fdn filesep fn],'wt');
        if fid<0
            mmm=msgbox(sprintf('Enable to create:\n\n%s\nProgram function not guaranteed.\n',fn),'file error','error','modal');
            uiwait(mmm);return;           
        else
            fprintf(fid,'%s\n%s\n%s\n%s\n','Default Data File Path',def_folder,'Default Result File Path',def_folder);
            fprintf(fid,'%s\n%s\n%s\n%s\n','Default Header File Path',def_folder,'Default Header File','Initial Data Config.csv');
            fprintf(fid,'%s\n%s\n%s\n%s\n','Default System Configuration File Path',def_folder,'Default System Configuration File','Initial System Config.ini');
            fprintf(fid,'%s\n%s','Default xml Data Path',def_folder);
            
            fclose(fid);
        end
        
     case 'Initial Data Config.csv'
         fid=fopen([fdn filesep fn],'wt');
         if fid<0
             mmm=msgbox(sprintf('Enable to create:\n\n%s\nIn\n%s\nProgram function not guaranteed.\n',fn,fdn),'file error','error','modal');
             uiwait(mmm);return;
         else
             l1=sprintf('%s,',vartitlesc{1:end-3});%all but QC, flag...
             l1=l1(1:end-1);%last comma out
             l2='YDay Calc,Date,Time,latitude,longitude,equ press,equ temp,H2O flow,licor flow,cell temp,cell press,Type,TSG Salinity V,TSG Temp Corr,CO2 fit,CO2 um/m,H2O mm/m,atm press,PC temp,DB temp,atm press,none';
             fprintf(fid,'%s\n%s',l1,l2);
             fclose(fid);
         end
     case 'Initial System Config.ini'
         fid=fopen([fdn filesep fn],'wt');
         if fid<0
             mmm=msgbox(sprintf('Enable to create:\n\n%s\nProgram function not guaranteed.\n',fn),'file error','error','modal');
             uiwait(mmm);return;
         else
             fprintf(fid,'Expocode\n');
             fprintf(fid,'XXXX\n');
             fprintf(fid,'Group\n');
             fprintf(fid,'Enter Group\n');
             fprintf(fid,'Ship\n');
             fprintf(fid,'Enter Ship\n');
             fprintf(fid,'Cruise ID\n');
             fprintf(fid,'cruise_ID\n');
             fprintf(fid,'Barometer Height (m)\n');
             fprintf(fid,'0\n');
             fprintf(fid,'EQU Pressure Differential(1-Yes 0-No)\n');
             fprintf(fid,'1\n');
             fprintf(fid,'Standard Values (1,2,3,..)\n');
             fprintf(fid,'0.00	200.00	300.00	400.00	500.00\n');
             fprintf(fid,'Use these Values(1-Yes 0-No)\n');
             fprintf(fid,'1\n');
             fprintf(fid,'Use Zero Std in fit(1-Yes 0-No)\n');
             fprintf(fid,'1\n');
             fprintf(fid,'min-max for gas flow, water flow, sst, sss, condsr, dt\n');
             fprintf(fid,'  40	 200	 0.3	   5	  -2	  32	 0.1	  45	   0	   5	-0.5	   1\n');
             fprintf(fid,'p min-max and t min-max for equ, licor, deck box\n');
             fprintf(fid,' 900	1040	  -2	  34	 900	1040	  10	  50	 900	1040	 -20	  50\n');
             fprintf(fid,'Abnormal values for sst and sss\n');
             fprintf(fid,'1	1\n');
             fprintf(fid,'Range Check Flags\n');
             fprintf(fid,'3	3	3	3	3	3	3	3	3	3	3	3	3	3\n');
             fprintf(fid,'Insitu Temperature Offset(minutes)\n');
             fprintf(fid,'0\n');
             fprintf(fid,'Correct xCO2 for water(1-Yes 0-No)\n');
             fprintf(fid,'0\n');
             fprintf(fid,'PI Names for this data\n');
             fprintf(fid,'Dupont,D.; Dupond,T.\n');
             fprintf(fid,'Save Options (1-Yes 0-No)(Expocode as col, Group as col, Ship as col, Save flags 4, Exclude ATMs)\n');
             fprintf(fid,'1	1	0	0	0');
             fclose(fid);
        end


 end


% --- Executes on button press in rdostdv.
function rdostdv_Callback(hObject, eventdata, handles)
global DataA vtype stdv vstdx Nonecol  
% hObject    handle to rdostdv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rdostdv

mmm=-1;stdv=[];
if ~get(handles.rdostdv,'Value') & (vstdx==Nonecol) %#ok<*AND2>
    mess=sprintf('Std Offsets cannot be calculated.\n\"%s\" is not assigned yet.\nCorrect and Re-Calc Fields','Std Value');
    stdv=[];cc='r'; set(handles.stdo_f,'BackgroundColor',cc,'TooltipString',mess);
elseif ~get(handles.rdostdv,'Value')
    if verLessThan('matlab','9.1'),stdr=find(~cellfun(@isempty,strfind(DataA(:,vtype),'STD')));
    else, stdr=find(contains(DataA(:,vtype),'STD'));end
    stds=unique(DataA(stdr,vtype));stdn=sort(unique(regexprep(stds,'\D','')));
    maxstd=max(str2num(char(stdn)));
%     cond1=(~cellfun(@isempty,regexp(DataA(2:end,vtype),'^STD\d\d?(?![sz])'))); %is a std but no zero/span
%     maxstd=max(cellfun(@str2num,strrep(strrep(DataA(find(cond1)+1,vtype),'STD',''),'-DRAIN',''))); %get max STD number
    for i=1:maxstd
        junk=find(strcmp(DataA(:,vtype),['STD' stdn{i}]),1);
        if ~isempty(junk)
            stdv=[stdv;str2num(DataA{junk,vstdx})];
            %check that all values in file agree
            stdx=find(~cellfun(@isempty,regexp(DataA(2:end,vtype),['^STD' num2str(i) '(?![sz])'])));
            if std(str2num(char(DataA(stdx+1,vstdx))))>1e-2
                strmsg=sprintf('Error in %s known values.\nStandard Deviation > 0.01\nCheck data!', ['STD' num2str(i)]);
                mmm=msgbox(strmsg,'Error','error','modal');
                CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
                cc='r';mess=strmsg; set(handles.stdo_f,'BackgroundColor',cc,'TooltipString',mess);
            end
        else %zero or span
            junk=find(strcmp(DataA(:,vtype),['STD' stdn{i} 'z']),1);
            if ~isempty(junk),stdv=[stdv;0];
            else
                stdv=[stdv;-999];
            end
        end
        if mmm~=-1,return;end
    end
    
end

if get(handles.rdostdv,'Value'),   stdv=get(handles.tblstdv,'Data'); end


% --- Executes during object creation, after setting all properties.
function list_999_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_999 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Check_for_999(handles)
global vflag vlat vlong vpatm vpeq vteq vsal vpamb vgflo vwflo vtin  vlicorx ;  
global DataA Nonecol;

hdl=[vlat vlong vpeq vsal vteq vtin vpatm vpamb vgflo vwflo vlicorx];
hdlt={'Latitude'; 'Longitude';'EQU Pressure'; 'Salinity';'EQU Temp';'in situ Temp';'ATM Pressure';'Ambient Pressure';'gas flow';'H2O flow';'licor xco2'};
junk=sprintf('%s\n','-999 found in...');

indx=find(hdl~=Nonecol);%If Data is assigned ==>
if ~isempty(indx)
    if ~isempty(find(strcmp(DataA(2:end,hdl(indx)),'-999')~=1 & strcmp(DataA(2:end,vflag),'4')~=1))
        indx2=[];
        for ix = indx
            if ~isempty(find(strcmp(DataA(2:end,hdl(ix)),'-999')==1 & strcmp(DataA(2:end,vflag),'4')~=1)) indx2=[indx2;ix]; end
        end
        if ~isempty(indx2)
            junk=sprintf('%s\n',hdlt{indx2});
            %mmm=questdlg(sprintf('You have data set to -999 which hasn''t been flagged 4 in:\n\n%s',junk),'-999 issue','STOP','CONTINUE','CONTINUE');
            %if strcmp(mmm,'STOP') return; end
            junk=sprintf('%s\n%s','-999 found in...',junk);
        end
    end
end
set(handles.list_999,'String',junk,'TooltipString',junk);



function Check_for_text(handles)
global DataA Nonecol DataAText;
global vdutc vtutc vlat vlong vpeq vteq vwflo  vgflo  vlicorcav  vpamb  vsal  vtin  vstdx  vlicorx  ...
    vlicorw  vpatm vtcpu vtdbox vpdbox vtcond;
hdl=[vlat vlong vpeq vteq vsal  vtin vpatm vpamb vgflo vwflo vlicorx vlicorw  vlicorcav    vstdx   vtcpu vtdbox vpdbox vtcond];
hdlt={'Latitude'; 'Longitude';'EQU Pressure';'EQU Temp'; 'Salinity';'in situ Temp';'ATM Pressure';'Ambient Pressure';...
    'gas flow';'H2O flow';'licor xco2';'licor xh2o';'licor cavity';'std xco2';'cpu temp';'deck box temp';'deck box pressure';...
    'condenser temp'};
textfound=0;
set(handles.list_text,'Visible','off','String','','TooltipString','');
junk=sprintf('%s\n','Text found in...');

indx2=[];DataAText=DataA(1,:);
indx=find(hdl~=Nonecol);%If Data is assigned ==>
if ~isempty(indx)
    if ~isempty(find(not(cellfun('isempty',regexp(DataA(2:end,hdl(indx)),'[^0-9\.-]')))))
        for ix = indx
            if ~isempty(find(not(cellfun('isempty',regexp(DataA(2:end,hdl(ix)),'[^0-9\.-]')))))
                indx2=[indx2;ix];
                DataAText=[DataAText;DataA(find(not(cellfun('isempty',regexp(DataA(2:end,hdl(ix)),'[^0-9\.-]'))))+1,:)];
            end
        end
        if ~isempty(indx2)
            junk=sprintf('%s%s\n',junk,hdlt{indx2});
            textfound=1;
        end
    end
    if ~isempty(find(cellfun('length',regexp(DataA(2:end,hdl(indx)),'[\.]'))>1))% finds more than one . in number
        for ix = indx
            if ~isempty(find(cellfun('length',regexp(DataA(2:end,hdl(ix)),'[\.]'))>1))
                indx2=[indx2;ix];
                DataAText=[DataAText;DataA(find(cellfun('length',regexp(DataA(2:end,hdl(ix)),'[\.]'))>1)+1,:)];
            end
        end
        if ~isempty(indx2)
            junk=sprintf('%s%s\n',junk,hdlt{indx2});
            textfound=1;
        end
    end
    
end

%now for date and time columns
hdldt=[vdutc vtutc ];
hdldtt={'Date'; 'Time'};

indx=find(hdldt~=Nonecol);%If Data is assigned ==>
if ~isempty(indx)
    if sum(ismember(['T','Z',':'],char(DataA(2:100,hdldt(indx)))))~=3
        if ~isempty(find(not(cellfun('isempty',regexp(DataA(2:end,hdldt(indx)),'[^0-9/:]')))))
            for ix = indx
                if ~isempty(find(not(cellfun('isempty',regexp(DataA(2:end,hdldt(ix)),'[^0-9/:]')))))
                    indx2=[indx2;ix];
                    DataAText=[DataAText;DataA(find(not(cellfun('isempty',regexp(DataA(2:end,hdldt(ix)),'[^0-9/:]'))))+1,:)];
                end
            end
            if ~isempty(indx2)
                junk=sprintf('%s%s\n',junk,hdldtt{indx2});
                textfound=1;
            end
        end
    end
end
if textfound
    set(handles.list_text,'Visible','on','String',junk,'TooltipString',junk);
    mmm=msgbox(junk,'Text found in data','modal');
    CenterWindow(handles,mmm,handles.figure1)
    uiwait(mmm);
end


function Display_DataText(hObject, eventdata, handles)
% hObject    handle to Data File frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DataAText ;

if strcmp(handles.figure1.SelectionType, 'alt') %right-click
    if ~exist('DataAText','var') || isempty(DataAText)
        mmm=msgbox('no text found in data array.','display Error','error','modal');
        CenterWindow(handles,mmm,handles.figure1)
        uiwait(mmm);
        return;
    end
    if ~exist('f','var'),f=figure('Name','Text in DataA Array');f.Position=f.Position*diag([1 1 0.5 0.5]);CenterWindow(handles,f,handles.figure1);end
    t=uitable('Parent',f,'Position',[20 20 f.Position(3)-40 f.Position(4)-40],'Data',DataAText);
    t.Units='normalized';
end


% --------------------------------------------------------------------
function Hide_Restore_Callback(hObject, eventdata, handles)
% hObject    handle to hrd1r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DataA ind1 ind2 hp hp2 x y Currentplot;
global sub_flag;
global vflag vsubf vsubfu vtype vlong newlong1 newlong2 vtin vtini;

% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';...
%         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';...
%         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'% sub_flag={'out of standard range';'questionable/interpolated SST';'questionable/interpolated EqT';'anomalous DT';... %         'questionable/interpolated SSS';'questionable/interpolated p';'Low LICOR Gas Flow';'questionable air value';'corrected with less than 3 standards';'other-see metadata';... %         'LICOR t';'water flow';'deck box t';'deck box p';'condenser t';'interpolated position';'interpolated using questionable SST'};  variable missing';'interpolated using questionable SST'};
    mmm=msgbox(sprintf('%s\n','This function is not ready yet'),' Issue','modal');CenterWindow(handles,mmm,handles.figure1) ; uiwait(mmm);
    return;



activity(1,handles); %drawnow;
equ_atm={'EQU';'ATM';'EQU-DRAIN';'ATM-DRAIN'};
equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};


gs=get(handles.popGselect,'String');
if strcmp(gs(1,:),'Graph 1')==1 & strcmp(gs(2,:),'Graph 2')==1
    mmm=msgbox('No data plotted','Point flagging','error','modal');
    CenterWindow(handles,mmm,handles.figure1) ;  uiwait(mmm);
    activity(0,handles);return;
end


ax_fl=get(h,'Tag');  %Tag is 'fl'+axe # + flag #
axstr=ax_fl(3);fl=ax_fl(4);
subflag=ax_fl(5:end);
if strcmp(axstr,'1')==1
    ax=handles.axes1;
    yax=handles.LeftY;
    ind=ind1;
    chp=hp;
    if exist('newlong1','var') newlong=newlong1; end
else
    ax=handles.axes2;
    yax=handles.RightY;
    ind=ind2;
    chp=hp2;
    if exist('newlong2','var') newlong=newlong2; end
end

xlimits = get(ax,'XLim');
ylimits = get(ax,'YLim');
xlimits1 = get(handles.axes1,'XLim');
ylimits1 = get(handles.axes1,'YLim');
xlimits2 = get(handles.axes2,'XLim');
ylimits2 = get(handles.axes2,'YLim');
xp = (xlimits(2)-xlimits(1))/100;
yp = (ylimits(2)-ylimits(1))/100;
type='';types='';
switch subflag
    case {'3','4','5','6','7'} % EqT, DT, SSS, pressures, low EQ flow
        type=equ;types='EQU';
    case {'0','1','2','9','10'}% no subflag, std range, SST, <3 std corr, Other
        type=equ_atm;types='EQU and ATM';
    case {'8'}% Q? Air value
        type=atm;types='ATM';
end

if get(yax,'Value')==vtin && str2num(fl)==4
    msg1='MEASURED SST cannot be flagged 4 due to the time offset interpolation.';
    msg2='Only INTERPOLATED SST (Tin interp) can be flagged 4.';
    msg3='Either INTERPOLATE or GENERATE from other source and flag 3 for Quest. SST';
    msg4='Or simply flag 3 for Quest. SST';
    mmm=msgbox(sprintf('%s\n%s\n%s\n%s',msg1,msg2,msg3,msg4),'Flag Issue','modal');CenterWindow(handles,mmm,handles.figure1) ; uiwait(mmm);
    return;
end
if get(yax,'Value')==vtini && str2num(fl)==3
    msg1='INTERPOLATED SST cannot be flagged 3 manually.';
    msg2='If MEASURED SST used to interpolate it is not questionable...';
    msg3='then simply flag 4.';
    mmm=msgbox(sprintf('%s\n%s\n%s',msg1,msg2,msg3),'Flag Issue','modal');CenterWindow(handles,mmm,handles.figure1) ; uiwait(mmm);
    return;
end

if ax_fl(1)~='r', return;end   %point only, not rectangle

    
%rectangle

rh=findobj(ax,'Type','Rectangle');
ptaallr=[];rm=0;%all rectangles
if isempty(rh)
    mmm=msgbox('No area selected','Area Flagging','error','modal');
    uiwait(mmm);activity(0,handles);return;
end

for r=1:length(rh)% for each rectangle
    rv=get(rh(r),'Position');
    xv=[rv(1) rv(1) rv(1)+rv(3) rv(1)+rv(3)];
    yv=[rv(2) rv(2)+rv(4) rv(2)+rv(4) rv(2)];
    pta1all=[];pta2all=[];ptaall=[];
    for i =1:length(chp)
        if (handles.CurrentX.Value==vlong & get(handles.rdodateline,'value'))
            pta1{i} = inpolygon(newlong{i},str2num(char(DataA(ind{i},get(yax,'Value')))),xv,yv);
        else
            pta1{i} = inpolygon(str2num(char(DataA(ind{i},handles.CurrentX.Value))),str2num(char(DataA(ind{i},get(yax,'Value')))),xv,yv);
        end
        %pta1 are points in polygon, pta2 are points of right type
        pta1{i}=find(pta1{i});pta1all=union(pta1all,ind{i}(pta1{i}));
        pta2{i}=[];
        if str2num(fl)==3 & ~strcmp(type,''),pta2{i}=find(ismember(DataA(ind{i},vtype),type));
        else pta2{i}=pta1{i};end
        pta2all=union(pta2all,ind{i}(pta2{i}));
        pta{i}=intersect(pta1{i},pta2{i});
    end
    % ptaall are points of right type in rectangle
    ptaall=intersect(pta1all,pta2all);
    % ptaallr are points of right type in all rectangles drawn
    ptaallr=union(ptaallr,ptaall);
    
    %for fl=3 if (all or some) points are of the wrong type (pta1all>ptaall) --> message
    %for fl=4, can flag anything
    if str2num(fl)==3 & ~isempty(pta1all) & (size(pta1all,1)>size(ptaall,1)), rm=1;end
end %end of each rectangle

ptaallr=unique(ptaallr); %to not select same points several times in case rectangles intersect each other
%if points of the wrong type were selected: message
if rm==1, mmm=msgbox(sprintf('Only %s will be flagged for this parameter.',types),'Flag Issue','modal');CenterWindow(handles,mmm,handles.figure1) ; uiwait(mmm);end

%     for i = 1:length(chp)%for each plot, search if flag present, and adds if not
%
%         if str2num(fl)~=3 %flags 2 and 4
%             DataA(ind{i}(pta{i}),vflag)={fl};
%             DataA(ind{i}(pta{i}),vsubf)=cellstr('');
%             DataA(ind{i}(pta{i}),vsubfu)=cellstr('');
%         else
%             if str2num(subflag)==0, Set_subflag(ind{i}(pta{i}),str2num('-1'),handles);  % Set_subflag sets flag to 4 if subflag=0
%             elseif str2num(subflag)<=10,Set_subflag(ind{i}(pta{i}),str2num(subflag),handles);
%             elseif str2num(subflag)>10, Set_subflag(ind{i}(pta{i}),10,handles);Set_subflaguser(ind{i}(pta{i}),str2num(subflag),handles);
%             end
%         end
%     end

if str2num(fl)~=3 %flags 2 and 4
    DataA(ptaallr,vflag)={fl};
    DataA(ptaallr,vsubf)=cellstr('');
    DataA(ptaallr,vsubfu)=cellstr('');
else
    if str2num(subflag)==0, Set_subflag(ptaallr,str2num('-1'),handles);  % Set_subflag sets flag to 4 if subflag=0
    elseif str2num(subflag)<=10,Set_subflag(ptaallr,str2num(subflag),handles);
    elseif str2num(subflag)>10, Set_subflag(ptaallr,10,handles);Set_subflaguser(ptaallr,str2num(subflag),handles);
    end
end

feval(Currentplot,1,1,handles,1);
set(handles.axes1,'XLim',xlimits1);set(handles.axes1,'YLim',ylimits1);
set(handles.axes2,'XLim',xlimits2);set(handles.axes2,'YLim',ylimits2);



clear pti1all pti2all ptiall pti1 pti2 pti;
clear pta1all pta2all ptaall pta1 pta2 pta ptaallr;

Stat_Update_Callback(0,0,handles);
Check_for_999(handles);
activity(0,handles);
