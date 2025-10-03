function varargout = Edit_xml(varargin)

% EDIT_XML MATLAB code for Edit_xml.fig
%      EDIT_XML, by itself, creates a new EDIT_XML or raises the existing
%      singleton*.
%
%      H = EDIT_XML returns the handle to a new EDIT_XML or the handle to
%      the existing singleton*.
%
%      EDIT_XML('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDIT_XML.M with the given input arguments.
%
%      EDIT_XML('Property','Value',...) creates a new EDIT_XML or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Edit_xml_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Edit_xml_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Edit_xml

% Last Modified by GUIDE v2.5 24-Apr-2015 13:35:59


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Edit_xml_OpeningFcn, ...
                   'gui_OutputFcn',  @Edit_xml_OutputFcn, ...
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




% --- Executes just before Edit_xml is made visible.
function Edit_xml_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Edit_xml (see VARARGIN)
global il maxtags; il=0;maxtags=15; %maxtags is maximum number of similar tags (like investigator, ports_of_call...)

% Choose default command line output for Edit_xml
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Edit_xml wait for user response (see UIRESUME)
% uiwait(handles.figure1);
initialize_gui(handles);
% Use system color scheme for figure:
set(hObject,'Color',get(0,'DefaultUicontrolBackgroundColor'));

movegui(hObject,'center');




% --- Outputs from this function are returned to the command line.
function varargout = Edit_xml_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in xml_list.
function xml_list_Callback(hObject, eventdata, handles)
% hObject    handle to xml_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xml_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xml_list
global  xml  xml_list fieldn pocs pt ppl rg crt instr maxtags reds spch3;
%fieldn(i) contains the complete field name of the structure corresponding to
%the value selected on row i. It's created in list_struct function
%in regexp: .*? means any number of characters    if enclosed in () then
%becomes a token (/\1 means "/token#1" which is first tag)

if isempty(xml_list),msgbox('No data listed...Import an xml file.','modal');uiwait;return;end

%determines closest 'red line'
lvalue=get(handles.xml_list,'Value');
[m,linen]=min(abs(reds-lvalue));
set(handles.redslider,'Value',linen);

line_s=get(hObject,'Value');if isempty(xml_list{line_s}),msgbox('Select a line with a value','modal');uiwait;return;end

if ~isempty(regexp(xml_list{line_s},'black.*?(?=[blured]+)'))
    %Selected from a list saved in xml.tsv.txt
    pocnum=findstr(char(fieldn(line_s)),'xml.x_tags.Cruise_Info.Experiment.Cruise.Ports_of_Call');
    isppl=regexpi(char(fieldn(line_s)),'(.User.|.Investigator)','match');
    isaddinfo=regexpi(char(fieldn(line_s)),'Additional_Information','match');
    isstdinfo=regexpi(char(fieldn(line_s)),'Manufacturer_of_Calibration_Gas','match');
    
    isplattype=regexp(char(fieldn(line_s)),'Platform_Type','match');
    iscrtype=regexp(char(fieldn(line_s)),'Experiment_Type','match');
    isinstr=regexp(char(fieldn(line_s)),'Co2_Instrument_Type','match');
    
    if ~isempty(pocnum)
        pocnum=str2num(fieldn{line_s}(length('xml.x_tags.Cruise_Info.Experiment.Cruise.Ports_of_Call')+2:end));
        [matc,tok]=regexp(xml_list{line_s},'<html>.*?<font color="[blured]+"><i>(.*?)</i></font>','match','tokens');
        newpocs={};newpocs=strtrim(inputdlg('Edit the list below...','Ports of Call',10,tok{:},'on'));
        del=[];if isempty(newpocs),return;end
        for i=1:size(newpocs{1},1),    if isempty(strtrim(newpocs{1}(i,:))),del=[del i];end,   end
        newpocs{1}(del,:)=[];
        if size(newpocs{1},1)==0,       pocs(pocnum)=[];
        elseif size(newpocs{1},1)==1,   pocs{pocnum}=newpocs{1};
        elseif size(newpocs{1},1)>1,    pocs{pocnum}=newpocs{1}(1,:);   for i=2:size(newpocs{1},1),pocs{end+1}=newpocs{1}(i,:);end
        end
        %Ports of Call
        xml=Update_Mult_Sim_tags(xml,'x_tags.Cruise_Info.Experiment.Cruise','Ports_of_Call',pocs,maxtags);
        
    elseif ~isempty(isppl)%User, Investigator
        if strcmp(char(isppl),'.Investigator')==1
            %invnum=str2num(fieldn{line_s}(length('xml.Investigator')+2:end));
            invnum=regexp(fieldn{line_s},'.*?_([0-9]*)[.].*?','tokens');
            invnum=str2num(char(invnum{1}));
            %[matc,tok]=regexp(xml_list{line_s},'<html>.*?<font color="blue"><i>(.*?)</i></font>','match','tokens');
            [s,ok] = listdlg('PromptString','Select one or more Investigators:', 'ListString',char(ppl.Name));
            if ~ok, return;end
            xml.x_tags=Set_Investigators(xml.x_tags,s,maxtags);
        else %User
            [s,ok] = listdlg('PromptString','Select a person:','SelectionMode','single', 'ListString',char(ppl.Name));
            if ~ok, return;end
            f=fieldnames(xml.x_tags.(isppl{1}(2:end-1)));
            for i=1:numel(f),xml.x_tags.(isppl{1}(2:end-1)).(f{i})=ppl(s).(f{i});end
        end
        
    elseif ~isempty(isplattype)%Platform Type
        [s,ok] = listdlg('PromptString','Select a platform type:','SelectionMode','single', 'ListString',char(pt));
        if ~ok, return;end
        eval([char(fieldn(line_s)) '=char(pt{s});']);
               
    elseif ~isempty(iscrtype)%Cruise Type
        [s,ok] = listdlg('PromptString','Select a cruise type:','SelectionMode','single', 'ListString',char(crt));
        if ~ok, return;end
        eval([char(fieldn(line_s)) '=char(crt{s});']);
        
    elseif ~isempty(isinstr)%CO2 Instr. Type
        [s,ok] = listdlg('PromptString','Select a CO2 instrument type:','SelectionMode','single', 'ListString',char(instr));
        if ~ok, return;end
        eval([char(fieldn(line_s)) '=char(instr{s});']);
        
    elseif ~isempty(isaddinfo)%Additional Info
        comments=xml.x_tags.Additional_Information;
        newinfo=strtrim(inputdlg('Edit the text below...','Additional Info',10,{comments},'on'));
        if ~isempty(newinfo), xml.x_tags.Additional_Information=char(newinfo);end
        
    elseif ~isempty(isstdinfo)%STD Info
        comments=xml.x_tags.Method_Description.CO2_Sensors.CO2_Sensor.Manufacturer_of_Calibration_Gas;
        newinfo=strtrim(inputdlg('Edit the text below...','Standards Info',10,{comments},'on'));
        if ~isempty(newinfo), xml.x_tags.Method_Description.CO2_Sensors.CO2_Sensor.Manufacturer_of_Calibration_Gas=char(newinfo);end
        
    else %normal string edit
        [matc,tok]=regexp(xml_list{line_s},'<html>.*? <.*?>(\w+)<.*?>: <.*?><i>(.*?)</i></font><.*?> /\1<.*?>','match','tokens');
        newstr=strtrim(inputdlg('Edit value...','Field Value',[1, length(char(tok{:}(2)))+10],tok{:}(2)));
        if isempty(newstr),return;end
        if ~isempty(newstr) & ~strcmp(newstr,tok{:}(2))
            newstr=regexprep(newstr,'''','''''');%replaces apostrophy by ''
            eval( [char(fieldn(line_s)) '='''  newstr{1} ''';']);
        end
    end
    
elseif ~isempty(regexp(xml_list{line_s},'[blured]+')) %items on more than 1 line (stds, addl info,...)
    %might have to add &#xD (newline char)
    switch char(fieldn(line_s))
        case 'xml.x_tags.Method_Description.CO2_Sensors.CO2_Sensor.Manufacturer_of_Calibration_Gas'
            comments=xml.x_tags.Method_Description.CO2_Sensors.CO2_Sensor.Manufacturer_of_Calibration_Gas;
            newinfo=strtrim(inputdlg('Edit the text below...','Standards Info',10,{comments},'on'));
            if ~isempty(newinfo), xml.x_tags.Method_Description.CO2_Sensors.CO2_Sensor.Manufacturer_of_Calibration_Gas=char(newinfo);end
            
        case 'xml.x_tags.Additional_Information'
            comments=xml.x_tags.Additional_Information;
            newinfo=strtrim(inputdlg('Edit the text below...','Additional Info',10,{comments},'on'));
            if ~isempty(newinfo), xml.x_tags.Additional_Information=char(newinfo);end
            
        case 'xml.x_tags.Cruise_Info.Experiment.Cruise.Geographical_Coverage.Geographical_Region'
            [s,ok] = listdlg('PromptString','Select one or more regions:', 'ListString',char(rg));
            if ~ok, return;end
            xml.x_tags.Cruise_Info.Experiment.Cruise.Geographical_Coverage.Geographical_Region=char(rg{s});
    end
else
    [matc,tok]=regexp(xml_list{line_s},'<html>.*?<font color="black">(.*?)</font>','match','tokens');
    msgbox(sprintf('You selected: %s\n\nSelect a line with a value',char(tok{:})),'modal');uiwait;return;
end
%end
    
lvl=1;list_struct(handles,xml);

xml_list{numel(xml_list)+1,1}='';
set(handles.xml_list,'String',xml_list);
set(handles.xml_list,'Value',lvalue);





% --- Executes during object creation, after setting all properties.
function xml_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xml_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --------------------------------------------------------------------
function initialize_gui(handles)
global xml_list  vcruiseid  pt ppl rg crt ef instr projDir  spch spch2 il ffn;
global datap resultp headp hfile sysinip sysinif where_xml;


vcruiseid=strtrim(vcruiseid);resultp=strtrim(resultp);
rg={};crt={};ef={};
btns=[handles.specialc;handles.specialc_txt;handles.tags_present;handles.tags_present_txt];
%replace special characters: see 'extract_lvlTag_from_text'
spch={'&' '&amp;';'<' '&lt;';'>' '&gt;';'oC' '&#176;C';'+-' '&#177;';'µ' '&#956;'};
spch2={'\[' '*b1*';'\]' '*b2*';'\{' '*c1*';'\}' '*c2*';'\(' '*p1*';'\)' '*p2*';'\=' '*eq*';'''' '*ap*';...
    '\,' '*cm*';'\;' '*sc*';'\:' '*cl*';'\%' '*pc*';'\!' '*ep*';'\@' '*at*';'&' '*and*';};
ppl=struct('Name','','Organization','','Address','','Phone','','Email',''); fppl=fieldnames(ppl);

%reads 'final.csv' file name to input as default for readme.
offn='';
if ~isempty(ffn),[junk,offn,junk]=fileparts(ffn);end
offn=regexprep(offn,'_Final','');
if ~isempty(offn),set(handles.xml_fn,'String',offn,'ToolTipString',offn);
elseif ~isempty(vcruiseid),set(handles.xml_fn,'String',vcruiseid,'ToolTipString',vcruiseid);end
if ~isempty(resultp),set(handles.xml_path,'String',resultp,'ToolTipString',resultp);end
set(btns,'visible','off');


if isempty(resultp),set(handles.xml_path,'String',projDir);end

%set highlight color of listbox
v  = version ('-release'); 
vn = str2num(v(1:4));

if (isempty(il) | ~il) & (vn>2010) %set highlight color of listbox
    jScrollPane  = findjobj(handles.xml_list); % get the scroll-pane object
    jListbox = jScrollPane.getViewport.getComponent(0);
    if vn<2014 
        set(jListbox, 'SelectionBackground',[220/255 220/255 220/255]); % option #2
        set(jListbox, 'SelectionForeground',[0 0 0]); % option #2
    else
        set(jListbox, 'SelectionBackground',java.awt.Color(220/255,220/255,220/255)); % option #2
        set(jListbox, 'SelectionForeground',java.awt.Color(0, 0, 0)); % option #2
    end
        il=1;
end

%Populate info lists
mots={'Platform Type Info';'People Info';'Cruise Type Info';'Edit_Fields Info';'CO2 Instr Info'};
exprs={'(.*)';'(.*)\t(.*)\t(.*)\t(.*)\t(.*)';'(.*)';'(\S*)';'(.*)'};tok={};mot=0;


fid=fopen([where_xml filesep 'xml.tsv.txt'],'r');lines={};
if fid<0
    msgbox(sprintf('Problem opening ''xml.tsv.txt'' file\nin %s.',where_xml));
else
    if ~isempty(where_xml),set(handles.txt_xml_path,'String',sprintf('Path of xml.tsv.txt file:\n%s',where_xml),'ToolTipString',sprintf('Path of xml.tsv.txt file:\n%s',where_xml));end

    while ~feof(fid)
        lines=textscan(fid,'%s',1,'delimiter','\n');lines{1}=regexprep(lines{1},'"','');
        if isempty(char(lines{1})),continue;end
        [matc,tok]=regexp(lines{1},'(.*? Info).*?','match','tokens'); %#ok<*ASGLU>
        if ~isempty(tok{1,1}),mot=find(ismember(mots,tok{1}{1}));if isempty(mot),mot=0;end,c=0;continue;end
        if mot>0 & ~isempty(strtrim(lines{1})),[matc,tok]=regexp(strtrim(lines{1}),exprs{mot},'match','tokens');c=c+1;end
        switch mot
            case 1 %Populates list for Platform Type
                pt{c}=tok{1}{1}{1};
            case 2 %Populates list for People
                for i= 1:numel(fppl), ppl(c).(fppl{i})=tok{1}{1}{i};  end
            case 3 %Cruise Type
                crt{c}=tok{1}{1}{1};
            case 4 %Fields to Edit (display in red)
                ef{c}=['xml.' tok{1}{1}{1}];
                for i=1:size(spch2,1),ef{c}=regexprep(ef{c},spch2{i,1},spch2{i,2});end
            case 5 %CO2 Instr Type
                instr{c}=tok{1}{1}{1};
        end
    end
    
    fclose(fid);
end
ppl(end+1).(fppl{1})='Reset Field';

set(handles.btnPop,'Visible','off');
if numel(xml_list)
    set(handles.xml_list,'String', xml_list);set(handles.btnPop,'Visible','on');

    %stores the position of red lines to scroll through them
    reds=~cellfun(@isempty,regexp(xml_list,'color="red"')); reds=find(reds);reds=flip(reds,1);
    set(handles.redslider,'Min',0);if numel(reds)>0, set(handles.redslider,'Max',numel(reds)); else set(handles.redslider,'Max',1);end
    if numel(reds)>0
        set(handles.redslider,'SliderStep',[1/numel(reds) 1]);%1st number is step when hitting arrow, second when hitting bar
        %Slider is inverted: max value = line 1 and 1 = line max
        %determines closest 'red line'
        lvalue=get(handles.xml_list,'Value');
        [m,linen]=min(abs(reds-lvalue));
        set(handles.redslider,'Value',linen);
    end
end
guidata(handles.figure1, handles);





function list_struct(handles,s)
global lvl line xml_txt xml_list fieldni fieldnic fieldn ef reds spch;

% recursive function
% if it's a structure, print fieldname on next line then expand next level.
% if not, print "field(black): value(blue or red)"
% fieldn is array containing structure names, fieldni is the complete
% fieldname at that time

txtspc=2;listspc=6;skip=1;
indentlvl=repmat(sprintf(' '),1,txtspc*(lvl-1));
indentlvlp1=repmat(sprintf(' '),1,txtspc*(lvl+1-1));
lindentlvlm1=repmat(sprintf('&nbsp '),1,listspc*(lvl-1));
lindentlvl=repmat(sprintf('&nbsp '),1,listspc*(lvl));

    if lvl==1, line=-1;xml_txt='';xml_list={};fieldn={};fieldni='xml';fieldnic=fieldni; end
    
    if ~isstruct(s)
        set([handles.tags_present;handles.tags_present_txt],'visible','on');
        return;
    end
    f=fieldnames(s);
    for i=1:numel(f)
        colorn='blue'; 
        %gets rid of '_number' if present until _999
        nf=regexprep(char(f{i}),'_[0-9][0-9]?[0-9]?$','');
        
        if i>1, fieldni=regexprep(fieldni,'.\w+$','');fieldnic=regexprep(fieldnic,'.\w+$','');end
        
        fieldni=[fieldni '.' f{i}]; fieldnic=[fieldnic '.' nf]; %fieldnic has no '_number'
        xml_txt=sprintf('%s%s%s',xml_txt,indentlvl,['<' nf '>']);
        line=line+1+skip;xml_list{line,1}=sprintf('<html>%s<font color="black">%s</font>',lindentlvlm1,[nf]);

        if isstruct(s.(f{i}))
            lvl=lvl+1; xml_txt=sprintf('%s\r\n',xml_txt);list_struct(handles,s.(f{i}));lvl=lvl-1;
            fieldni=regexprep(fieldni,'.\w+$','');fieldnic=regexprep(fieldnic,'.\w+$','');
            xml_txt=sprintf('%s%s%s\r\n',xml_txt,indentlvl,['</' nf '>']);
            line=line+1+skip;xml_list{line,1}=sprintf('<html>%s<font color="black">%s</font>',lindentlvlm1,['/' nf]);
        else  %display values of field
            %fieldnic
            if sum(strcmpi(fieldnic,ef))>0,colorn='red';end
            if size(s.(f{i}),1)<=1 % will be 0 if field has no value
                %nt=regexprep(char(s.(f{i})),'[\<]','&lt;');   nt=regexprep(nt,'[\>]','&gt;');
                nt=char(s.(f{i}));
                % spch={'&' '&amp;';'<' '&lt;';'>' '&gt;';'oC' '&#176;C';'+-' '&#177;';'µ' '&#956;'};
                for j=1:size(spch,1),nt=regexprep(nt,spch{j,1},spch{j,2});end
                %xml_txt=sprintf('%s%s%s\r\n',xml_txt,s.(f{i}),['</' nf '>']);
                xml_txt=sprintf('%s%s%s\r\n',xml_txt,nt,['</' nf '>']);
                %xml_list{line,1}=sprintf('%s: <font color="%s"><i>%s</i></font><font color="black"> %s</font>',xml_list{line},colorn,nt,['/' nf]);
                xml_list{line,1}=sprintf('%s: <font color="%s"><i>%s</i></font><font color="black"> %s</font>',xml_list{line},colorn,char(s.(f{i})),['/' nf]);
                fieldn{line}=fieldni;
            elseif size(s.(f{i}),1)>1
                xml_list{line,1}=sprintf('%s:',xml_list{line});
                for j=1:size(s.(f{i}),1)
                    %nt=regexprep(char(s.(f{i})(j,:)),'[\<]','&lt;');   nt=regexprep(nt,'[\>]','&gt;');
                    nt=char(s.(f{i})(j,:));
                    % spch={'&' '&amp;';'<' '&lt;';'>' '&gt;';'oC' '&#176;C';'+-' '&#177;';'µ' '&#956;'};
                    for jj=1:size(spch,1),nt=regexprep(nt,spch{jj,1},spch{jj,2});end
                    %xml_txt=sprintf('%s\r\n%s%s',xml_txt,indentlvlp1,s.(f{i})(j,:));
                    xml_txt=sprintf('%s\r\n%s%s',xml_txt,indentlvlp1,nt);
                    %line=line+1+skip;xml_list{line,1}=sprintf('<html>%s<font color="%s"><i>%s</i></font>',lindentlvl,colorn,nt);
                    line=line+1+skip;xml_list{line,1}=sprintf('<html>%s<font color="%s"><i>%s</i></font>',lindentlvl,colorn,char(s.(f{i})(j,:)));
                    fieldn{line}=fieldni;
                end
                xml_txt=sprintf('%s\r\n%s%s\r\n',xml_txt,indentlvl,['</' nf '>']);
                line=line+1+skip;xml_list{line,1}=sprintf('<html>%s<font color="black">%s</font>',lindentlvlm1,['/' nf]);
            end
        end
    end
    
    if lvl==1
        %stores the position of red lines to scroll through them
        reds=~cellfun(@isempty,regexp(xml_list,'color="red"')); reds=find(reds);reds=flipdim(reds,1);
        set(handles.redslider,'Min',0);if numel(reds)>0, set(handles.redslider,'Max',numel(reds)); else set(handles.redslider,'Max',1);end
        if numel(reds)>0
            set(handles.redslider,'SliderStep',[1/numel(reds) 1]);%1st number is step when hitting arrow, second when hitting bar
            %Slider is inverted: max value = line 1 and 1 = line max
            %determines closest 'red line'
            lvalue=get(handles.xml_list,'Value');
            [m,linen]=min(abs(reds-lvalue));
            set(handles.redslider,'Value',linen);
        end
    end 
%     set(handles.redslider,'Value',numel(reds));
%     redslider_Callback(handles.redslider, -1, handles);

    


function s=pop_xml_from_file(s)
global vship vcruiseid expot expod pocs ppl rg crt maxtags;
global DataA vlat vlong vflag vyday vdutc vtype stdv where_xml;


if isempty(DataA), msgbox('No Data to Import.');return;end
if ~isstruct(s), msgbox('import xml file first.');return;end
if ~isfield(s,'x_tags')
    mmm=questdlg(sprintf('Loaded xml doesn''t seem to have the proper fields.\nContinue Anyway?'),'Populating File Error','Continue','Stop','Continue');
    if strcmp(mmm,'Stop')>0,return;end
end
fid=fopen([where_xml filesep 'xml.tsv.txt'],'r');lines='';

mots={'Standards Info';'Platform Info'};
exprs={'(.*)\t(.*)\t(.*)';'(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)'};tok={};mot=0;

stdstr={};
for i=1:length(stdv), stdstr{i}='Unknown';end
%Determines time between each standard (median of all time diff between
%each STD
non_zero_stds=numel(find(stdv~=0));

%if STDxz exists, finds which is used and if not listed in stds, adds line for std x+1
if sum(~cellfun(@isempty,(regexpi(unique(DataA(:,vtype)),'STD.z'))))
    
    zeros= regexpi(unique(DataA(:,vtype)),'STD(.)z','tokens');
    zeroi= find(~cellfun(@isempty,zeros));zeroi=zeros{zeroi(1)};
    if numel(find(stdv==0))==0 %if no std has 0 concentration
       
      % determines frequency of measurement of that STD in hours
       stdii=find(~cellfun(@isempty,regexpi(DataA(:,vtype),['STD' char(zeroi{1})])));%indices of that STD
       dt=str2num(char(DataA(stdii(3:end),vyday)))-str2num(char(DataA(stdii(2:end-1),vyday)));%time between each
       dt2=dt(dt>0.5/24);%don't count measurement within 1/2 hour of each other.
       dti=round(median(dt2)*24*2,0)/2;% *2 and /2 to round to nearest 0.5

       stdstr{i+1}=regexprep(sprintf('Std %d: %s, %.2f ppm, owned by %s, used every ~%.1f hours.',str2num(char(zeroi{1})),'', 0, 'AOML', dti),' , ',' ');
    end
end
while ~feof(fid)
    lines=textscan(fid,'%s',1,'delimiter','\n');lines{1}=regexprep(lines{1},'"','');
    if isempty(char(lines{1})),continue;end
    [matc,tok]=regexp(lines{1},'(.*? Info).*?','match','tokens'); %#ok<*ASGLU>
    if ~isempty(tok{1,1}),mot=find(ismember(mots,tok{1}{1}));if isempty(mot),mot=0;end,c=0;continue;end
    if mot>0 & ~isempty(strtrim(lines{1})),[matc,tok]=regexp(strtrim(lines{1}),exprs{mot},'match','tokens');c=c+1;end
    switch mot
        case 1 %Stds
            stdi=find(stdv==str2num(tok{1}{1}{1}));
      % determines frequency of measurement of that STD in hours
            stdii=find(~cellfun(@isempty,regexpi(DataA(:,vtype),['STD' num2str(stdi)])));%indices of that STD
            dt=str2num(char(DataA(stdii(3:end),vyday)))-str2num(char(DataA(stdii(2:end-1),vyday)));%time between each
            dt2=dt(dt>0.5/24);%don't count measurement within 1/2 hour of each other.
            dti=round(median(dt2)*24*2,0)/2;% *2 and /2 to round to nearest 0.5
            if ~isempty(stdi),stdstr{stdi}=regexprep(sprintf('Std %d: %s, %.2f ppm, owned by %s, used every ~%.1f hours.',stdi,char(tok{1}{1}{2}), stdv(stdi), char(tok{1}{1}{3}), dti),' , ',' ');end
        case 2 %Platform
            if strcmp(expot,strtrim(tok{1}{1}{1}))
                s.x_tags.Cruise_Info.Vessel.Vessel_ID=expot;
                s.x_tags.Cruise_Info.Vessel.Vessel_Name=tok{1}{1}{2};
                s.x_tags.Cruise_Info.Experiment.Platform_Type=tok{1}{1}{3};
                s.x_tags.Cruise_Info.Vessel.Vessel_Owner=tok{1}{1}{4};
                s.x_tags.Cruise_Info.Experiment.Experiment_Type=tok{1}{1}{5};
                s.x_tags.Cruise_Info.Experiment.Cruise.Cruise_Info=tok{1}{1}{6};
            end
    end
end

fclose(fid);

if  ismember('Unknown',stdstr), msgbox('Error Searching for Standards Info in   xml.tsv.txt   file.');end
stdstr=char(stdstr);
s.x_tags.Method_Description.CO2_Sensors.CO2_Sensor.Manufacturer_of_Calibration_Gas=stdstr;
s.x_tags.Method_Description.CO2_Sensors.CO2_Sensor.No_Of_Non_Zero_Gas_Stds=strtrim(num2str(non_zero_stds,'%d'));

% %Ports of Call
% s=Update_Mult_Sim_tags(s,'x_tags.Cruise_Info.Experiment.Cruise','Ports_of_Call',pocs,maxtags);

equ={'EQU';'EQU-DRAIN'};
atm={'ATM';'ATM-DRAIN'};
equ_atm=cat(1,equ,atm);
    
indf1=find(ismember(DataA(:,vtype),equ_atm));
indf2=find(strcmp(DataA(2:end,vlat),'-999')~=1 & strcmp(DataA(2:end,vlong),'-999')~=1 & strcmp(DataA(2:end,vyday),'-999')~=1  & strcmp(DataA(2:end,vflag),'4')~=1)+1;
indf=intersect(indf1,indf2);

fs=char(filesep);
if isempty(expot), expot='00XX';end
if isempty(expod), expod=datestr(now,'yyyymmdd');end
if isempty(vcruiseid), vcruiseid='vcruiseid' ;end
%if isempty(vship), vship='vcruiseid' ;end

s.x_tags.Cruise_Info.Experiment.Experiment_Name=vcruiseid;
s.x_tags.Cruise_Info.Experiment.Cruise.Cruise_ID=[expot expod];

A=str2num(char(DataA(indf,vlat)));
s.x_tags.Cruise_Info.Experiment.Cruise.Geographical_Coverage.Bounds.Northernmost_Latitude=num2str(max(A)+.05,'%.1f');
s.x_tags.Cruise_Info.Experiment.Cruise.Geographical_Coverage.Bounds.Southernmost_Latitude=num2str(min(A)-.05,'%.1f');
A=str2num(char(DataA(indf,vlong)));
s.x_tags.Cruise_Info.Experiment.Cruise.Geographical_Coverage.Bounds.Westernmost_Longitude=num2str(min(A)-.05,'%.1f');
s.x_tags.Cruise_Info.Experiment.Cruise.Geographical_Coverage.Bounds.Easternmost_Longitude=num2str(max(A)+.05,'%.1f');

% A=str2num(char(DataA(indf,vyday)));
% ind=find(A==min(A),1);[ysave, msave, dsave, junk, junk1, junk2]=datevec(DataA(indf(ind),vdutc),'mm/dd/yy');
% s.x_tags.Cruise_Info.Experiment.Cruise.Temporal_Coverage.Start_Date=[num2str(ysave,'%04i') num2str(msave,'%02i') num2str(dsave,'%02i')];
% ind=find(A==max(A),1);[ysave, msave, dsave, junk, junk1, junk2]=datevec(DataA(indf(ind),vdutc),'mm/dd/yy');
% s.x_tags.Cruise_Info.Experiment.Cruise.Temporal_Coverage.End_Date=[num2str(ysave,'%04i') num2str(msave,'%02i') num2str(dsave,'%02i')];

A=datenum(DataA(indf,vdutc));
ind=find(A==min(A),1);
s.x_tags.Cruise_Info.Experiment.Cruise.Temporal_Coverage.Start_Date=datestr(DataA(indf(ind),vdutc),'yyyymmdd');
ind=find(A==max(A),1);
s.x_tags.Cruise_Info.Experiment.Cruise.Temporal_Coverage.End_Date=datestr(DataA(indf(ind),vdutc),'yyyymmdd');



function s=Update_Mult_Sim_tags(s,pf,struct,AValues,maxv)

%if ~(isfield(s, pf)),return;end   doesn't detect the number

for i=1:maxv
    if i<=size(AValues,2), eval(['s.' pf '.' struct '_' num2str(i) '=''' AValues{i} ''';']);
    else
        if (eval(['isfield(s.' pf ',''' struct '_' num2str(i) ''');']))
            eval(['s.' pf '=rmfield(s.' pf ',' '''' struct '_' num2str(i) ''');']);
        end
    end
end


function s=Set_Investigators(s,AValues,maxv)
%full fieldname = struct.pf.f
global ppl;
dflt={'Name';'Organization';'Address';'Phone';'email'};

%if ~isempty(pf), structstr=['struct.' pf ]; else structstr=['struct'];end
f=fieldnames(s.Investigator_1);
    for i=1:maxv        
        if ismember(size(ppl,2),AValues)
            if i==1, for j=1:size(f), s.Investigator_1.(f{j})=char(dflt{j});end
            else
                if (eval(['isfield(s,''Investigator_' num2str(i) ''');']))
                eval(['s=rmfield(s, ''Investigator_' num2str(i) ''');']);
                end
            end
        elseif i<=size(AValues,2)
            for j=1:size(f), eval(['s.Investigator_' num2str(i) '.(f{j})=ppl(' num2str(AValues(i)) ').(f{j});']);end
        else %deletes Investigator_x when x > number of Investigators chosen
            if (eval(['isfield(s,''Investigator_' num2str(i) ''');']))
                %pf=regexp(['s.Investigator_'],'(.*?)[.](.*)$','tokens');
                eval(['s=rmfield(s, ''Investigator_' num2str(i) ''');']);
            end
        end
    end
s=order_struct(s);




% --- Executes on button press in btnPop.
function btnPop_Callback(hObject, eventdata, handles)
global lvl xml xml_list reds;

xml=pop_xml_from_file(xml);
lvl=1;list_struct(handles,xml);
xml_list{numel(xml_list)+1,1}='';
set(handles.xml_list,'String',xml_list);

% reds=~cellfun(@isempty,regexp(xml_list,'color="red"')); reds=find(reds);reds=flipdim(reds,1);
% set(handles.redslider,'Min',0);set(handles.redslider,'Max',numel(reds));
% set(handles.redslider,'SliderStep',[1/numel(reds) 1]);%1st number is step when hitting arrow, second when hitting bar
% %Slider is inverted: max value = line 1 and 1 = line max
% set(handles.redslider,'Value',numel(reds));redslider_Callback(handles.redslider, eventdata, handles);



% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
global xml_txt pretxt spch2 hdr;

fid=fopen([get(handles.xml_path,'String') filesep get(handles.xml_fn,'String') '_Readme.xml'],'w');
if fid<0,msgbox(sprintf('Problem Opening the xml File...\nMake Sure Path is Correct.'));return;end

%reprint template header
for i=1:size(hdr,2),fprintf(fid,'%s\n',char(hdr{i}));end

%special characters
%put back the special characters that were replaced in 'extract_lvlTag_from_text'
% spch2={'\[' '*b1*';'\]' '*b2*';'\{' '*c1*';'\}' '*c2*';'\(' '*p1*';'\)' '*p2*';'\=' '*eq*';'''' '*ap*';...
%         '\,' '*cm*';'\;' '*sc*';'\:' '*cl*';'\%' '*pc*';'\!' '*ep*';'\@' '*at*';'&' '*and*';};
for i=1:size(spch2,1),xml_txt=regexprep(xml_txt,spch2{i,2},spch2{i,1});end

fprintf(fid,'%s\n',char(pretxt{1}));
nbytes=fprintf(fid,'%s',xml_txt);
if nbytes<size(xml_txt,2), msgbox(sprintf('\n%s\n\n%s\n','Writing to File INCOMPLETE !','Contact Denis Pierrot (denis.pierrot@noaa.gov)'));end
fclose(fid);



% --- Executes on button press in dir_change.
function dir_change_Callback(hObject, eventdata, handles)

newdir=uigetdir(get(handles.xml_path,'String'),'Choose Directory');
if isdir(newdir), set(handles.xml_path,'String',newdir,'ToolTipString',newdir);end


% --- Executes on button press in import_xml.
function import_xml_Callback(hObject, eventdata, handles)
global xml xml_list line lvl xml_txt special hdr;
xml=struct;line=-1;xml_txt='';xml_list={};special=0;

[fullxml fullxml_p]=uigetfile('*.xml','Choose xml file');
if fullxml==0, return;end
if ~exist([fullxml_p fullxml],'file')==2, return;end
fid=fopen([fullxml_p fullxml],'r');

% searches for first tag
tag1=0;hdr='';
while(tag1==0 & ~feof(fid))
    hdr_in=textscan(fid,'%s',1,'delimiter','\n');
    hdr_in1 = strfind(hdr_in{1}, '<x_tags');
    tag1=~cellfun('isempty', hdr_in1);
    if ~tag1,hdr=[hdr hdr_in];end
end
hdr_in=regexprep(hdr_in{1},'<x_tags xml:space="preserve">','<x_tags>');

while ~feof(fid),fullxml_t=textscan(fid,'%s',10000,'delimiter','\n');end
fullxml_t=cat(1,hdr_in,fullxml_t{1});%adds first tag to rest of xml
fullxml_t=char(fullxml_t);
fclose(fid);

fullxml_t2='';for i=1:size(fullxml_t,1), fullxml_t2=strcat(fullxml_t2, fullxml_t(i,:));end

extract_lvlTag_from_text('xml',fullxml_t2);

if special, set([handles.specialc;handles.specialc_txt],'visible','on');end
lvl=1;list_struct(handles,xml);
xml_list{numel(xml_list)+1,1}='';
set(handles.xml_list,'String',xml_list);set(handles.btnPop,'Visible','on');



 
 
function extract_lvlTag_from_text(f,taxt)
global xml special pretxt spch spch2 pocs;
ff={};tags_1={'Investigator';'Ports_of_Call';'Variable'};%tags to add '_1' to, no matter what

if strcmp(f,'xml')>0, j=1;end

% if ~isempty(strfind(f,'Manufacturer_of'))
%     taxt
% end

%next line will find most outer tags (will skip one tag if same one is
%encountered after...e.g. Cruise_Info

[matc, tags]=regexp(taxt,'<([^ \f\n\r\t\v\<\>]+)?>(.*)</\1>','match','tokens'); %'?>' matches first '>'
%[matc, tags]=regexp(taxt,'<(\w+).*?>(.*)</\1>','match','tokens'); %'?>' matches first '>'
%*******************************************************************************
%*  NOTE: no ? inbetween tags because CDIAC xml has <Cruise_Info> tag inside a *
%*  <Cruise_Info> tag. Need to take outer one. This create a problem with      *
%*  repeating tags (like the ones in tags_1 variable). Taking the outer one    *
%*  would not separate them hence the expression in regexp is different for    *
%*  those (see below).                                                         *
%*****************************************************************************
  
%tags{1,1}{1},tags{1,1}{2}
if exist('j','var'),if (j==1), pretxt=regexp(taxt,['(.*)<' tags{1}{1} '>.*'],'tokens'); end,end

for i=1:size(tags,2)
    %Check for dots in tag name
    if (regexp(char(tags{i}{1}),'\.')<0),msgbox(sprintf('Dots are not allowed in Tags (%s).\nCorrect and Restart.',char(tags{i}{1})),'modal');uiwait;return;end 
    %replace special characters. Will be put back when saving.
%     spch2={'\[' '*b1*';'\]' '*b2*';'\{' '*c1*';'\}' '*c2*';'\(' '*p1*';'\)' '*p2*';'\=' '*eq*';'''' '*ap*';...
%         '\,' '*cm*';'\;' '*sc*';'\:' '*cl*';'\%' '*pc*';'\!' '*ep*';'\@' '*at*';'&' '*and*';};
    fld0=char(tags{i}{1});fld=fld0;    %  %replaces special characters (except dot) to build structure
    for ii=1:size(spch2,1),fld=regexprep(fld,spch2{ii,1},spch2{ii,2});end
    if ~strcmp(fld0,fld), special=1;end
    if ismember(fld,tags_1)  %'Investigator';'Ports_of_Call';'Variable'  tags that can have multiple instances
        %[matc2, tags2]=regexp(matc{i},'<(\w+).*?>(.*?)</\1>','match','tokens');
        [matc2, tags2]=regexp(matc{i},'<([^ \f\n\r\t\v\<\>]+)?>(.*?)</\1>','match','tokens');
        for ii=1:size(tags2,2)
            %rename field_ii
            if ii==1,fld0=fld;end
            fld=[fld0 '_' num2str(ii)];
            if strcmp(fld0,tags_1{2}),pocs{ii}=tags2{ii}{2};end 
            %put text in tag
            if ii==1, tags{i}{1}=fld;tags{i}{2}=tags2{ii}{2};
            else tags{end+1}{1}=fld;tags{end+1}{2}=tags2{ii}{2}; end
            eval([f '.' fld '='''';']);
            extract_lvlTag_from_text([f '.' fld],char(tags2{ii}{2}));
        end
        continue;
    else
        if eval(['isstruct(' f ')']), eval(['ff=fieldnames(' f ');']);end
        if ~isempty(ff)
            %wherei=~cellfun(@isempty,strcmp(ff,fld));
            wherei=strcmp(ff,fld);
            if sum(wherei)%if > 0, field fld already exists
                if j==1 %add suffix '_1' and remove old
                    eval(['[' f '.' fld '_1]=' f '.' fld ';']);
                    eval([f ' = rmfield(' f ',''' fld ''');']);
                end
                j=j+1;fld=[fld '_' num2str(j)];%add suffix '_i'
            end
        else
            j=1;
        end
    end
    eval([f '.' fld '='''';']);
    extract_lvlTag_from_text([f '.' fld],char(tags{i}{2}));
end
if size(tags,2)==0  %value of field assigned
    taxt=regexprep(taxt,'''','''''');%replaces any ' with '''
    %spch={'&' '&amp;';'<' '&lt;';'>' '&gt;';'oC' '&#176;C';'+-' '&#177;';'µ' '&#956;'};
    for i=1:size(spch,1),taxt=regexprep(taxt,spch{i,2},spch{i,1});end
    eval([f '=''' taxt ''';']);
end


function s=order_struct(s)
global maxtags;

f=fieldnames(s);
inv=[];

%find where investigator_i fields are...positions are put in 'inv' array ex: inv=[2,7,8]
for i=1:maxtags
    wherei=~cellfun(@isempty,strfind(f,['Investigator_' num2str(i)]));
    if ~isempty(find(wherei)),inv(i)=find(wherei);end
end

perm=1:numel(f);
%perm is array giving new order of previous elements(pe)  ex: perm=[pe1,%pe2, pe5, pe3, pe4] relocates element 5 in 3rd position.
for i=2:size(inv,2),  perm=cat(2,perm(1:inv(1)+i-2),perm(inv(i)), perm(inv(1)+i-1:inv(i)-1), perm(inv(i)+1:end) );end

s=orderfields(s,perm);



% --- Executes on slider movement.
function xml_list_font_change_Callback(hObject, eventdata, handles)
% hObject    handle to xml_list_font_change (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
lf=get(hObject,'value');%is either -1 or +1
if lf<0, lf=-1; else lf=1;end
set(hObject,'value',0);
pc=get(handles.xml_list_font,'string');pc=pc(1:end-2);
plf=str2num(pc);ofs=get(handles.xml_list,'FontSize')*100/plf;
%lfs=sprintf('%.0f %%',str2num(get(handles.xml_list_font,'string')));
lf=(plf+10*lf); lf=max([10,lf]);lf=min([200,lf]);lfs=sprintf('%.0f %%',lf);%keep in bounds
set(handles.xml_list_font,'String',lfs,'ToolTipString',lfs);
set(hObject,'Value',0);%resets so that arrows work again

set(handles.xml_list,'FontSize',lf/100*ofs);drawnow;


% --- Executes during object creation, after setting all properties.
function xml_list_font_change_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xml_list_font_change (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function mems = Explode_Struct(X,sname,lvl,fname,mems)
% function fns = expstruc(X,sname,lvl,fname,mems)
%
% Recursively explode a structure to find all members.
%
% INPUTS
% ------
% X;        structure to be exploded
% sname;    structure name (generated during recursive call, not needed to start)
% lvl;       structure level (1= top)
% fname;    name of the field (string, completely expanded)
% mems;     structure containing field names and data types
%
% written by Andy Clifton, October 2010.

if nargin == 1
    % generate us an empty name vector
    sname = inputname(1);
    fname = [];
    lvl = 1;
    mems = struct('name',{''},'text',{''});
end

switch class(X)
    case {'struct'} % still something to expand
        fns = fieldnames(X);
        for f = 1:numel(fns)
            % keep expanding the structure
            switch lvl
                case 1
                    mems = Explode_Struct(X.(fns{f}),sname,lvl+1, [sname '.' char(fns{f})],mems);
                otherwise
                    mems = Explode_Struct(X.(fns{f}),sname,lvl+1, [fname '.' char(fns{f})],mems);
            end
        end
    otherwise % then no more fields to expand
        %fprintf('...Level (%i): %s\n', lvl, fname)
        mems.name{end+1} = fname(length(sname)+2:end);%gets rid of original structure name
        mems.text{end+1} = X;
        return
end    


% --- Executes on slider movement.
function redslider_Callback(hObject, eventdata, handles)
% hObject    handle to redslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global xml_list reds;

if numel(xml_list)==0, return;end

%linen=get(hObject,'Value'); set(handles.linenum,'String',sprintf('Line %d',(get(hObject,'Max')+get(hObject,'Min'))-round(linen)));
linen=get(hObject,'Value');
if linen<1, linen=1;set(hObject,'Value',linen);end
if ~isempty(reds)
    set(handles.linenum,'String',sprintf('Line %d',reds(round(linen))));
    set(handles.xml_list,'Value',reds(round(linen)));
end
