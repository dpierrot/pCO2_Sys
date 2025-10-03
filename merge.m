function varargout = merge(varargin)

% MERGE M-file for merge.fig
%      MERGE, by itself, creates a new MERGE or raises the existing
%      singleton*.
%
%      H = MERGE returns the handle to a new MERGE or the handle to
%      the existing singleton*.
%
%      MERGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MERGE.M with the given input arguments.
%
%      MERGE('Property','Value',...) creates a new MERGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before merge_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to merge_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help merge

% Last Modified by GUIDE v2.5 22-Jun-2009 14:32:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @merge_OpeningFcn, ...
                   'gui_OutputFcn',  @merge_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
               
try
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end
    
    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
catch
    disp(lasterr);
    if strcmp(lasterr,'Attempt to reference field of non-structure array.')~=1
        msgbox(sprintf('Error has occurred\n(%s)',lasterr),'Error','error','modal');
    uiwait;
    end
end

% End initialization code - DO NOT EDIT

%    try
%         [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
%     catch
%         disp(lasterr);
%         if strcmp(lasterr,'Attempt to reference field of non-structure array.')~=1
%             msgbox(sprintf('Error has occurred\n(%s)',lasterr),'Error','error','modal');
%             uiwait;
%         end
%     end
% --- Executes just before merge is made visible.
function merge_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to merge (see VARARGIN)

% Choose default command line output for merge
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes merge wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%Denis Code starts here..................
global IMM2 sData listf vyday pimport Nonecol;%listf_file listf_calc cmh cmhf cmhc cmhn Nonecol;
global monitors Dmonitor;

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

listf=strtrim(listf);
set(handles.merged, 'String',listf);if (vyday~= Nonecol), set(handles.merged,'value',vyday);end
set(handles.importd,'visible','off');
set(handles.merge_select,'visible','off');
set(handles.savefile,'visible','off');
set(handles.merge_limit_t,'visible','off');
set(handles.merge_limit,'visible','off');
set(handles.rdoref,'visible','off');
set(handles.rdomin,'visible','off');
set(handles.uimethod,'visible','off');
set(handles.rdolinear,'visible','off');
set(handles.rdonearest,'visible','off');
pimport=[];
if exist('sData','var') & isstruct(sData)
    if ~isempty(sData.data)
        mmm=questdlg('Keep Previously Imported Data?','Merge','Yes','No','No');
        
        if strcmp(mmm,'Yes')
            pimport=1;  IMM2=[];
            set(handles.importd,'String',sData.hdr(:,:));
            set(handles.importd,'visible','on');
            set(handles.importdtxt,'visible','on');
            set(handles.merge_select,'visible','on');
            set(handles.savefile,'visible','on');
            set(handles.merge_limit_t,'visible','on');
            set(handles.merge_limit,'visible','on');
            set(handles.rdoref,'visible','on');
            set(handles.rdomin,'visible','on');
            set(handles.uimethod,'visible','on');
            set(handles.rdolinear,'visible','on');
            set(handles.rdonearest,'visible','on');
        else
            sData.hdr={};sData.data=[];pimport=[]; IMM2=[];
        end
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = merge_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




function select_merge_file_Callback(hObject, eventdata, handles)
global sData ifpath datap  pimport var0 var1;

if ~isempty(pimport)% some data has already been imported and will be reused
    mmm=questdlg(sprintf('This will replace the data already in memory.\n\nContinue anyway?'),'Merge File Selection','Yes','No','Yes');
    if strcmp(mmm,'Yes'),pimport=[];sData.hdr={};sData.data=[];end
end
if isempty(pimport)% some data has already been imported and will be reused
    [ifname,ifpath] = uigetfile('*.txt; *.csv; *.tsv; *.mat','Select Data File(s)','MultiSelect','off',datap);
    if (ifpath==0),return;end
    [~,fn,fe]=fileparts(ifname);
    fid=fopen([ifpath,ifname],'r');
    line1=fgetl(fid);line2=strtrim(line1);mmm='';
    if length(line1)~=length(line2), mmm=questdlg(sprintf('File seems to have extra spaces at the end of each line.\n\nDo you want the program to correct it?'),'Import Data','Yes','No','Yes');end
    if strcmp(mmm,'Yes')
        fidout=fopen([ifpath,fn,'out',fe],'w');
        fprintf(fidout,'%s\n',line2);line1=line2;
        while ~feof(fid)
            linein=fgetl(fid); linein=strtrim(linein);fprintf(fidout,'%s\n',linein);
        end
        fclose(fidout);    fclose(fid);
        movefile([ifpath,fn,'out',fe],[ifpath,ifname]);
    else
        fclose(fid);
    end
   % line1 = char(textread([ifpath,ifname], '%[^\n]', 1));
    queststr=sprintf(['Does the file have a HEADER LINE?\n\n'...
                      'First line is:\n%s'],line1);
    hdrline= strcmp('Yes',questdlg(queststr,'INPUT','Yes','No','Yes'));
    if importfilenew([ifpath  ifname],hdrline),return;end

end


set(handles.importd,'String',sData.hdr(:,:));

set(handles.importd,'visible','on');
set(handles.importdtxt,'visible','on');
set(handles.merge_select,'visible','on');
set(handles.savefile,'visible','on');
set(handles.merge_limit_t,'visible','on');
set(handles.merge_limit,'visible','on');
set(handles.rdoref,'visible','on');
set(handles.rdomin,'visible','on');
set(handles.uimethod,'visible','on');
set(handles.rdolinear,'visible','on');
set(handles.rdonearest,'visible','on');


function merge_select_Callback(hObject, eventdata, handles)
global sData DataA M IMM2 s oxi ixi del ixiv oxiv;

imethod=strvcat('linear','nearest');

[s,v] = listdlg('PromptString','Select variables to merge:',...
                'SelectionMode','multiple',...
                'ListString',sData.hdr);
if ~v, return; end %exits if Cancel
set(handles.plotl,'value',1);
set(handles.plotl,'String',sData.hdr(s(:)));
drawnow;

mlimit=get(handles.merge_limit,'String');
mlimit=str2num(mlimit);
if get(handles.rdomin,'value')
    mlimit=mlimit/60./24.;
end

M=[];
%oxi = out x index    ixi= in x index
oxi=str2num(char(DataA(2:end,get(handles.merged,'value'))));ixi=str2num(char(sData.data(:,get(handles.importd,'value'))));
oxiv = isfinite(oxi);oxiv=find(oxiv);oxi=oxi(oxiv);%gets rid of NaN and Inf
ixiv = isfinite(ixi);ixiv=find(ixiv);ixi=ixi(ixiv);

if numel(oxi)~= numel(unique(oxi))
    mmm=msgbox(sprintf('%s: The Reference used to merge the data has non-unique values. \nCheck the reference data in the %s.','pCO2 files'),'Import Data','modal');
    CenterWindow(handles,mmm,handles.figure1) ;uiwait(mmm);
    erreur=1;
    return ;
end
if numel(ixi)~= numel(unique(ixi))
    msgbox(sprintf('%s: The Reference used to merge the data has non-unique values. \nCheck the reference data in the %s.','merging file'),'Import Data','modal');
    CenterWindow(handles,mmm,handles.figure1) ;uiwait;
    erreur=1;
    return ;
end
M(2:length(oxi)+1,1)=oxi;
% M(2:length(oxi)+1,length(s)+2)=1;

M(2:length(oxi)+1,length(s)+2)=mean(str2num(char(sData.data(ixiv,s(1)))));
if isempty(ixi) | isempty(oxi)
    msgbox('Select Reference Variables to do the merging','Merge Error');
    CenterWindow(handles,mmm,handles.figure1) ;uiwait;
    return;
end
[k,d] = dsearchn(ixi,oxi);  %returns the indices k of the closest points in ixi for each point in oxi. d = the distances d to the closest points.
im=get(handles.rdonearest,'Value')+1;
del=[];
for i=2:length(s)+1
    y=interp1(ixi,str2num(char(sData.data(ixiv,s(i-1)))),oxi,imethod(im,:));
    if sum(isnan(y))==length(y) % no common data to import
       del=[del i]; %will delete column to import
       M(2:end,i)=NaN;
    else    
        M(2:length(y)+1,i)=y;
        ol=find((d>mlimit)& ~isnan(M(2:end,i)));%out of limit
        M(ol+1,i)=NaN;
        oobi=find((d<=mlimit) & (oxi<min(ixi) | oxi>max(ixi)));% out of bounds indices where interpolation doesn't occur and still within limits
        M(oobi+1,i)=str2num(char(sData.data(k(oobi),s(i-1)))); % k(oobi) are the indices of NEAREST points
    end
end
hwb=waitbar(1,'Plotting Data...Please wait!');
plot(handles.maxes1,ixi,str2num(char(sData.data(ixiv,s(1)))),'ok',oxi,M(2:end,2),'.r',oxi,M(2:end,length(s)+2),'.b');
legend(handles.maxes1,'All imported Data','Merged Data','pCO2 Data Range','Location','Best');
close(hwb);
    
% M(:,del)=[];%deletes columns with no data
% s(del-1)=[];
if length(s)~=length(del)
    % creates cell array IMM2 containing merged data, including headers, and sets it to varargout
    IMM2=0; IMM2=cell(size(M,1),length(s));
    for i= 1:length(s),    IMM2(1,i)=strtrim(sData.hdr(s(i)));end
    for j= 2:size(M,1), for i= 1:length(s),    IMM2(j,i)={strtrim(num2str(M(j,i+1),'%-1.3f'))};end, end
    %empty values and '-9' will be replaced by '-999'
    IMM2(find(strcmp(IMM2,'-9.000')==1))={'-999'};
    IMM2=strrep(IMM2,'NaN','-999');
end


function plotl_Callback(hObject, eventdata, handles)
global M s sData oxi ixi ixiv ;

pv=get(handles.plotl,'value');
M(2:length(oxi)+1,length(s)+2)=mean(str2num(char(sData.data(ixiv,s(pv)))));

plot(handles.maxes1,ixi,str2num(char(sData.data(ixiv,s(pv)))),'ok',oxi,M(2:end,1+pv),'.r',oxi,M(2:end,length(s)+2),'.b');
legend(handles.maxes1,'All imported Data','Merged Data','pCO2 Data Range','Location','Best');


function savefile_Callback(hObject, eventdata, handles)
global   sData DataA M s; %#ok<NUSED>
global ifpath

if (ifpath==0), ifpath='c:\'; end
[sfname,sfpath] = uiputfile('*.csv','Save Merged Data In',ifpath);
% [sfname,sfpath] = uigetfile('*.txt; *.mat','Select Data File(s)','MultiSelect','off',ifpath);
if (sfpath==0),return;end
% if exist([sfpath,sfname],'file')
%     mmm=questdlg('File already exists...Overwrite?','Saving','Yes','Cancel','Cancel');
%     if strcmp(mmm,'Cancel'), return; end
% end

spr=', ';
% hdr(1,:)=char(DataA(1,get(handles.merged,'value')));
% for i=2:length(s)+1
%     hdr(i)=sData.hdr(s(i-1),:);
% end
hdr=strvcat(char(DataA(1,get(handles.merged,'value'))),strtrim(char(sData.hdr(s))));
% hdr(2,:)='hder2';hdr(3,:)='hder3';hdr(4,:)='hder4';hdr(5,:)='hder5';
catstr0='MM=cat(2,M(2:end,1)';
fstr0='%-1.5f';
cathdr0='HH=[hdr(1,:)';
catstr1=');';
cathdr1='];';
for i=2:length(s)+1
    catstr0=[catstr0 ',M(2:end,'  num2str(i)  ')'];
    fstr0=[fstr0 ', %-1.3f'];
    cathdr0=[cathdr0 ' spr hdr('  num2str(i)  ',:)'];
end
eval([catstr0 catstr1]);
eval([cathdr0 cathdr1]);
fstr0=[fstr0 '\n'];

fid=fopen ([sfpath,sfname],'w');

% dlmwrite([sfpath,sfname],HH,'');
% dlmwrite([sfpath,sfname],MM,'-append','delimiter',',','precision',2);
fprintf(fid, '%s\n',HH);
fprintf(fid,fstr0,MM');
fclose(fid);

function erreur=importfile(fileToRead,HEADERLINES)
global sData;
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read
%  Auto-generated by MATLAB on 21-Aug-2008 17:55:01
erreur=0;
fid=fopen(fileToRead,'r');finfo=dir(fileToRead);
if fid>0,sData=struct('hdr','','data',[]);end%resets sData
lines=fgetl(fid); %reads header of 1st file
lines2=fgetl(fid); %reads 1st line of data
fclose(fid);
ncol0=length(regexp(lines,','));  ncol1=length(regexp(lines,'\t'));
if ncol0>ncol1, DELIMITER=sprintf(',');delimtxt='Comma'; else DELIMITER=sprintf('\\t');delimtxt='Tab';end
ncol0=length(regexp(lines,DELIMITER)); ncol2=length(regexp(lines2,DELIMITER)); 
if ncol0~=ncol2
    if HEADERLINES>0, junk='Header';junk2='1st'; else junk='1st';junk2='2nd';end
    msgbox(sprintf('%s line has different number of columns than %s line...\nCheck number of delimiters (%s)',junk,junk2,delimtxt),'Import Data','modal');
    CenterWindow(handles,mmm,handles.figure1) ;uiwait;
    erreur=1;
    return ;
end
%HEADERLINES = 0;

% Import the file
    fid=fopen(fileToRead,'r');

    if HEADERLINES
        line1=fgetl(fid);line2=fgetl(fid);
        line2=regexprep(line2,[DELIMITER DELIMITER],[DELIMITER '-999' DELIMITER]);
        line2=regexprep(line2,[DELIMITER DELIMITER],[DELIMITER '-999' DELIMITER]);  %needed when more than 2 consecutive delim      
        if line2(end)==','   line2=[line2 '-999'];  end
        ins=strsplit(line2,DELIMITER);  nums=find(cellfun(@isempty,regexp(strtrim(ins),'[^0-9\.-+]'))); %indices of numeric columns
        ins=strsplit(line1,DELIMITER);  sData.hdr=cellstr(strtrim(char(ins(nums))));
        ins=strsplit(line2,DELIMITER);  sData.data(1,:)=str2num(char(ins(nums)))';
    else
        line1=fgetl(fid);
        line1=regexprep(line1,[DELIMITER DELIMITER],[DELIMITER '-999' DELIMITER]);
        line1=regexprep(line1,[DELIMITER DELIMITER],[DELIMITER '-999' DELIMITER]);  %needed when more than 2 consecutive delim        
        if line1(end)==','   line1=[line1 '-999'];  end
        ins=strsplit(line1,DELIMITER);  nums=find(cellfun(@isempty,regexp(strtrim(ins),'[^0-9\.-+]'))); %indices of numeric columns
        sData.data(1,:)=str2num(char(ins(nums)))';
        sData.hdr={};
        for i=1:size(nums,2),sData.hdr(i,:)={['Column ' num2str(nums(i))]};end        
    end
    
    ln=2;hwb=waitbar(0,'Loading Data...');
    while ~feof(fid)
        if mod(ln,100)==0
            waitbar(ftell(fid)/finfo.bytes,hwb,sprintf('Loading Data...%d lines.',ln));drawnow;
        end
        linein=fgetl(fid);
        linein=regexprep(linein,[DELIMITER DELIMITER],[DELIMITER '-999' DELIMITER]);
        linein=regexprep(linein,[DELIMITER DELIMITER],[DELIMITER '-999' DELIMITER]);
        if linein(end)==','   linein=[linein '-999'];  end
        ins=strsplit(linein,DELIMITER);
        sData.data(ln,:)=str2num(char(ins(nums)));ln=ln+1;
    end
    close(hwb);
    fclose(fid);


% sData = importdata(fileToRead, DELIMITER, HEADERLINES);
 
% Create new variables in the base workspace from those fields.
% vars = fieldnames(sData);
% for i = 1:length(vars)
%     assignin('base', vars{i}, sData.(vars{i}));
% end

function erreur=importfilenew(fileToRead,HEADERLINES)
global sData;
%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/denis.pierrot/Documents/SOOP CO2/Daryin cruises to check/HB2403Leg3.4txt.txt
%
% Auto-generated by MATLAB on 07-Dec-2023 17:05:34
erreur=0;

fid=fopen(fileToRead,'r');finfo=dir(fileToRead);
if fid>0,sData=struct('hdr','','data',[]);end%resets sData
line1=fgetl(fid); %reads header of 1st file
line2=fgetl(fid); %reads 1st line of data

ncol0=length(regexp(line1,','));  ncol1=length(regexp(line1,'\t'));
if ncol0>ncol1, DELIMITER=sprintf(',');delimtxt='Comma'; else DELIMITER=sprintf('\\t');delimtxt='Tab';end
ncol0=length(regexp(line1,DELIMITER)); ncol2=length(regexp(line2,DELIMITER)); 
if ncol0~=ncol2
    if HEADERLINES>0, junk='Header';junk2='1st'; else junk='1st';junk2='2nd';end
    msgbox(sprintf('%s line has different number of columns than %s line...\nCheck number of delimiters (%s)',junk,junk2,delimtxt),'Import Data','modal');
    CenterWindow(handles,mmm,handles.figure1) ;uiwait;
    erreur=1;
    return ;
end

%create header

if HEADERLINES
    line2=regexprep(line2,[DELIMITER DELIMITER],[DELIMITER '-999' DELIMITER]);
    line2=regexprep(line2,[DELIMITER DELIMITER],[DELIMITER '-999' DELIMITER]);  %needed when more than 2 consecutive delim
    if line2(end)==','   line2=[line2 '-999'];  end
    ins=strsplit(line2,DELIMITER);  nums=find(cellfun(@isempty,regexp(strtrim(ins),'[^0-9\.-+]'))); %indices of numeric columns
    ins=strsplit(line1,DELIMITER);  sData.hdr=cellstr(strtrim(char(ins(nums))));
    % ins=strsplit(line2,DELIMITER);  sData.data(1,:)=str2num(char(ins(nums)))';
else
    line1=regexprep(line1,[DELIMITER DELIMITER],[DELIMITER '-999' DELIMITER]);
    line1=regexprep(line1,[DELIMITER DELIMITER],[DELIMITER '-999' DELIMITER]);  %needed when more than 2 consecutive delim
    if line1(end)==','   line1=[line1 '-999'];  end
    ins=strsplit(line1,DELIMITER);  nums=find(cellfun(@isempty,regexp(strtrim(ins),'[^0-9\.-+]'))); %indices of numeric columns
    % sData.data(1,:)=str2num(char(ins(nums)))';
    sData.hdr={};
    for i=1:size(nums,2),sData.hdr(i,:)={['Column ' num2str(nums(i))]};end
end

fclose(fid);

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", ncol0+1);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = DELIMITER;

% Specify column names and types
%opts.VariableNames = sData.hdr;


%opts.VariableTypes = ["double", "char", "char", "double", "double", "double", "double", "double", "double", "double", "double", "double", "char", "char", "double", "double", "char", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.VariableNamingRule="preserve";
% Specify variable properties
% opts.WhitespaceRule="preserve";
% opts.EmptyFieldRule="auto";

% opts = setvaropts(opts, sData.hdr, "EmptyFieldRule", "auto");
% opts = setvaropts(opts, sData.hdr, "WhitespaceRule", "preserve");

% opts = setvaropts(opts, ["Date", "Time", "TSGExternalTemp", "TSGInternalTemp", "TSGSalinity"], "WhitespaceRule", "preserve");
% opts = setvaropts(opts, ["Date", "Time", "TSGExternalTemp", "TSGInternalTemp", "TSGSalinity"], "EmptyFieldRule", "auto");

% Import the data
hwb=waitbar(0,'Loading Data...');
%tablea = readtable(fileToRead, opts);
tablea = readtable(fileToRead, opts);

sData.data=tablea{:,nums};
close(hwb);
sData.data(find(cellfun('isempty',sData.data)))= cellstr('NaN');

%% Clear temporary variables
clear opts line1 line2 nums ncol0 ncol1 ncol2 ins finfo fid delimtxt 



function close_cancel(hObject, eventdata, handles)
global IMM IMM2;
    IMM=[];
    if handles.donebut.Value,  IMM=IMM2;end
    close all;

