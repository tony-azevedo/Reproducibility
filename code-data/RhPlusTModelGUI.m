function varargout = RhPlusTModelGUI(varargin)
% RHPLUSTMODELGUI M-file for RhPlusTModelGUI.fig
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RhPlusTModelGUI

% Last Modified by GUIDE v2.5 06-Oct-2011 12:49:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RhPlusTModelGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RhPlusTModelGUI_OutputFcn, ...
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


% --- Executes just before RhPlusTModelGUI is made visible.
function RhPlusTModelGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RhPlusTModelGUI (see VARARGIN)

condition_selector_Callback(handles.condition_selector, eventdata, handles)

% set Default values
data = get(handles.params,'Data');
data{1,1} = 30;
data{2,1} = 100;
data{3,1} = 150;
data{3,2} = 350;
data{3,3} = 100;
data{4,1} = 1;
data{5,1} = 4;
data{5,2} = 5;
data{5,3} = 1;
data{6,1} = 0;
data{7,1} = 0;
data{8,1} = .1;
data{8,2} = .5;
data{8,3} = .1;
data{9,1} = 100;
data{9,2} = 300;
data{9,3} = 100;
data{10,1} = 0;
data{10,2} = .01;
data{10,3} = .005;
data{11,1} = 0;
data{12,1} = 0;
data{13,1} = 0;
set(handles.params,'Data',data);

handles = guidata(hObject);
% Choose default command line output for RhPlusTModelGUI
handles.output = hObject;

handles.user.scmp = java.util.HashMap;

% Update handles structure
guidata(hObject, handles);
params_CellEditCallback(hObject, eventdata, handles);

% --- Outputs from this function are returned to the command line.
function varargout = RhPlusTModelGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ------------------ Interactive functions ------------------
function condition_selector_Callback(hObject, eventdata, h) %#ok<*DEFNU>
load LinearityResps

conditions = get(hObject,'String');
condition = conditions{get(hObject,'Value')};
tokens = regexp(condition,'\s','split');
geno = tokens{1};
temp = tokens{end};
if length(tokens)>2
    geno = [geno,'Het'];
end

eval(sprintf('load %sSingles%s',geno,temp));
%eval(sprintf('h.user.target = %s%sSingles;',geno,temp))

target.MeanTPeakRatio = Target.MeanTPeakRatio;
target.SEMTPeakRatio = Target.SEMTPeakRatio;
target.MeanVarWidth = Target.MeanVarWidth;
target.SEMVarWidth = Target.SEMVarWidth;
target.MeanCVArea = Target.MeanCVArea;
target.SEMCVArea = Target.SEMCVArea;

h.user.target = target;
%eval(sprintf('h.user.munorm = %sAveNormSingle%sC;',geno,temp))
eval(sprintf('h.user.munorm = Target.AveSingle/max(Target.AveSingle);'))
eval(sprintf('h.user.varnorm = %sVarNormSingle%sC;',geno,temp))
h.user.time = 1:length(h.user.munorm);

rn = {'CV','TPeak','VarW','Like'};
data = {h.user.target.MeanCVArea;h.user.target.MeanTPeakRatio;h.user.target.MeanVarWidth} ;
set(h.outputdata,'rowname',rn,'data',data);

guidata(hObject, h);

function run_rh_Callback(run, evnt, h)
simfunc = @RhCascadeModelTA;
runSimulations(run,evnt,h,simfunc);

function run_t_Callback(run, evnt, h)
simfunc = @RhPlusTransducinCascadeModelTA;
runSimulations(run,evnt,h,simfunc);

function plotlist_Callback(hObject, eventdata, h)
h = guidata(hObject);
if ~isfield(h.user,'simulationfunction')
    return
end
contents = get(hObject,'String');
func = contents{get(hObject,'Value')};
if strcmp(func,'ShowExamples')
    if isfield(h,'togglepanel'),set(h.togglepanel,'visible','off');end
    if isfield(h,'togglepanel'),set(h.plotchecks,'visible','on');end
else
    if isfield(h,'togglepanel'),set(h.togglepanel,'visible','on');end
    if isfield(h,'togglepanel'),set(h.plotchecks,'visible','off');end
end
eval(sprintf('%s(hObject,[],h);',func));
h = guidata(gcbf);
populateOutput(h.outputdata,eventdata,h);

function generic_slider_callback(hObject,eventdata)
h = guidata(gcbf);
sldrs = findobj('parent',h.sliderpanel,'style','slider');
cnt = zeros(size(sldrs));
simloopsnames = h.user.simloopsnames;

for c = 1:length(sldrs)
    vec = get(sldrs(c),'user');
    val = get(sldrs(c),'value');
    [~, ind] = min(abs(val-vec));
    cnt(strcmp(simloopsnames,get(sldrs(c),'tag'))) = ind;
end
h = makeKeysFromSliders(h,cnt);
guidata(hObject,h);
plotlist_Callback(h.plotlist, [], h);

function toggles_callback(hObject,eventdata) %#ok<*INUSD>
h = guidata(hObject);
plotlist_Callback(h.plotlist,[],h);

function generic_check_Callback(hObject, eventdata, h) %#ok<*INUSL>
plotlist_Callback(h.plotlist, [], h);

% ---- Params for the database are defined here -----
function params_CellEditCallback(hObject, eventdata, h)
% find the Inputs to the simulation

data = get(h.params,'data');
rn = get(h.params,'rowname');
cn = get(h.params,'columnname');

% find those that have three entries
params = nan(size(data));
loops = false(size(rn));
for r = 1:length(rn)
    for c = 1:length(cn)
        if ~isempty(data{r,c})
            if c==3 && data{r,c}==0
                error('Third Column entries must be NaN')
            end
            params(r,c) = data{r,c};
        end
    end
    if ~sum(isnan(params(r,:)))
        loops(r) = true;
    end
end

if isfield(h.user,'simulationfunction')
    simfunc = h.user.simulationfunction;
else
    simfunc = @RhPlusTransducinCascadeModelTA;
end
switch func2str(simfunc)
    case 'RhCascadeModelTA'
        rownames = {'Freq';'# Resp';'Rh Tau';'Rh Decay';'Rh Steps';'Rh Cmpr';'I Cmpr'};
        paramnames = {'FreqCutoff';    'NumResponses';   'RhShutoffFact';    'RhDecayFact';    'RhSteps';    'RhodopsinCompression';    'ResponseCompression'};
    case 'RhPlusTransducinCascadeModelTA'
        rownames = {'Freq';'# Resp';'Rh Tau';'Rh Decay';'Rh Steps';'Rh Cmpr';'Rh Dtrm';'T Rate';'T Decay';'T Cmpr';'T* Cmpr';'T Dtrm';'I Cmpr'};
        paramnames = {'FreqCutoff';'NumResponses';'RhShutoffFact';'RhDecayFact';'RhSteps';'RhodopsinCompression';'DeterministicRhModel';'TransRate';'TransTConst';'TransCompression';'TransRateCompression';'DeterministicTrModel';'ResponseCompression'};
end

loopinds = loops*1;
simloopsnames = {};
for r = 1:length(rn)
    n = find(strcmp(rownames,rn(r)));
    if ~isempty(n)
        sims.(paramnames{n}) = params(r,1);
        if loops(r)
            sims.(paramnames{n}) = params(r,1):params(r,3):params(r,2);
            loopinds(r) = length(sims.(paramnames{n}));
            simloopsnames{end+1} = paramnames{n}; %#ok<AGROW>
        end
    end
end
loopinds(loopinds<=1) = 0;
h.user.loopinds = loopinds;
h.user.sims = sims;
h.user.simloopsnames = simloopsnames;
maxcnt = loopinds(loopinds>0);

set(h.repetitions,'string',sprintf('%d',prod(maxcnt,1)));

guidata(hObject,h);

function verbose_DB_Callback(obj, eventdata, h)
switch get(h.verbose_DB,'Value')
     case 1, set(h.verbose_DB,'String','DB print off?')
     case 0, set(h.verbose_DB,'String','Database print?')
end
guidata(obj,h);

function most_likely_Callback(hObject, eventdata, handles)
h = guidata(gcbf);

key = h.user.simulations.mostLikely;
pnames = fieldnames(h.user.simulations.parameters);

sldrs = findobj('parent',h.sliderpanel,'style','slider');
for c = 1:length(sldrs)
    p = strcmp(pnames,get(sldrs(c),'tag'));
    ind = key(p);
    set(sldrs(c),'value',h.user.simulations.parameters.(pnames{p})(ind));
end
guidata(hObject,h);
generic_slider_callback(hObject,h);


% ------------------ Interactive functions ------------------


% -------------------- Simulate across all parameters --------------------

function runSimulations(run,evnt,h,simfunc)
h = rmfield(h,'user');
condition_selector_Callback(h.condition_selector, evnt, h);
h = guidata(run);

h.user.simulationfunction = simfunc;

h.user.scmp = java.util.HashMap;

field = {'sliderpanel','togglepanel'};
for name = 1:length(field)
    if isfield(h,field{name})
        if ishandle(h.(field{name}))
            delete(h.(field{name}));
        end
        h = rmfield(h,field{name});
    end
end

guidata(run,h);
h = guidata(run);

params_CellEditCallback(h.params,[],h);
h = guidata(run);

simloopsnames = h.user.simloopsnames;
sims = h.user.sims;
loopinds = h.user.loopinds;

% ---- Load or create a simulations database
simulations = rhTrSimulations(sims,h.user.simulationfunction);
switch get(h.verbose_DB,'value')
    case 0, simulations = simulations.verboseOff;
    case 1, simulations = simulations.verboseOn;
end
simulations.isVerbose;
%

clear SimulateCondition;
% create simulation skeleton, first parameters that change
field = fieldnames(sims);
for name = 1:length(field)
    SimulateCondition.(field{name}) = sims.(field{name})(1);
end
% add to simulation skeleton parameters that don't change
SimulateCondition.SamplingInterval = 1;         % ms
SimulateCondition.EpochPts = length(h.user.munorm);
SimulateCondition.NumPtsToPeak = 100;
SimulateCondition.Resp = h.user.munorm;
SimulateCondition.target = h.user.target;
SimulateCondition.SimFunc = simfunc;         % run rh or rhplustr function

% loop over possible loops
maxcnt = loopinds(loopinds>0)';
if isempty(maxcnt), maxcnt = 1; end

total = prod(maxcnt);
waitstr = sprintf('Running %d Simulations',total);
distribute = get(h.distribute,'value');
if distribute
    tic, fprintf('Creating Job: '); 
    jm = findResource('scheduler', 'configuration', 'rieke-server');
    job = createJob(jm, 'FileDependencies',...
        {'DoOneCascadeModelTA.m','RhCascadeModelTA.m','RhPlusTransducinCascadeModelTA.m'});
    waitstr = sprintf('Creating %d Tasks',total);
    toc, fprintf('Creating %d Tasks: ',total);
end

cnt = ones(size(maxcnt));
x = 0; wb = waitbar(x/total,waitstr);
tic
while ~isempty(cnt)
    % specify new parameter values
    for sl = 1:length(simloopsnames)
        SimulateCondition.(simloopsnames{sl}) = sims.(simloopsnames{sl})(cnt(sl));
    end
    
    if ~simulations.simulationExists(SimulateCondition)
        
        % Run Simulations
        if distribute
            createTask(job, ...                     % the job to add the task to
                @DoOneCascadeModelTA, ...           % function to be run
                1, ...                              % number of output arguments returned
                {SimulateCondition});             % input arguments

        else 
            features = DoOneCascadeModelTA(SimulateCondition);
            simulations = simulations.putSimulation(features.SimulateCondition,features);
        end
        
    else
        % return the features of a particular simulation
        simulations.getFeatures(SimulateCondition)
    end
    
    guidata(run,h);

    cnt = increasecnt(cnt,maxcnt);
    x = x+1;
    waitbar(x/total,wb);
end
toc
delete(wb);
% submit the job to the cluster
if distribute
    tic, fprintf('Job Submitted: ');
    submit(job);    
    waitForState(job);    
    SimReturns = getAllOutputArguments(job);
    destroy(job);
    toc
    waitstr = sprintf('Putting %d Simulations in DB',total);
    wb = waitbar(0/total,waitstr);
    tic, fprintf('Loading database: ');
    for i = 1:length(SimReturns)
        features = SimReturns{i};
        simulations = simulations.putSimulation(features.SimulateCondition,features);
        waitbar(i/total,wb);
    end
    toc
    delete(wb);
end

simulations = simulations.writeSimulations;
h.user.simulations = simulations;

% Setup meta plot lists
plist = {'TPeakRatioVs';'VarWidthVs';'CVVs';'LikelihoodVs';'LikelihoodSurfs';'ShowExamples'};
set(h.plotlist,'string',plist);
set(h.plotlist,'value',1);

guidata(run,h);
h = guidata(run);
h = makeSliders(run,[],h);
h = makeVariableToggles(run,[],h);
guidata(run,h);
generic_slider_callback(h.plotlist);

% -------------------- Simulate across all parameters --------------------



% ------------------ Internal functions ------------------
function h = makeSliders(hObject,eventdata,h)

simloopsnames = h.user.simloopsnames;

% Set initial conditions (best from the TPeakRatio, VarWidth, and CV)
% change this when you find out how to make the best condition;
likely = h.user.simulations.mostLikely;

% Setup sliders for each loop (max, means, step)
if ~isempty(simloopsnames);
    h.sliderpanel = uipanel('parent',gcbf,'units','normalized','position',[.63, .81, .3255, .18],'tag','Adjust Values');
    set(h.sliderpanel,'units','pixels')
    pos = get(h.sliderpanel,'position');
    set(h.sliderpanel,'units','normalized')
    for sn = 1:length(simloopsnames)
        vec = h.user.sims.(simloopsnames{sn});
        pn = strcmp(fieldnames(h.user.simulations.parameters),simloopsnames{sn});
        user = vec(2)-vec(1);
        dist = vec(end)-vec(1);
        step = user/dist-eps;
        uicontrol('parent',h.sliderpanel,...
            'style','slider',...
            'units','pixels',...
            'position',[70,pos(4)-10-(sn)*20, pos(3)-70-20,20 ],...
            'min',min(vec),...
            'max',max(vec),...
            'sliderstep',[step,step],...
            'user',vec,...
            'value',vec(likely(pn)),...
            'tag',simloopsnames{sn},...
            'Callback', {@generic_slider_callback});
        % write text
        uicontrol('parent',h.sliderpanel,...
            'style','text',...
            'units','pixels',...
            'position',[10,pos(4)-10-(sn)*20, 60,20 ],...
            'String',simloopsnames{sn},...
            'Fontsize',7);
    end
end
guidata(hObject,h);

function h = makeVariableToggles(hObject,eventdata,h)
simloopsnames = h.user.simloopsnames;

% Setup radio for each looped variable (to use as xaxis in some plots)
if ~isempty(simloopsnames);
    set(h.canvas,'units','pixels')
    pos = get(h.canvas,'position');
    set(h.sliderpanel,'units','normalized')
    h.togglepanel = uibuttongroup('parent',gcbf,...
        'units','pixels',...
        'position',[300, pos(2)+pos(4),pos(3)-87,27],...
        'tag','Variables',...
        'bordertype','none',...
        'SelectionChangeFcn',{@toggles_callback});
    for sn = 1:length(simloopsnames)
        uicontrol('parent',h.togglepanel,...
            'style','togglebutton',...
            'units','normalized',...
            'position',[(sn-1)*.1,0,.1,1 ],...
            'string',simloopsnames{sn},...
            'tag',simloopsnames{sn});
    end
end
guidata(hObject,h);

function TPeakRatioVs(hObject,evnt,h)
feature = 'TPeakRatio';
compare = h.user.target.MeanTPeakRatio;
plotMeasureVs(hObject,h,feature,compare);

function VarWidthVs(hObject,evnt,h)
feature = 'VarWidth';
compare = h.user.target.MeanVarWidth;
plotMeasureVs(hObject,h,feature,compare);

function CVVs(hObject,evnt,h)
feature = 'CV';
% Only if there is a CV 
if isfield(h.user.target,'MeanCVArea') 
    compare = h.user.target.MeanCVArea;
    plotMeasureVs(hObject,h,feature,compare);
else
    plotMeasureVs(hObject,h,feature);
end    

function LikelihoodVs(hObject,evnt,h)
feature = 'Likelihood';
plotMeasureVs(hObject,h,feature);

function LikelihoodSurfs(hObject,evnt,h)
% make grid of plots
delete(get(h.canvas,'children'));
ps = length(h.user.simloopsnames)-1;
hgrids = ceil(sqrt(ps));
if hgrids == 0
    return % for now
end
vgrids = ceil(ps/hgrids);
gridmax = [vgrids,hgrids];
gridpnt = [1,1];

hoff = .025; hbrd = .04;
voff = .03; vbrd = .05;
wdth = (1-hoff)/hgrids;
hght = (1-voff)/vgrids;

key = makeKey(h.user.currentkey);

xname = get(get(h.togglepanel,'SelectedObject'),'tag');
ynames = h.user.simloopsnames(~strcmp(h.user.simloopsnames,xname));

xind = find(strcmp(h.user.simloopsnames,xname));
x = h.user.sims.(xname);

%     end
%     [MaxVal, MaxLoc] = max(Likelihood(:));
%     MaxLoc = MaxLoc-1;
%     Index1 = floor(MaxLoc / (length(RhDecayFact) * length(RhTCon)));
%     Index2 = floor((MaxLoc - Index1 * (length(RhDecayFact) * length(RhTCon))) / length(RhDecayFact));
%     Index3 = MaxLoc - Index1 * (length(RhDecayFact) * length(RhTCon)) - Index2 * length(RhDecayFact);
%
%     fprintf(1, '%d %d %d\n', Index1+1, Index2+1, Index3+1);
%     fprintf(1, '%d %d %d\n', NumSteps(Index1+1), RhTCon(Index2+1), RhDecayFact(Index3+1));
%
% end
%
% ContourLevels = [-4:0.5:0] + log10(max(Likelihood(:)));
% figure(1);clf
%
% for step = 1:length(NumSteps)
%     subplot(1, length(NumSteps), step);
%     contour(RhTCon, RhDecayFact, log10(Likelihood(:, :, step)), sort(ContourLevels))
%     temp = strcat(num2str(NumSteps(step)), ' steps');
%     text(RhTCon(2), RhDecayFact(2), temp);
% end
% subplot(1, length(NumSteps), 1)
% ylabel('Rh decay factor');
% subplot(1, length(NumSteps), ceil(length(NumSteps)/2))
% xlabel('Rh* duration');

for n = 1:length(ynames)
    sp = axes('units','normalized',...
        'position',[hoff+(gridpnt(2)-1)*wdth+hbrd, voff+vbrd+(gridmax(1) - gridpnt(1))*hght, wdth-2*hbrd, hght-2*vbrd],...
        'parent',h.canvas,...
        'ButtonDownFcn',{@zoomIn,h},'tag',ynames{n});
    xlabel(sp,xname);
    ylabel(ynames{n});
    gridpnt = increasecnt(gridpnt,gridmax);

    yind = find(strcmp(h.user.simloopsnames,ynames{n}));
    y = h.user.sims.(ynames{n});
    
    % matrices are surfed as colums (y) vs rows (x)
    lkhd = zeros(length(y),length(x));
    lkhdkey = key;
    for i = 1:length(x)
        for j = 1:length(y)
            lkhdkey(xind) = i;
            lkhdkey(yind) = j;
            lkhd(j,i) = h.user.likelihoodmp.get(makeKey(lkhdkey));
        end
    end
    colormap hsv
    contour(sp,x,y,lkhd);
    
    xlabel(sp,xname);
    ylabel(sp,ynames{n});
    zlabel(sp,'Likelihood');

end

function zoomIn(sp,evnt,h)
delete(findobj('parent',h.canvas,'-not','tag',get(sp,'tag')));
hoff = .025; hbrd = .04;
voff = .03; vbrd = .05;
set(sp,'position',[hoff+hbrd, voff+vbrd, (1-hoff)-2*hbrd, (1-voff)-2*vbrd],...
    'ButtonDownFcn',{@plotLikelihoodSurfs,h});

function plotMeasureVs(hObject,h,feature,varargin)
delete(get(h.canvas,'children'));

xname = get(get(h.togglepanel,'SelectedObject'),'tag');
x = h.user.sims.(xname);
y = zeros(size(x));
key = h.user.pkey;
pnames = fieldnames(h.user.simulations.parameters)';
for i = 1:length(x)
    key(strcmp(pnames,xname)) = i;
    features = h.user.simulations.featuresByIndices(key);
    y(i) = features.(feature);
end
border = .1;
sp = axes('position',[border,border,1-2*border,1-2*border],'parent',h.canvas);
if ~isempty(varargin);
    line(x,ones(size(x))*varargin{1},...
        'parent',sp,'Linestyle','--','Color',.7*[1 1 1]);
end
line(x,y,...
    'parent',sp,...
    'Linestyle','--','Color',[0 0 1],...
    'Marker','o','MarkerFaceColor','none','MarkerEdgeColor',[0 0 1]);

xlabel(xname);
ylabel(feature);
a = h.user.simulations.getFeatureLims(feature);
axis tight
if strcmp(feature,'Likelihood')
    guidata(hObject,h);
    set(sp,'yscale','log')
    popxaxis(sp);
    ylims = log(get(sp,'ylim'));  dist = ylims(2)-ylims(1); 
    ylims = ylims + .1*dist*[-1 1]; set(sp,'ylim',exp(ylims));
    return
end
if ~isempty(varargin)
    ylims = [min([a,varargin{1}]),max([a,varargin{1}])]; set(sp,'ylim',ylims);
end
popaxis(sp);
guidata(hObject,h);

function ShowExamples(hObject,evnt,h)
h = guidata(hObject);
key = h.user.pkey;

% Check if the sim has been run
if ~isfield(h.user,'SimExample') || sum(key~=h.user.SimExample.pkey)
    sim = h.user.simulations.simByIndices(key);
    
    % add to simulation skeleton parameters that don't change
    sim.FreqCutoff = 20;
    sim.NumResponses = 200;
    sim.SamplingInterval = 1;         % ms
    sim.EpochPts = length(h.user.munorm);
    sim.NumPtsToPeak = 100;
    sim.Resp = h.user.munorm;
    sim.target = h.user.target;
    sim.SimFunc = h.user.simulationfunction;         % run rh or rhplustr function
    sim.storeExamples = 1;
    
    features = DoOneCascadeModelTA(sim);
    h.user.SimExample = features.SimulateCondition;
    h.user.SimExample.pkey = key;
end

guidata(hObject,h);
plotCheckedCharts(hObject,h)


function populateOutput(ot,eventdata,h)
h = guidata(ot);

features = h.user.simulations.featuresByIndices(h.user.pkey);
fn = fieldnames(features);
data = get(ot,'data');
rn = get(ot,'rowname');
r = 0;
for i = 1:length(rn)
    for j = 1:length(fn)
        if ~isempty(strfind(fn{j},rn{i}))
            data{i,2} = features.(fn{j});
            r = r+1;
        end
    end
end
slkey = h.user.slkey;
simloopsnames = h.user.simloopsnames;
for sln = 1:length(simloopsnames);
    data{r+sln,2} = h.user.simulations.parameters.(simloopsnames{sln})(slkey(sln));
    rn{r+sln} = h.user.simloopsnames{sln}(1:6);
end

set(ot,'data',data);
set(ot,'rowname',rn);
set(ot,'columnname',{'best','curr'});
guidata(ot,h);
    
% ------------------ Internal functions ------------------



% ------------------ Data plots checks ------------------
function plotCheckedCharts(hObject,h)
% only run if a simulation has been run and if SimConditions is plotted
h = guidata(hObject);
conditions = get(h.plotlist,'string');
if ~isfield(h,'user') ||...
        ~isfield(h.user,'simulationfunction') ||...
        ~strcmp(conditions{get(h.plotlist,'value')},'ShowExamples')||...
        ~isfield(h.user,'SimExample');
    return
end

% clear the canvas
delete(get(h.canvas,'children'));

checks = get(h.plotchecks,'children');
checked = zeros(2,3);
dataplot = 0; filterplot = 0;
traces = false;
for c = 1:length(checks);
    if get(checks(c),'value')
        switch get(checks(c),'tag')
            case 'rhstar', checked(1,1) = checks(c);
            case 'rhstar_eg', checked(2,1) = checks(c);
            case 'tstar', checked(1,2) = checks(c);
            case 'tstar_eg', checked(2,2) = checks(c);
            case 'current', checked(1,3) = checks(c);
            case 'current_eg', checked(2,3) = checks(c);
            case 'mu_and_var', dataplot = get(checks(c),'value');
            case 'filter', filterplot = get(checks(c),'value');
        end
    end
end

if strcmp(func2str(h.user.simulationfunction),'RhCascadeModelTA');
    checked(:,2) = false;
end

if strcmp(func2str(h.user.simulationfunction),'RhPlusTransducinCascadeModelTA')
    checked(:,1) = false;
end

if sum(sum(checked))
    traces = true;
end

border = .04;
if traces && (dataplot||filterplot)
    if (dataplot && filterplot)
        sp = axes('position',[.5+border,.5+border,.5-2*border,.5-2*border],'Parent',h.canvas);
        h = plotMuAndVar(sp,h);
        sp = axes('position',[.5+border,border,.5-2*border,.5-2*border],'Parent',h.canvas);
        h = plotFilter(sp,h);
    else
        sp = axes('position',[.5+border,border,.5-2*border,1-2*border],'parent',h.canvas);
        if dataplot
            h = plotMuAndVar(sp,h);
        elseif filterplot
            h = plotFilter(sp,h);
        end
    end
else
    if (dataplot && filterplot)
        sp = axes('position',[border,.5+border,1-2*border,.5-2*border],'Parent',h.canvas);
        h = plotMuAndVar(sp,h);
        sp = axes('position',[border,border,1-2*border,.5-2*border],'Parent',h.canvas);
        h = plotFilter(sp,h);
    elseif (dataplot || filterplot)
        sp = axes('position',[border,border,1-2*border,1-2*border],'parent',h.canvas);
        if dataplot
            h = plotMuAndVar(sp,h);
        elseif filterplot
            h = plotMuAndVar(sp,h);
        end
    end
end    
        
dataplot = (dataplot||filterplot);

w1 = 1/sum(traces+dataplot);
w2 = w1/sum(sum(checked)>0);
l2 = 1/(1+(sum(checked(2,:))>0));
topleft = 0;
bottomleft = topleft;
topbottom = 1-l2;
if checked(1,1)
    sp = axes('position',[topleft+border,topbottom+border,w2-2*border,l2-2*border],...
        'Parent',h.canvas);
    h = plotRhData(sp,h);
    topleft = w2;
end
if checked(2,1)
    sp = axes('position',[bottomleft+border,border,w2-2*border,l2-2*border],...
        'Parent',h.canvas);
    h = plotRhEGData(sp,h);
    topleft = w2;
end
bottomleft = topleft;

if checked(1,2)
    sp = axes('position',[topleft+border,topbottom+border,w2-2*border,l2-2*border],...
        'Parent',h.canvas);
    h = plotTData(sp,h);
    topleft = topleft+w2;
end
if checked(2,2)
    sp = axes('position',[bottomleft+border,border,w2-2*border,l2-2*border],...
        'Parent',h.canvas);
    h = plotTEGData(sp,h);
    if bottomleft==topleft
        topleft = topleft+w2;
    end
end
bottomleft = topleft;

if checked(1,3)
    sp = axes('position',[topleft+border,topbottom+border,w2-2*border,l2-2*border],...
        'Parent',h.canvas);
    h = plotIData(sp,h);
end
if checked(2,3)
    sp = axes('position',[bottomleft+border,border,w2-2*border,l2-2*border],...
        'Parent',h.canvas);
    h = plotIEGData(sp,h);
end


guidata(hObject,h)

function h = plotRhData(sp,h)
line(h.user.time,h.user.SimExample.MeanStochTimeCourse,'color',[0 0 0],'parent',sp);
axis tight, xlim([0 1000]), popaxis(sp);

function h = plotRhEGData(sp,h)
if isfield(h.user.SimExample,'RhExamples')
    line(h.user.time,h.user.SimExample.TExamples,'color',[0 0 0],'parent',sp);
end
axis tight, xlim([0 1000]), popaxis(sp);


function h = plotTData(sp,h)
line(h.user.time,h.user.SimExample.MeanStochTimeCourse,'color',[0 0 0],'parent',sp);
axis tight, xlim([0 1000]), popaxis(sp);


function h = plotTEGData(sp,h)
if isfield(h.user.SimExample,'TExamples')
    line(h.user.time,h.user.SimExample.TExamples,'color',[0 0 0],'parent',sp);
end
axis tight, xlim([0 1000]), popaxis(sp);

function h = plotIData(sp,h)
line(h.user.time,h.user.SimExample.AverageResponse/max(h.user.SimExample.AverageResponse),'color',[0 0 0],'Linestyle','-','parent',sp);
line(h.user.time,h.user.munorm,'color',[0 0 0],'Linestyle','--','parent',sp);
axis tight, xlim([0 1000]), popaxis(sp);

function h = plotIEGData(sp,h)
if isfield(h.user.SimExample,'IExamples')
    line(h.user.time,h.user.SimExample.IExamples,'color',[0 0 0],'parent',sp);
end
axis tight, xlim([0 1000]), popaxis(sp);

function h = plotFilter(sp,h)
%line(h.user.time,h.user.SimExample.Filter,'color',[0 0 0],'parent',sp);
line((1:length(h.user.SimExample.Filter)),h.user.SimExample.Filter,'color',[0 0 0],'parent',sp);
axis tight, popaxis(sp);

function h = plotMuAndVar(sp,h)
avesqsim = h.user.SimExample.AverageResponse;
[MaxVal, WTMeanTPeak] = max(avesqsim);
avesqsim = (avesqsim / MaxVal).^2;
avesqsim = avesqsim / 10;
varsim = h.user.SimExample.VarianceResponse / MaxVal^2;

line(h.user.time,avesqsim,'color',[0 0 1],'parent',sp);
line(h.user.time,varsim,'color',[0 0 0],'parent',sp);

[MaxVal, MeanTPeak] = max(h.user.munorm);

line(h.user.time,h.user.munorm.^2/ MaxVal^2 / 10,'color',[0 0 1],'Linestyle','--','parent',sp);
line((1:length(h.user.varnorm)),h.user.varnorm/ MaxVal^2,'color',[0 0 0],'Linestyle','--','parent',sp);

xlabel('normalized time');
ylabel('normalize amp');
ylim([-.05 .15])
axis tight, popaxis(sp);
% ------------------ Data plots checks ------------------




% ------------------ Create functions ------------------

function plotlist_CreateFcn(hObject, eventdata, h)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
plist = {'Plots'};
set(hObject,'string',plist);


% Wild type, GCAP KO, etc...
function condition_selector_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function repetitions_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ------------------ Create functions ------------------




% ------------------ Utilities ------------------
function cnt = increasecnt(cnt,maxcnt)
cnt(end) = cnt(end)+1;
if cnt(end) > maxcnt(end)
    if length(cnt) == 1;
        cnt = [];
        return
    end
    n = increasecnt(cnt(1:end-1),maxcnt(1:end-1));
    if isempty(n);
        cnt = [];
        return
    end
    cnt = [n,1];
end

function ms = structsAndMaps(sm)
if isstruct(sm)
    ms = java.util.HashMap;
    field = fieldnames(sm);
    for name = 1:length(field)
        if ~isnumeric(sm.(field{name}));
            continue
        end
        ms.put(field{name},sm.(field{name}));
    end
elseif strcmp(class(sm),'java.util.HashMap')
    keys = sm.keySet.toArray;
    for i = 1:length(keys)
        ms.(keys(i)) = sm.get(keys(i));
    end
end

function out = makeKey(in)
if isnumeric(in)
    out = '';
    for ind = 1:length(in)
        out = sprintf('%s,%d',out,in(ind));
    end
    out = out(2:end);
elseif ischar(in)
    out = regexp(in,'\,','split');
    out = str2double(out);
end                           

function h = makeKeysFromSliders(h,slkey)
slnames = h.user.simloopsnames;
pnames = fieldnames(h.user.simulations.parameters);

pkey = ones(1,length(pnames));
for sln = 1:length(slnames)
    pkey(ismember(pnames,slnames{sln})) = slkey(sln);
end
h.user.pkey = pkey;
h.user.slkey = slkey;

function popaxis(sp)
popxaxis(sp)
popyaxis(sp)

function popxaxis(sp)
xlims = get(sp,'xlim');  dist = xlims(2)-xlims(1); 
xlims = xlims + .1*dist*[-1 1]; set(sp,'xlim',xlims);

function popyaxis(sp)
ylims = get(sp,'ylim');  dist = ylims(2)-ylims(1); 
ylims = ylims + .1*dist*[-1 1]; set(sp,'ylim',ylims);


% --- Executes on button press in toIgor.
function toIgor_Callback(hObject, eventdata, handles)
figData = guidata(hObject);
plotCanvas = figData.canvas;
ch = get(plotCanvas,'children');
moving_ch = [];

for i=1:length(ch)
    if ~isprop(ch(i),'style') || ~strfind(get(ch(i),'style'), 'button')
        moving_ch = [moving_ch ch(i)];
    end
end

for c = 1:length(moving_ch)
    prompt{c} = sprintf('Axis %d',c);
    def{c} = get(moving_ch(c),'tag');
end
dlg_title = 'Names for Axes';

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answer = inputdlg(prompt,dlg_title,1,def,options);

for c = 1:length(moving_ch)
    makeAxisStructLocal(moving_ch(c),answer{c});
end


% --- Executes on button press in rhstar.
function rhstar_Callback(hObject, eventdata, handles)
% hObject    handle to rhstar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rhstar
