function compare_EWCSPH_streamfunction_contours(folderList)
%COMPARE_EWCSPH_INTERACTIVE Responsive, synchronized EWCSPH comparison.
%   compare_EWCSPH_matched_colours() opens multi-folder chooser (Ctrl/Shift+click).
%   compare_EWCSPH_matched_colours({"C:\\Case1","C:\\Case2"}) uses supplied folders.
%   MATLAB R2021b+ recommended (clim, shared axes controls).
%
% All panels use IDENTICAL x/D, y/D, colour limits, and axes aspect.
% The timeline is the FIRST case's frames in the shared time interval.
% Other cases display their nearest saved timestep at each selected time.
%
%% USER SETTINGS
plotVariable = "vorticity"; % vorticity | velocity | u | v | rho | pressure
viewX_D = [-10 10];          % identical dimensionless x-view in every panel
viewY_D = [-1.3 1.3];      % identical dimensionless y-view in every panel
manualCLim = []; % [] = shared 5th-95th percentile auto scale; or enter e.g. [-4.5 4.5]
colorPercentiles = [5 95]; % Match original single-case MATLAB viewer
colorSampleFrames = 5; % evenly spaced real snapshots per case used for auto limits
playbackFPS = 60;   % Target playback speed; disk read/streamline tracing may limit actual FPS
playbackFrameStride = 1; % Fast preview: advance 0.1 s per rendered frame; slider/video still use 0.1 s
streamGridSize = 300;       % Cartesian streamline integration resolution
streamSeedSpacingD = 0.15; % Transverse spacing of incoming streamlines (in D)
wakeSeedX_D = [0.85 1.10 1.40]; % A few recirculation seeds behind cylinder
wakeSeedY_D = [0.08 -0.08 0.05];
streamStep = 0.25;         % streamline integrator step (grid cells)
streamMaxVertices = 1500;
% Streamline rendering options:
% "streamfunction" uses YOUR EWCSPH velocity snapshots, reconstructing
% psi via a polar Poisson solution; "particle" retains original stream2.
streamlineMethod = "streamfunction"; % "streamfunction" | "particle"
streamfunctionSpacingD = 0.015;  % Delta psi/(U_inf*D); smaller = denser, GLOBAL
streamfunctionNearWallContours = true;
% Additional low-|psi| levels apply THROUGHOUT the domain, including wake.
streamfunctionNearWallLevelsD = [0.0005 0.001 0.002 0.004 0.008 0.015 0.03 0.05];
streamfunctionVorticitySource = "velocity"; % "velocity" | "solver"

initialShowStreamlines = true;
showCylinder = true;
frameCacheSize = 60;   % Cache recently viewed frames per case (reduces slider/loop disk reads)
startWithFloatingController = true;
startWithCaseList = true; % Show a separate, hideable folder key (not part of exports)
comparisonTimeStep = 0.1; % s, timeline spacing; actual frames are never invented
maximumFrames = 1201; % cap on shared animation frames

%% SELECT MULTIPLE SIMULATION FOLDERS IN A SINGLE DIALOG
% On desktop MATLAB the Java directory chooser supports Shift/Ctrl+click.
% Fallback: select folders individually if Java GUI is unavailable.
if nargin<1 || isempty(folderList)
    folderList = chooseSimulationFolders();
end
folderList=string(folderList(:));
folderList=folderList(strlength(folderList)>0);
if isempty(folderList), disp('No folders selected.'); return; end
nCases=numel(folderList);

% User-requested orientation refers to the actual FIELD OF VIEW:
% landscape => one column, portrait => one row, square => grid.
viewAspect=diff(viewX_D)/diff(viewY_D);
if viewAspect>=1.3
    nRows=nCases; nCols=1;
elseif viewAspect<=1/1.3
    nRows=1; nCols=nCases;
else
    nCols=ceil(sqrt(nCases)); nRows=ceil(nCases/nCols);
end

%% INDEX DATASETS; READ FIRST FRAME OF EACH ONLY
S=cell(nCases,1);
for ic=1:nCases
    folder=folderList(ic);
    metadataFile=fullfile(folder,'SimulationConstants.csv');
    if ~isfile(metadataFile), error('Missing %s',metadataFile); end
    meta=readcell(metadataFile,'Delimiter',',');
    Ntheta=round(metaNumber(meta,'NTheta'));
    D=metaNumber(meta,'Dcylinder');
    Uinf=metaNumber(meta,'inletVelocity');
    rho0=metaNumber(meta,'rho0');
    Re=metaNumber(meta,'Reynolds');
    if Uinf<=0 || rho0<=0, error('Invalid inlet velocity / rho0 in %s',folder); end
    if Ntheta<3 || D<=0, error('Invalid NTheta or Dcylinder in %s',folder); end

    f=dir(fullfile(folder,'EWCSPH_hdp_*_t_*.csv'));
    if isempty(f), error('No EWCSPH timestep files in %s',folder); end
    tt=nan(numel(f),1); hh=nan(numel(f),1);
    numberPattern='[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?';
    for k=1:numel(f)
        tok=regexp(f(k).name,['EWCSPH_hdp_(' numberPattern ')_t_(' numberPattern ')\.csv$'],'tokens','once');
        if ~isempty(tok), hh(k)=str2double(tok{1}); tt(k)=str2double(tok{2}); end
    end
    good=isfinite(tt)&isfinite(hh);
    f=f(good); tt=tt(good); hh=hh(good);
    if isempty(f), error('Cannot parse saved timestep names in %s',folder); end
    allH=unique(hh);
    if numel(allH)>1
        [choice,ok]=listdlg('PromptString',sprintf('Select h/dp for folder %d',ic), ...
            'SelectionMode','single','ListString',cellstr(compose('%.12g',allH)));
        if ~ok, return; end
        h=allH(choice);
    else
        h=allH(1);
    end
    take=abs(hh-h)<=1e-10*max(1,abs(h));
    f=f(take); tt=tt(take);
    [tt,ord]=sort(tt); f=f(ord);
    % Avoid duplicate saved times in navigation.
    [tt,keep]=unique(tt,'stable'); f=f(keep);

    T=readParticleCSV(fullfile(f(1).folder,f(1).name));
    typ=lower(strtrim(string(T.Type)));
    fluidRows=find(typ=="fluid"); wallRows=find(typ=="boundary");
    if isempty(fluidRows)||mod(numel(fluidRows),Ntheta)~=0
        error('Fluid particle order/count is inconsistent with NTheta in %s',folder);
    end
    nr=numel(fluidRows)/Ntheta;
    if nr<2, error('At least two fluid rings are required in %s',folder); end
    ringRows=reshape(fluidRows,Ntheta,nr);
    surfaceRows=[ringRows;ringRows(1,:)];
    if ~isempty(wallRows)
        xc=mean(T.x(wallRows)); yc=mean(T.y(wallRows));
    else
        xc=mean(T.x(ringRows(:,1))); yc=mean(T.y(ringRows(:,1)));
    end
    xD=(T.x-xc)/D; yD=(T.y-yc)/D;
    [~,nm]=fileparts(char(folder));
    theta=mod(atan2(yD(ringRows(:,1)),xD(ringRows(:,1))),2*pi);
    [theta,angleOrder]=sort(theta);
    radii=mean(hypot(xD(ringRows),yD(ringRows)),1);
    if any(diff(radii)<=0), error('Fluid rings do not increase in radius: %s',folder); end
    d=struct();
    d.folder=folder;d.name=nm;d.Re=Re;d.D=D;d.Uinf=Uinf;d.rho0=rho0;d.times=tt;d.files=f;
    d.nTheta=Ntheta;d.fluidRows=fluidRows;d.surfaceRows=surfaceRows;
    d.ringRows=ringRows;d.wallRows=wallRows;d.xD=xD;d.yD=yD;
    d.surfaceX=reshape(xD(surfaceRows(:)),size(surfaceRows));
    d.surfaceY=reshape(yD(surfaceRows(:)),size(surfaceRows));
    d.wallX=xD(wallRows);d.wallY=yD(wallRows);
    d.theta=theta;d.angleOrder=angleOrder;d.radii=radii;
    % Precompute polar FFT/tridiagonal coefficients ONCE per case.
    d.psiSolver=[];
    if ~isempty(wallRows)
        wallRadiusD=mean(hypot(xD(wallRows),yD(wallRows)));
        d.psiSolver=ewcsphCreatePolarPsiSolver(theta,wallRadiusD,radii);
    end

    d.currentFile=1;d.table=T;
    d.frameCacheIndices=1; d.frameCacheTables={T};
    S{ic}=d;
    fprintf('Case %d: %s | Re=%.5g | h/dp=%.5g | %d frames, t=[%.5g %.5g]\n', ...
        ic,nm,Re,h,numel(tt),tt(1),tt(end));
end

% Shared interval: don't extrapolate outside a case's saved time range.
tStart=max(cellfun(@(q)q.times(1),S));
tEnd=min(cellfun(@(q)q.times(end),S));
if tEnd<tStart, error('No overlapping simulation time interval across the selected cases.'); end
% A regular 0.1-s comparison timeline, capped at 1201 frames.
% Note: each case still displays its closest REAL CSV frame. A 0.2-s
% saved cadence cannot be turned into 0.1-s physical solution data.
firstTick=ceil((tStart-1e-9)/comparisonTimeStep);
lastTick=floor((tEnd+1e-9)/comparisonTimeStep);
if lastTick>=firstTick
    tickIndices=(firstTick:min(lastTick,firstTick+maximumFrames-1))';
    masterTimes=tickIndices*comparisonTimeStep;
else
    masterTimes=tStart;
end
if isempty(masterTimes), error('No saved frames within the shared time interval.'); end
nFrames=numel(masterTimes);

%% SHARED NONDIMENSIONAL COLOUR LIMITS
[fieldLabel,~]=fieldMetadata(plotVariable);
if isempty(manualCLim)
    manualCLim=initialSharedLimits();
end

%% RESPONSIVE FIGURE AND EQUAL-SIZE PHYSICAL PLOT REGIONS
% A common, fixed X/Y span plus equal pixel plot boxes makes each cylinder
% visually identical in diameter, even with different MATLAB tile shapes.
% True 16:9 figure canvas for PowerPoint; scale to available screen
% while preserving the exact ratio (rather than setting width/height separately).
screen=get(groot,'ScreenSize');
fitScale=min([1, (screen(3)-90)/1600, (screen(4)-110)/900]);
fitScale=max(0.45,fitScale);
figW=max(720,16*floor((1600*fitScale)/16));
figH=figW*9/16;
fig=figure('Color','w','Name','EWCSPH | 16:9 fast comparison + MP4 export', ...
    'NumberTitle','off','Units','pixels','Position',[40 65 figW figH]);
mainPanel=uipanel(fig,'Units','pixels','BorderType','none','BackgroundColor','w');
AX=gobjects(nCases,1);SURF=gobjects(nCases,1);STREAM=gobjects(nCases,1);
WALL=gobjects(nCases,1);TITLE=gobjects(nCases,1);
for ic=1:nCases
    d=S{ic};
    AX(ic)=axes('Parent',mainPanel,'Units','pixels','Color','w');
    firstC=fieldOfND(d.table,plotVariable,d);
    firstZ=reshape(firstC(d.surfaceRows(:)),size(d.surfaceRows));
    SURF(ic)=surf(AX(ic),d.surfaceX,d.surfaceY,zeros(size(firstZ)),firstZ, ...
        'FaceColor','interp','EdgeColor','none','HitTest','off');
    hold(AX(ic),'on');
    WALL(ic)=scatter(AX(ic),d.wallX,d.wallY,18,[0.12 0.12 0.12],'filled', ...
        'HitTest','off');
    WALL(ic).Visible=onoff(showCylinder);
    STREAM(ic)=line(AX(ic),nan,nan,nan,'Color',[0 0 0], ...
        'LineWidth',0.65,'Visible','off','HitTest','off');
    view(AX(ic),2);
    set(AX(ic),'SortMethod','childorder', 'XLim',viewX_D,'YLim',viewY_D, ...
        'XLimMode','manual','YLimMode','manual', ...
        'DataAspectRatio',[1 1 1],'DataAspectRatioMode','manual');
    clim(AX(ic),manualCLim);
    grid(AX(ic),'off');box(AX(ic),'on');
    xlabel(AX(ic),'x/D');ylabel(AX(ic),'y/D');
    TITLE(ic)=title(AX(ic),'','Interpreter','none','FontSize',9);
end
% Link all axes: zoom/pan in one applies to every case.
if nCases>1, linkaxes(AX,'xy');end
setComparisonColormap();
cb=colorbar(AX(end),'eastoutside');
cb.Label.String=fieldLabel;
cb.Units='pixels';
if initialShowStreamlines
    set(STREAM,'Visible','on');
end
% Keep field surface and colourbar visible beneath the streamline overlay.
set(SURF,'Visible','on');
cb.Visible='on';

%% SHARED INTERACTIVE CONTROLS
controls=uipanel(fig,'Units','pixels','BorderType','none', ...
    'BackgroundColor',[0.96 0.96 0.96]);
makeBtn('|<',[.012 .61 .045 .27],@(~,~)setFrame(1,false));
makeBtn('<',[.062 .61 .045 .27],@(~,~)setFrame(currentFrame-1,false));
playBtn=makeBtn('Play',[.112 .61 .070 .27],@togglePlay);
makeBtn('>',[.187 .61 .045 .27],@(~,~)setFrame(currentFrame+1,false));
makeBtn('>|',[.237 .61 .045 .27],@(~,~)setFrame(nFrames,false));
frameSlider=uicontrol(controls,'Style','slider','Units','normalized', ...
    'Position',[.300 .64 .237 .19],'Min',1,'Max',max(2,nFrames), ...
    'Value',1,'SliderStep',[1/max(1,nFrames-1) min(1,10/max(1,nFrames-1))], ...
    'Callback',@(~,~)sliderCommitted());
if nFrames==1,frameSlider.Enable='off';end
% Standard uicontrol Callback normally fires only on release. Listen to
% ContinuousValueChange as well so dragging previews the actual flow.
sliderListener=[];
try
    sliderListener=addlistener(frameSlider,'ContinuousValueChange', ...
        @(~,~)sliderDragging());
catch ME
    warning('Continuous slider updates unavailable: %s',ME.message);
end
videoButton=makeBtn('Save MP4',[.544 .61 .093 .27],@saveAnimation);
videoButton.TooltipString='Export all comparison time steps (0.1 s) as a 16:9 MP4';
streamCheck=makeCheck('Streamlines',[.645 .61 .105 .27],initialShowStreamlines,@toggleStreamlines);
makeBtn('Zoom', [.754 .61 .057 .27],@toggleZoom);
makeBtn('Pan',  [.815 .61 .055 .27],@togglePan);
makeBtn('Reset view',[.875 .61 .109 .27],@resetView);

% Second row: field, colour range, view limits, layout, and FPS.
makeLabel('Field',[.010 .30 .045 .18]);
fieldItems={'Vorticity','Velocity magnitude','u velocity','v velocity','Density','Pressure'};
fieldKeys=["vorticity","velocity","u","v","rho","pressure"];
fieldMenu=uicontrol(controls,'Style','popupmenu','Units','normalized', ...
    'Position',[.056 .24 .133 .25],'String',fieldItems,'Value',find(fieldKeys==plotVariable,1), ...
    'Callback',@changeField);
makeLabel('Colour',[.195 .30 .046 .18]);
cMinEdit=makeEdit(manualCLim(1),[.244 .24 .060 .25]);
cMaxEdit=makeEdit(manualCLim(2),[.308 .24 .060 .25]);
makeLabel('x/D',[.372 .30 .033 .18]);
xMinEdit=makeEdit(viewX_D(1),[.408 .24 .055 .25]);
xMaxEdit=makeEdit(viewX_D(2),[.466 .24 .055 .25]);
makeLabel('y/D',[.524 .30 .033 .18]);
yMinEdit=makeEdit(viewY_D(1),[.560 .24 .055 .25]);
yMaxEdit=makeEdit(viewY_D(2),[.618 .24 .055 .25]);
makeBtn('Apply',[.679 .24 .058 .25],@applySettings);
makeLabel('Layout',[.742 .30 .048 .18]);
layoutMenu=uicontrol(controls,'Style','popupmenu','Units','normalized', ...
    'Position',[.792 .24 .090 .25], ...
    'String',{'Auto','Column','Row','Grid'},'Value',1,'Callback',@(~,~)layoutPanels());
makeLabel('FPS',[.886 .30 .034 .18]);
fpsEdit=makeEdit(playbackFPS,[.923 .24 .060 .25]);

timeLabel=uicontrol(controls,'Style','text','Units','normalized', ...
    'Position',[.02 .02 .74 .14],'String','', ...
    'HorizontalAlignment','left','BackgroundColor',controls.BackgroundColor);
makeLabel('Go to t(s):',[.765 .025 .080 .14]);
timeEdit=uicontrol(controls,'Style','edit','Units','normalized', ...
    'Position',[.850 .012 .091 .18],'String',sprintf('%.3f',masterTimes(1)), ...
    'BackgroundColor','w','Callback',@jumpToTime);

currentFrame=1;playing=false;showStreams=initialShowStreamlines;
lastDragClock=tic; isUpdatingFrame=false;
controllerFig=gobjects(1);controllerFloating=false;
caseListFig=gobjects(1);
lastStreamFrame=-1;
% Export commands and keyboard navigation, in addition to on-screen controls.
exportMenu=uimenu(fig,'Label','Export');
uimenu(exportMenu,'Label','Save comparison PNG','Callback',@saveSnapshot);
uimenu(exportMenu,'Label','Save synchronized MP4','Callback',@saveAnimation);
controllerMenu=uimenu(fig,'Label','Controller');
uimenu(controllerMenu,'Label','Float / move controller','Callback',@(~,~)floatController());
uimenu(controllerMenu,'Label','Dock controller at bottom','Callback',@(~,~)dockController());
caseMenu=uimenu(fig,'Label','Cases');
uimenu(caseMenu,'Label','Show / hide folder key','Callback',@(~,~)toggleCaseList());
set(fig,'WindowKeyPressFcn',@keyControl);
set(fig,'CloseRequestFcn',@closeMainFigure);
set(fig,'SizeChangedFcn',@(~,~)layoutPanels());
if startWithFloatingController, floatController();end
layoutPanels();
setFrame(1,false);
if startWithCaseList,showCaseList();end

%% NESTED FUNCTIONS
    function toggleCaseList()
        if isgraphics(caseListFig)
            hideCaseList();
        else
            showCaseList();
        end
    end

    function showCaseList()
        if isgraphics(caseListFig)
            figure(caseListFig);
            return;
        end
        % Standalone reference window: does not occupy plot space and is
        % deliberately absent from both exportgraphics(fig) and getframe(fig).
        mainPos=get(fig,'Position');
        keyW=480;
        keyH=min(440,max(125,62+29*nCases));
        keyX=max(30,mainPos(1)+mainPos(3)-keyW);
        keyY=max(40,mainPos(2)+mainPos(4)-keyH-30);
        caseListFig=figure('Name','EWCSPH | Case folder key', ...
            'NumberTitle','off','MenuBar','none','ToolBar','none', ...
            'Color',[0.98 0.98 0.98],'Units','pixels', ...
            'Position',[keyX keyY keyW keyH], ...
            'Resize','on','CloseRequestFcn',@(~,~)hideCaseList());
        uicontrol(caseListFig,'Style','text','Units','normalized', ...
            'Position',[.025 .84 .95 .12], ...
            'BackgroundColor',[0.98 0.98 0.98], ...
            'HorizontalAlignment','left','FontWeight','bold', ...
            'String','Panel order: top to bottom, left to right');
        caseRows=cell(nCases,2);
        for ci=1:nCases
            caseRows{ci,1}=sprintf('Panel %d',ci);
            caseRows{ci,2}=char(S{ci}.folder);
        end
        uitable(caseListFig,'Units','normalized', ...
            'Position',[.025 .07 .95 .76], ...
            'Data',caseRows,'ColumnName',{'Plot','Simulation folder (full path)'}, ...
            'ColumnEditable',[false false], 'RowName',[], ...
            'ColumnWidth',{75 355});
    end

    function hideCaseList()
        if isgraphics(caseListFig)
            h=caseListFig;
            caseListFig=gobjects(1);
            set(h,'CloseRequestFcn','');
            delete(h);
        end
    end

    function keyControl(~,event)
        switch event.Key
            case {'rightarrow','d'}
                setFrame(currentFrame+1,false);
            case {'leftarrow','a'}
                setFrame(currentFrame-1,false);
            case 'home'
                setFrame(1,false);
            case 'end'
                setFrame(nFrames,false);
        end
    end
    function saveSnapshot(~,~)
        [name,folder]=uiputfile({'*.png','PNG image'},'Save comparison snapshot', ...
            sprintf('EWCSPH_comparison_t_%g.png',masterTimes(currentFrame)));
        if isequal(name,0),return;end
        wasFloating=controllerFloating;
        if ~wasFloating,floatController();end
        oldPosition=get(fig,'Position');
        enforceExportAspect();
        try
            drawnow;
            exportgraphics(fig,fullfile(folder,name),'Resolution',180);
        catch err
            errordlg(sprintf('PNG export failed: %s',err.message),'Export');
        end
        if isgraphics(fig)
            set(fig,'Position',oldPosition);
            if ~wasFloating,dockController();end
        end
    end
    function saveAnimation(~,~)
        [name,folder]=uiputfile({'*.mp4','MP4 video'},'Save full 0.1-s comparison animation (16:9)', ...
            'EWCSPH_comparison.mp4');
        if isequal(name,0),return;end
        previousFrame=currentFrame;
        previousPlay=playing;
        wasFloating=controllerFloating;
        if ~wasFloating,floatController();end % keep MP4 plot-only, true 16:9
        oldPosition=get(fig,'Position');
        enforceExportAspect();
        playing=false;
        fps=str2double(fpsEdit.String);
        if ~isfinite(fps)||fps<=0,fps=playbackFPS;end
        writer=[];
        try
            writer=VideoWriter(fullfile(folder,name),'MPEG-4');
            writer.FrameRate=10;writer.Quality=90;
            open(writer);
            for f=1:nFrames
                if ~isgraphics(fig),break;end
                setFrame(f,false);
                drawnow;
                img=getframe(fig);
                writeVideo(writer,img);
            end
            close(writer);
            writer=[];
        catch err
            if ~isempty(writer)
                try,close(writer);catch,end
            end
            errordlg(sprintf('Video export failed: %s',err.message),'Export');
        end
        if isgraphics(fig)
            set(fig,'Position',oldPosition);
            setFrame(previousFrame,false);
            playing=previousPlay && false;
            if ~wasFloating,dockController();end
            if isgraphics(playBtn),playBtn.String='Play';end
        end
    end
    function h=makeCheck(label,pos,val,cbfun)
        h=uicontrol(controls,'Style','checkbox','Units','normalized', ...
            'Position',pos,'String',label,'Value',double(val), ...
            'BackgroundColor',controls.BackgroundColor,'Callback',cbfun);
    end
    function h=makeEdit(v,pos)
        h=uicontrol(controls,'Style','edit','Units','normalized', ...
            'Position',pos,'String',sprintf('%.6g',v),'BackgroundColor','w');
    end
    function h=makeLabel(v,pos)
        h=uicontrol(controls,'Style','text','Units','normalized', ...
            'Position',pos,'String',v,'BackgroundColor',controls.BackgroundColor);
    end
    function s=onoff(tf)
        if tf,s='on';else,s='off';end
    end
    function setComparisonColormap()
        % Matches EWCSPH_MATLAB_ViscousSinglePlot_v2_streamline_repaired.
        if plotVariable == "vorticity"
            anchors = [0.10 0.25 0.80; 1 1 1; 0.80 0.10 0.10];
            cmap = interp1([0 0.5 1],anchors,linspace(0,1,257));
        else
            cmap = turbo(256);
        end
        colormap(fig,cmap);
    end

    function limits=initialSharedLimits()
        % Common 5th/95th percentiles over a representative selection of
        % real saved frames from EVERY case (bounded I/O and memory).
        values=[];
        for j=1:nCases
            d=S{j};
            sampleIndices=unique(round(linspace(1,numel(d.times), ...
                min(colorSampleFrames,numel(d.times)))));
            for si=sampleIndices(:)'
                if si==d.currentFile
                    sampleT=d.table;
                else
                    sampleT=readParticleCSV(fullfile(d.files(si).folder,d.files(si).name));
                end
                v=fieldOfND(sampleT,plotVariable,d);
                v=v(d.fluidRows);
                v=v(isfinite(v));
                % Subsample evenly if very large, to limit colour-limit memory.
                if numel(v)>100000
                    v=v(round(linspace(1,numel(v),100000)));
                end
                values=[values;v]; %#ok<AGROW>
            end
        end
        if isempty(values),error('No finite colour data in selected cases.');end
        values=sort(values);
        n=numel(values);
        lo=max(1,min(n,ceil(colorPercentiles(1)/100*n)));
        hi=max(1,min(n,ceil(colorPercentiles(2)/100*n)));
        limits=[values(lo) values(hi)];
        if plotVariable == "vorticity"
            m=max(abs(limits));
            if m==0 || ~isfinite(m),m=1;end
            limits=[-m m];
        elseif plotVariable == "velocity"
            limits(1)=0;
        end
        if limits(1)>=limits(2)
            extra=max(1,abs(limits(1)))*1e-6;
            limits=limits+[-extra extra];
        end
    end
    function enforceExportAspect()
        % Always produce a 16:9 plot-only canvas even if the user has
        % manually resized the comparison figure to another aspect ratio.
        if ~isgraphics(fig),return;end
        p=get(fig,'Position');
        w=max(720,16*floor(p(3)/16));
        h=w*9/16;
        set(fig,'Position',[p(1) p(2) w h]);
        layoutPanels();
    end
    function floatController()
        % A uipanel is tied to its parent MATLAB figure. A separate figure
        % gives it a genuine OS title bar, making the controller draggable.
        if controllerFloating && isgraphics(controllerFig)
            figure(controllerFig); return;
        end
        controllerFig=figure('Color',[.96 .96 .96], ...
            'Name','EWCSPH | Animation controller (drag this window)', ...
            'NumberTitle','off','MenuBar','none','ToolBar','none', ...
            'Units','pixels','Position',[80 80 1130 150], ...
            'Resize','on','CloseRequestFcn',@(~,~)dockController());
        cm=uimenu(controllerFig,'Label','Controller');
        uimenu(cm,'Label','Dock to plot window','Callback',@(~,~)dockController());
        set(controls,'Parent',controllerFig,'Units','normalized', ...
            'Position',[0 0 1 1]);
        controllerFloating=true;
        set(controllerFig,'SizeChangedFcn',@(~,~)resizeFloatingControls());
        set(controllerFig,'WindowKeyPressFcn',@keyControl);
        set(controllerFig,'WindowButtonUpFcn',@(~,~)sliderCommitted());
        layoutPanels();
    end
    function resizeFloatingControls()
        if controllerFloating && isgraphics(controllerFig) && isgraphics(controls)
            set(controls,'Units','normalized','Position',[0 0 1 1]);
        end
    end
    function dockController()
        if ~controllerFloating,return;end
        floatingHandle=controllerFig;
        controllerFloating=false;
        controllerFig=gobjects(1);
        if isgraphics(fig) && isgraphics(controls)
            set(controls,'Parent',fig,'Units','pixels');
        end
        if isgraphics(floatingHandle)
            set(floatingHandle,'CloseRequestFcn','');
            delete(floatingHandle);
        end
        if isgraphics(fig),layoutPanels();end
    end
    function closeMainFigure(~,~)
        playing=false;
        hideCaseList();
        if isgraphics(controllerFig)
            floatingHandle=controllerFig;
            set(floatingHandle,'CloseRequestFcn','');
            delete(floatingHandle);
        end
        if isgraphics(fig),delete(fig);end
    end
    function layoutPanels()
        if ~isgraphics(fig),return;end
        pos=get(fig,'Position');fw=pos(3);fh=pos(4);
        if controllerFloating && isgraphics(controllerFig)
            cp=get(controllerFig,'Position');
            set(controls,'Units','pixels','Position',[0 0 cp(3) cp(4)]);
            controlsH=0;
        else
            controlsH=125;
            set(controls,'Units','pixels','Position',[0 0 fw controlsH]);
        end
        set(mainPanel,'Position',[0 controlsH fw max(100,fh-controlsH)]);
        usableH=max(100,fh-controlsH);
        % Mode Auto follows the plotted flow-region aspect, not desktop size.
        switch layoutMenu.Value
            case 1
                aspect=diff(viewX_D)/diff(viewY_D);
                if aspect>=1.3,nr=nCases;nc=1;
                elseif aspect<=1/1.3,nr=1;nc=nCases;
                else,nc=ceil(sqrt(nCases));nr=ceil(nCases/nc);end
            case 2,nr=nCases;nc=1;
            case 3,nr=1;nc=nCases;
            otherwise,nc=ceil(sqrt(nCases));nr=ceil(nCases/nc);
        end
        % Reserve a dedicated strip for the SHARED colourbar.
        cbW=58;rightMargin=cbW+20;leftMargin=12;
        tileW=max(100,(fw-leftMargin-rightMargin)/nc);
        tileH=max(65,(usableH-12)/nr);
        wantedAspect=diff(viewX_D)/diff(viewY_D);
        % Identical plot-box width/height in every tile, centred in each tile.
        leftPad=51;rightPad=12;topPad=30;bottomPad=44;
        fitW=max(10,tileW-leftPad-rightPad);
        fitH=max(10,tileH-topPad-bottomPad);
        boxW=min(fitW,fitH*wantedAspect);
        boxH=boxW/wantedAspect;
        for j=1:nCases
            r=floor((j-1)/nc);c=mod(j-1,nc);
            tileX=leftMargin+c*tileW;
            tileTop=usableH-6-r*tileH;
            x=tileX+leftPad+(fitW-boxW)/2;
            y=tileTop-topPad-(fitH-boxH)/2-boxH;
            set(AX(j),'Units','pixels','Position',[x y boxW boxH]);
            set(AX(j),'DataAspectRatio',[1 1 1], ...
                'DataAspectRatioMode','manual','PlotBoxAspectRatioMode','auto');
        end
        set(cb,'Units','pixels','Position',[fw-53 max(42,usableH*.15) 17 max(55,usableH*.70)]);
    end
    function changeField(~,~)
        plotVariable=fieldKeys(fieldMenu.Value);
        [fieldLabel,~]=fieldMetadata(plotVariable);
        manualCLim=initialSharedLimits();
        cMinEdit.String=sprintf('%.6g',manualCLim(1));
        cMaxEdit.String=sprintf('%.6g',manualCLim(2));
        cb.Label.String=fieldLabel;
        setComparisonColormap();
        for j=1:nCases
            d=S{j};z=fieldOfND(d.table,plotVariable,d);
            set(SURF(j),'CData',reshape(z(d.surfaceRows(:)),size(d.surfaceRows)));
            clim(AX(j),manualCLim);
        end
    end
    function applySettings(~,~)
        vals=[str2double(xMinEdit.String) str2double(xMaxEdit.String) ...
              str2double(yMinEdit.String) str2double(yMaxEdit.String) ...
              str2double(cMinEdit.String) str2double(cMaxEdit.String)];
        if any(~isfinite(vals)) || vals(1)>=vals(2) || vals(3)>=vals(4) || vals(5)>=vals(6)
            errordlg('Enter finite increasing min/max values.','Invalid view');return;
        end
        viewX_D=vals(1:2);viewY_D=vals(3:4);manualCLim=vals(5:6);
        set(AX,'XLim',viewX_D,'YLim',viewY_D);
        for j=1:nCases,clim(AX(j),manualCLim);end
        layoutPanels();
    end
    function resetView(~,~)
        set(AX,'XLim',viewX_D,'YLim',viewY_D);
        zoom(fig,'off');pan(fig,'off');
    end
    function toggleZoom(~,~)
        pan(fig,'off');z=zoom(fig);
        if strcmp(z.Enable,'on'),z.Enable='off';else,z.Enable='on';end
    end
    function togglePan(~,~)
        zoom(fig,'off');p=pan(fig);
        if strcmp(p.Enable,'on'),p.Enable='off';else,p.Enable='on';end
    end
    function h=makeBtn(label,pos,cbfun)
        h=uicontrol(controls,'Style','pushbutton','Units','normalized', ...
            'Position',pos,'String',label,'Callback',cbfun);
    end

    function sliderDragging()
        if isUpdatingFrame || ~isgraphics(frameSlider),return;end
        k=max(1,min(nFrames,round(frameSlider.Value)));
        % The time readout updates immediately; throttle expensive CSV reads.
        if isgraphics(timeLabel)
            timeLabel.String=sprintf('Seeking: t = %.3f s | frame %d/%d', ...
                masterTimes(k),k,nFrames);
        end
        if toc(lastDragClock)>0.10 && k~=currentFrame
            lastDragClock=tic;
            setFrame(k,false);
        end
    end
    function sliderCommitted()
        if isUpdatingFrame || ~isgraphics(frameSlider),return;end
        k=max(1,min(nFrames,round(frameSlider.Value)));
        if k~=currentFrame,setFrame(k,false);end
    end
    function jumpToTime(~,~)
        t=str2double(timeEdit.String);
        if ~isfinite(t),timeEdit.String=sprintf('%.3f',masterTimes(currentFrame));return;end
        [~,k]=min(abs(masterTimes-t));
        setFrame(k,false);
    end

    function togglePlay(~,~)
        if playing
            playing=false;set(playBtn,'String','Play');return;
        end
        playing=true;set(playBtn,'String','Pause');
        if currentFrame>=nFrames, setFrame(1,false); end
        while playing && isgraphics(fig) && currentFrame<nFrames
            before=tic;
            setFrame(min(nFrames,currentFrame+playbackFrameStride),true);
            drawnow limitrate; % permits responsive Play/Pause and other callbacks
            fps=str2double(fpsEdit.String);
            if ~isfinite(fps)||fps<=0, fps=playbackFPS; end
            pause(max(0,1/fps-toc(before)));
        end
        playing=false;
        if isgraphics(playBtn), set(playBtn,'String','Play'); end
    end

    function toggleStreamlines(~,~)
        showStreams=logical(streamCheck.Value);
        for j=1:nCases
            if isgraphics(STREAM(j))
                STREAM(j).Visible=onoff(showStreams);
            end
            % Streamlines are an overlay; toggling them must not hide the field.
            if isgraphics(SURF(j)), SURF(j).Visible='on'; end
        end
        if isgraphics(cb), cb.Visible='on'; end
        if showStreams, lastStreamFrame=-1;setFrame(currentFrame,false); end
    end

    function setFrame(frameNumber,isPlayback)
        if ~isgraphics(fig), return; end
        isUpdatingFrame=true;
        currentFrame=max(1,min(nFrames,round(frameNumber)));
        target=masterTimes(currentFrame);
        % Trace streamlines from the velocity field of EVERY displayed frame,
        % including playback and MP4 export (not just every 60th frame).
        updateStream=showStreams;
        actualTimes=nan(nCases,1);
        for j=1:nCases
            d=S{j};
            [~,k]=min(abs(d.times-target));
            actualTimes(j)=d.times(k);
            snapshotChanged = (k~=d.currentFile);
            if snapshotChanged
                [hit,cacheLocation]=ismember(k,d.frameCacheIndices);
                if hit
                    nextT=d.frameCacheTables{cacheLocation};
                else
                    nextT=readParticleCSV(fullfile(d.files(k).folder,d.files(k).name));
                    if height(nextT)~=numel(d.xD)
                        error('Particle count changed: %s',d.files(k).name);
                    end
                end
                % Small bounded LRU cache: most recently used entry first.
                if hit
                    d.frameCacheIndices(cacheLocation)=[];
                    d.frameCacheTables(cacheLocation)=[];
                end
                d.frameCacheIndices=[k d.frameCacheIndices];
                d.frameCacheTables=[{nextT} d.frameCacheTables];
                if numel(d.frameCacheIndices)>frameCacheSize
                    d.frameCacheIndices=d.frameCacheIndices(1:frameCacheSize);
                    d.frameCacheTables=d.frameCacheTables(1:frameCacheSize);
                end
                d.table=nextT;d.currentFile=k;
                S{j}=d;
            end
            % Reuse the previously rendered field when the 0.1-s master
            % timeline points to the SAME saved CSV snapshot. This matters
            % when C++ output is saved every 0.2 s.
            if snapshotChanged || ~isPlayback
                z=fieldOfND(d.table,plotVariable,d);
                C=reshape(z(d.surfaceRows(:)),size(d.surfaceRows));
                set(SURF(j),'CData',C);
                set(TITLE(j),'String',sprintf('Re = %.0f | t = %.5g s', ...
                    d.Re,d.times(k)));
            end
            if updateStream
                [xx,yy,nPaths]=computeStreamlines(d,d.table);
                set(STREAM(j),'XData',xx,'YData',yy, ...
                    'ZData',0.01*d.radii(end)*ones(size(xx)));
                uistack(STREAM(j),'top');
                if nPaths==0
                    warning('No streamline paths for %s at t=%g',d.name,d.times(k));
                end
            end
        end
        if updateStream, lastStreamFrame=currentFrame; end
        if isgraphics(frameSlider), set(frameSlider,'Value',currentFrame); end
        if isgraphics(timeEdit),timeEdit.String=sprintf('%.3f',target);end
        set(timeLabel,'String',sprintf('Comparison t = %.5g s  |  frame %d/%d  |  actual times: %s', ...
            target,currentFrame,nFrames,strjoin(compose('%.4g',actualTimes),', ')));
        isUpdatingFrame=false;
        drawnow limitrate;
    end

    function [allX,allY,nPaths]=computeStreamlines(d,T)
        if streamlineMethod == "streamfunction" && ~isempty(d.psiSolver)
            % Uses SPH velocities of the actual selected CSV snapshot.
            % Conversion U/U_inf, V/U_inf and r/D gives psi/(U_inf*D).
            U=reshape(T.u(d.fluidRows),d.nTheta,[]);
            V=reshape(T.v(d.fluidRows),d.nTheta,[]);
            U=U(d.angleOrder,:)/d.Uinf;
            V=V(d.angleOrder,:)/d.Uinf;
            if streamfunctionVorticitySource == "solver"
                omega=reshape(T.vorticity(d.fluidRows),d.nTheta,[]);
                omega=omega(d.angleOrder,:)*d.D/d.Uinf;
            else
                omega=ewcsphPolarCurl(d.psiSolver,U,V);
            end
            psi=ewcsphSolvePolarPsi(d.psiSolver,U,V,omega);
            step=streamfunctionSpacingD;
            if ~isfinite(step) || step<=0
                error('streamfunctionSpacingD must be a positive finite number.');
            end
            zmax=max(abs(psi(:)));
            nLevel=min(800,ceil(zmax/step));
            levels=(-nLevel:nLevel)*step;
            if streamfunctionNearWallContours
                extra=abs(streamfunctionNearWallLevelsD(:)');
                extra=extra(isfinite(extra)&extra>0);
                levels=sort(unique([levels -extra extra]));
            end
            [allX,allY]=ewcsphPsiContourXY(d.psiSolver,psi,levels);
            nPaths=sum(isnan(allX));
            return;
        end
        % Original stream2() mode, unchanged.
        % Native polar-ring interpolation matches the user's viewer.
        U=reshape(T.u(d.fluidRows),d.nTheta,[]);
        V=reshape(T.v(d.fluidRows),d.nTheta,[]);
        U=U(d.angleOrder,:); V=V(d.angleOrder,:);
        angleExtended=[d.theta(end)-2*pi;d.theta;d.theta(1)+2*pi];
        U=[U(end,:);U;U(1,:)];V=[V(end,:);V;V(1,:)];
        Fu=griddedInterpolant({angleExtended,d.radii},U,'linear','none');
        Fv=griddedInterpolant({angleExtended,d.radii},V,'linear','none');
        outer=d.radii(end);inner=d.radii(1);
        [gx,gy]=meshgrid(linspace(-outer,outer,streamGridSize));
        qr=hypot(gx,gy); qa=mod(atan2(gy,gx),2*pi);
        gu=Fu(qa,qr);gv=Fv(qa,qr);
        invalid=qr<inner | qr>outer;
        gu(invalid)=NaN;gv(invalid)=NaN;
        % =========================================================
        % Paper-style streamlines: equally spaced upstream seeds,
        % with a few additional seeds in the recirculation region.
        % All locations are nondimensionalised by cylinder diameter.
        % =========================================================
        direction=[1 0];
        try
            cc=readcell(fullfile(d.folder,'SimulationConstants.csv'),'Delimiter',',');
            phi=metaNumber(cc,'flowAngleDeg')*pi/180;
            direction=[cos(phi),sin(phi)];
        catch
            % Default incident flow is along positive x.
        end
        normal=[-direction(2),direction(1)];

        % Regularly spaced seeds across the visible transverse extent.
        % Place them on a circular arc inside the radial computational domain.
        sr=0.90*outer;
        transverse=viewY_D(1):streamSeedSpacingD:viewY_D(2);
        transverse=transverse(abs(transverse)<sr-2*outer/(streamGridSize-1));
        longitudinal=-sqrt(max(0,sr^2-transverse.^2));
        sx=longitudinal*direction(1)+transverse*normal(1);
        sy=longitudinal*direction(2)+transverse*normal(2);

        % Supplement a small number of wake/recirculation streamlines.
        wakeX=wakeSeedX_D*direction(1)+wakeSeedY_D*normal(1);
        wakeY=wakeSeedX_D*direction(2)+wakeSeedY_D*normal(2);
        wr=hypot(wakeX,wakeY);
        validWake=wr>inner+2*outer/(streamGridSize-1) & ...
                  wr<outer-2*outer/(streamGridSize-1);
        wakeX=wakeX(validWake);wakeY=wakeY(validWake);
        validWake=isfinite(Fu(mod(atan2(wakeY,wakeX),2*pi), ...
                               hypot(wakeX,wakeY))) & ...
                  isfinite(Fv(mod(atan2(wakeY,wakeX),2*pi), ...
                               hypot(wakeX,wakeY)));
        wakeX=wakeX(validWake);wakeY=wakeY(validWake);

        opts2=[streamStep streamMaxVertices];
        pUpstream=stream2(gx,gy,gu,gv,sx,sy,opts2);
        pWake=stream2(gx,gy,gu,gv,wakeX,wakeY,opts2);
        % Backward tracing only for sparse wake seeds, so closed
        % recirculating patterns can be seen without crowding the field.
        pWakeBackward=stream2(gx,gy,-gu,-gv,wakeX,wakeY,opts2);
        allPaths=[pUpstream(:);pWake(:);pWakeBackward(:)];
        allX=[];allY=[];nPaths=0;
        for kk=1:numel(allPaths)
            p=allPaths{kk};
            if size(p,1)<2, continue; end
            good=isfinite(p(:,1))&isfinite(p(:,2));
            p=p(good,:);
            if size(p,1)<2, continue; end
            allX=[allX;p(:,1);NaN]; %#ok<AGROW>
            allY=[allY;p(:,2);NaN]; %#ok<AGROW>
            nPaths=nPaths+1;
        end
        allX=allX(:);allY=allY(:);
    end
end

function folders=chooseSimulationFolders()
% Desktop MATLAB: choose several sibling folders at once with Ctrl/Shift.
folders=strings(0,1);
try
    if usejava('desktop')
        chooser=javaObjectEDT('javax.swing.JFileChooser');
        chooser.setDialogTitle('Select simulation folders (Ctrl/Shift+click; Open to finish)');
        chooser.setFileSelectionMode(1); % JFileChooser.DIRECTORIES_ONLY
        chooser.setMultiSelectionEnabled(true);
        chooser.setCurrentDirectory(java.io.File(pwd));
        result=chooser.showOpenDialog([]);
        if result~=0,return;end % cancelled
        selected=chooser.getSelectedFiles();
        for k=1:numel(selected)
            folders(end+1,1)=string(char(selected(k).getAbsolutePath())); %#ok<AGROW>
        end
        if isempty(folders)
            one=chooser.getSelectedFile();
            if ~isempty(one),folders=string(char(one.getAbsolutePath()));end
        end
        return;
    end
catch ME
    warning('Multi-folder chooser unavailable (%s). Using sequential selection.',ME.message);
end
while true
    p=uigetdir(pwd,sprintf('Select comparison folder %d (Cancel = done)',numel(folders)+1));
    if isequal(p,0),break;end
    if ~any(folders==string(p)),folders(end+1,1)=string(p);end %#ok<AGROW>
end
end

function value=metaNumber(cells,key)
keys=strtrim(string(cells(2:end,1)));
idx=find(keys==string(key),1);
if isempty(idx),error('Missing SimulationConstants.csv value: %s',key);end
value=str2double(string(cells{idx+1,2}));
if ~isfinite(value),error('Invalid SimulationConstants.csv value: %s',key);end
end

function T=readParticleCSV(filename)
opts=delimitedTextImportOptions('NumVariables',16);
opts.DataLines=[7 Inf];opts.Delimiter=',';
opts.VariableNames=["ID","x","y","empty1","rho","drhodt","pressure", ...
    "empty2","u","v","velocity","empty3","dudt","dvdt","vorticity","Type"];
opts.VariableTypes=["double","double","double","string","double","double", ...
    "double","string","double","double","double","string", ...
    "double","double","double","string"];
T=readtable(filename,opts);
end

function [label,unit]=fieldMetadata(name)
switch lower(name)
    case 'vorticity',label='Vorticity, ω D/U_inf';unit='';
    case 'velocity',label='Speed, |u|/U_inf';unit='';
    case 'u',label='u/U_inf';unit='';
    case 'v',label='v/U_inf';unit='';
    case 'rho',label='rho/rho_0';unit='';
    case 'pressure',label='Pressure coefficient, C_p';unit='';
    otherwise,error('Unsupported plot variable: %s',name);
end
end
function z=fieldOfND(T,name,d)
switch lower(name)
    case 'vorticity',z=T.vorticity*d.D/d.Uinf;
    case 'velocity',z=T.velocity/d.Uinf;
    case 'u',z=T.u/d.Uinf;
    case 'v',z=T.v/d.Uinf;
    case 'rho',z=T.rho/d.rho0;
    case 'pressure',z=T.pressure/(0.5*d.rho0*d.Uinf^2);
    otherwise,error('Unsupported plot variable: %s',name);
end
end

function map=blueWhiteRed(n)
n1=floor(n/2);n2=n-n1;
map=[linspace(0,1,n1)' linspace(0,1,n1)' ones(n1,1); ...
     ones(n2,1) linspace(1,0,n2)' linspace(1,0,n2)'];
end

function S = ewcsphCreatePolarPsiSolver(theta,wallRadius,fluidRadii)
% Sparse-free spectral polar Poisson solver: FFT in theta, tridiagonal
% finite-volume operator in r. Factorisation is precomputed ONCE.
    theta=theta(:);
    r=[wallRadius;fluidRadii(:)];
    nt=numel(theta); nr=numel(r);
    if nt<4 || nr<3 || any(diff(r)<=0)
        error('Invalid polar ring geometry for streamfunction reconstruction.');
    end
    dtheta=2*pi/nt;
    if max(abs(diff([theta;theta(1)+2*pi])-dtheta))>0.05*dtheta
        error('Streamfunction reconstruction requires near-uniform angular spacing.');
    end
    i=2:nr-1;
    rm=r(i-1); ri=r(i); rp=r(i+1);
    drCell=0.5*(rp-rm);
    a=(0.5*(ri+rm))./(ri.*drCell.*(ri-rm));
    c=(0.5*(ri+rp))./(ri.*drCell.*(rp-ri));
    % Eigenvalues of the second-order periodic angular Laplacian.
    modes=(0:nt-1)';
    angularEigenvalue=4*sin(pi*modes/nt).^2/dtheta^2;
    diagCoeffs=-(a'+c')-angularEigenvalue*(1./ri'.^2);
    nw=numel(i);
    invPiv=zeros(nt,nw);
    upper=zeros(nt,nw);
    invPiv(:,1)=1./diagCoeffs(:,1);
    if nw>1
        upper(:,1)=c(1)*invPiv(:,1);
    end
    for k=2:nw
        invPiv(:,k)=1./(diagCoeffs(:,k)-a(k)*upper(:,k-1));
        if k<nw,upper(:,k)=c(k)*invPiv(:,k);end
    end
    % Signed Fourier indices for reconstructing psi at the outer boundary.
    freq=(0:nt-1)';
    freq(freq>floor(nt/2))=freq(freq>floor(nt/2))-nt;
    S=struct('theta',theta,'r',r,'nt',nt,'nr',nr, ...
             'a',a,'c',c,'invPiv',invPiv,'upper',upper, ...
             'frequency',freq,'dtheta',dtheta);
end

function omega=ewcsphPolarCurl(S,U,V)
% curl_z = (1/r)*d(r*u_theta)/dr - (1/r)*du_r/dtheta
% Reconstructs vorticity from the actual particle velocity data.
    theta=S.theta; rr=S.r(2:end)';
    ct=cos(theta);st=sin(theta);
    ur=U.*ct+V.*st;
    ut=-U.*st+V.*ct;
    A=ut.*rr;
    dA=zeros(size(A));
    dA(:,1)=(A(:,2)-A(:,1))/(rr(2)-rr(1));
    dA(:,end)=(A(:,end)-A(:,end-1))/(rr(end)-rr(end-1));
    for j=2:numel(rr)-1
        hm=rr(j)-rr(j-1); hp=rr(j+1)-rr(j);
        dA(:,j)=(-hp/(hm*(hm+hp)))*A(:,j-1) ...
             +((hp-hm)/(hm*hp))*A(:,j) ...
             +(hm/(hp*(hm+hp)))*A(:,j+1);
    end
    dUrTheta=(circshift(ur,-1,1)-circshift(ur,1,1))/(2*S.dtheta);
    omega=(dA-dUrTheta)./rr;
end

function psi=ewcsphSolvePolarPsi(S,U,V,omega)
% Delta psi=-omega; psi=0 on cylinder; outer Dirichlet psi derived from
% SPH radial velocity, not imposed from a potential-flow model.
    ct=cos(S.theta);st=sin(S.theta);
    urOuter=U(:,end).*ct+V(:,end).*st;
    urOuter=urOuter-mean(urOuter,'omitnan');
    if any(~isfinite(urOuter(:))) || any(~isfinite(omega(:)))
        error('Non-finite SPH velocities/vorticity cannot reconstruct psi.');
    end
    freq=S.frequency;
    fluxHat=fft(S.r(end)*urOuter);
    psiOuterHat=zeros(S.nt,1);
    nonzero=freq~=0;
    if rem(S.nt,2)==0
        nonzero(S.nt/2+1)=false; % Nyquist derivative is ambiguous
    end
    psiOuterHat(nonzero)=fluxHat(nonzero)./(1i*freq(nonzero));
    psiOuter=real(ifft(psiOuterHat));

    rhsHat=fft(-omega(:,1:end-1),[],1);
    rhsHat(:,end)=rhsHat(:,end)-S.c(end)*fft(psiOuter);
    nw=size(rhsHat,2);
    dp=zeros(size(rhsHat));
    dp(:,1)=rhsHat(:,1).*S.invPiv(:,1);
    for j=2:nw
        dp(:,j)=(rhsHat(:,j)-S.a(j)*dp(:,j-1)).*S.invPiv(:,j);
    end
    sol=zeros(size(rhsHat));
    sol(:,end)=dp(:,end);
    for j=nw-1:-1:1
        sol(:,j)=dp(:,j)-S.upper(:,j).*sol(:,j+1);
    end
    psiHat=[zeros(S.nt,1),sol,fft(psiOuter)];
    psi=real(ifft(psiHat,[],1));
end

function [xx,yy]=ewcsphPsiContourXY(S,psi,levels)
% Generate complete constant-psi curves on the native SPH polar rings.
% Uniform and near-zero levels share one global contour solve: no crop.
    r=S.r(:)';
    theta=[S.theta;S.theta(1)+2*pi];
    psiPeriodic=[psi;psi(1,:)];
    if numel(levels)<2
        xx=NaN;yy=NaN;return;
    end
    C=contourc(r,theta,psiPeriodic,levels);
    % MATLAB's contour matrix stores a header then vertex coordinates.
    xParts=cell(0,1); yParts=cell(0,1);
    pos=1;
    while pos<=size(C,2)
        n=round(C(2,pos));
        if ~isfinite(n)||n<2||pos+n>size(C,2),break;end
        rc=C(1,pos+1:pos+n);
        tc=C(2,pos+1:pos+n);
        xParts{end+1,1}=[rc(:).*cos(tc(:));NaN]; %#ok<AGROW>
        yParts{end+1,1}=[rc(:).*sin(tc(:));NaN]; %#ok<AGROW>
        pos=pos+n+1;
    end
    if isempty(xParts)
        xx=NaN;yy=NaN;
    else
        xx=vertcat(xParts{:});yy=vertcat(yParts{:});
    end
end
