% close all
clear all
clearvars -except mysys

%% set parameters

% sizes in pixels
grid1Size = 400;
grainSize = 5;
SLMsize = 20;

if rem(grid1Size,grainSize*SLMsize) ~=0
    error(['grid1Size/grainSize/SLMsize should be integer and is currently equal to ' num2str(grid1Size/grainSize/SLMsize)])
end

Nact = SLMsize^2; % SLM has Nact pixels
n_step=2*Nact;
upd_bar = textprogressbar21(n_step, 'barlength', 20, 'updatestep', 1, 'startmsg', ['Optimising - Step:'],'endmsg', ' Done.', ...
    'showbar', false, 'showremtime', true, 'showactualnum', true, 'barsymbol', '=','emptybarsymbol', '-');

n_phase=3;
s=zeros(grid1Size,grid1Size,1); %stocke tous les speckles successifs
NL_coef=1;

diskFlag = 0;
gaussFlag = 1;
squareFlag = 0;

padCoord = floor(0.5*grid1Size*(1-1/grainSize)):floor(0.5*grid1Size*(1+1/grainSize));
padCoord = padCoord(2:grid1Size/grainSize+1);

%% define physical parameters
lambda = 1000e-9; % in m
pixelSize = lambda/2/grainSize; % in m, assuming diffraction limited grains (likely even smaller in 3D)

fUS = 100e6; % in Hz
NAUS = 5;
wUS = 1500/fUS/NAUS; % in m, ~width of transducer profile
wUS_pix = wUS/pixelSize;


%% define optimisation target

if diskFlag
    radius = 20;
    N = 0.5*grid1Size;
    [X,Y] = meshgrid(-N+1:N,-N+1:N);
    [theta,r] = cart2pol(X,Y);
    % area = pi*radius^2;
    idx = double(r<=radius);
    figure(1), imagesc(idx), axis image
end

if gaussFlag
    N = 0.5*grid1Size;
    [X,Y] = meshgrid(-N+1:N,-N+1:N);
    radius=wUS_pix;
    idx=exp(-X.^2/(2*radius^2)-Y.^2/(2*radius^2));
end

if squareFlag
    l=120;
    x=floor(grid1Size*0.5);
    y=floor(grid1Size*0.5);
    RTO_x=(x-0.5*l):(x+0.5*l);
    RTO_y=(y-0.5*l):(y+0.5*l);
    % RTO_x=x;
    % RTO_y=y;
    idx = zeros*idx;
    idx(RTO_y,RTO_x) = 1;
    figure(1), imagesc(idx), axis image
end

%% Optimisation

n_med=1;
optVal = zeros(n_phase,n_step,n_med);
phase_ramp=((1:n_phase)-1)/n_phase*2*pi; % just a helpful vector for fitting a cosine(phi) function
C=transpose(cos(phase_ramp));
S=transpose(sin(phase_ramp));
Hada=hadamard(2^nextpow2(SLMsize.^2));
Hada=double(Hada>0); % zeros and ones
ss=zeros(grid1Size,grid1Size);

for imed=1:n_med
    disp(['imed = ' num2str(imed)])
    
    [s(:,:,1,imed),phase_mask]=gen_speckle_controlsize(grid1Size,grainSize);
    Delta_n = zeros(SLMsize,SLMsize,n_step);
    pad=zeros(grid1Size);
    
    for k=2:n_step
        % pixchoice=rand(SLMsize)>0.5; % random pixel selection
        pixchoice=reshape(Hada(mod(k,2^nextpow2(SLMsize.^2))+1,1:(SLMsize^2)),SLMsize,SLMsize).;
        
        for n=1:n_phase
            Delta_n(:,:,k) = mod(Delta_n(:,:,k-1) + pixchoice*(n-1),n_phase); % ramp the phase of chosen pixels only
            SLM_phase=kron(exp(1i*2*pi*Delta_n(:,:,k)/n_phase),ones(grid1Size/grainSize/SLMsize));
            
            pad(padCoord,padCoord)=SLM_phase.*phase_mask(padCoord,padCoord);
            
            if k==2
                mask(:,:,n)= Delta_n(:,:,k);
                ss(:,:,n)=abs(fftshift(fft2(fftshift(pad)))).^2;
            end
            sd=abs(fftshift(fft2(fftshift(pad)))).^2;
            optVal(n,k,imed)=mean2((sd.*idx).^NL_coef);
            peakVal(n,k,imed)=max(max((sd.*idx).^NL_coef));
            
        end
        
        I_to_optimize=optVal(:,k,imed);
        cos_amp=sum(I_to_optimize.*C);
        sin_amp=sum(I_to_optimize.*S);
        max_phase(k)=angle(cos_amp+1i*sin_amp);
        maxInd=(max_phase(k)/(2*pi)*n_phase)+1;
        
        Delta_n(:,:,k) = mod(Delta_n(:,:,k-1) + pixchoice*(maxInd-1),n_phase); % ramp the phase of chosen pixels only
        SLM_phase=kron(exp(1i*2*pi*Delta_n(:,:,k)/n_phase),ones(grid1Size/grainSize/SLMsize));
        
        pad(padCoord,padCoord)=SLM_phase.*phase_mask(padCoord,padCoord);
        sd=abs(fftshift(fft2(pad))).^2;
        
        upd_bar(k)
        
    end
end
% return
%%
% circlep=(g>=0.45*max(g(:)));
% circlem=(g<=0.55*max(g(:)));
% c=(circlep+circlem)==2;

figure(1),
subplot(221), imagesc(s(:,:,1)), axis image,colormap jet,set(gca,'visible','off')

h = drawcircle('Center',.5*[size(sd,1), size(sd,2)],'Radius',radius,...
    'linewidth',1,'StripeColor','white','InteractionsAllowed','none');

subplot(222), imagesc(sd/max(sd(:)),[0 1]), axis image ,colormap jet,set(gca,'visible','off')
% subplot(133), imagesc(sdd/max(sd(:)),[0 1]), axis image,set(gca,'visible','off')
% plot(max(optVal)/mean2(s1),'linewidth',3), axis tight, ylim([0 max(optVal(:))/mean2(s1)+10])
colormap hot

xlabel('Steps')
ylabel('Enhancement factor')

MoptVal=max(optVal);
MpeakVal=max(peakVal);

subplot(2,2,[3,4]), plot(MoptVal/MoptVal(2))
hold on
plot(MpeakVal/MpeakVal(2))
for k = 1:n_step/Nact
    plot([k*Nact k*Nact],ylim,'color',[0 0 0 .2])
end
hold off
legend({'Mean value in ROI (optimized)', 'Peak value in ROI'},'location','northwest')

Ngrains = radius^2/(grainSize/2)^2; % there is probably a sqrt(2) because of definition of grain size
EF_th = pi/4*Nact/Ngrains; % see Popoff 2011 NJP for pre-factor justification
title(['Theoretical enhancement factor = ' num2str(EF_th)])


%% Auxiliary functions

function upd = textprogressbar21(n, varargin)
% UPD = TEXTPROGRESSBAR(N) initializes a text progress bar for monitoring a
% task comprising N steps (e.g., the N rounds of an iteration) in the
% command line. It returns a function handle UPD that is used to update and
% render the progress bar. UPD takes a single argument i <= N which
% corresponds to the number of tasks completed and renders the progress bar
% accordingly.
%
% TEXTPROGRESSBAR(...,'barlength',L) determines the length L of the
% progress bar in number of characters (see 'barsymbol' option). L must be
% a positive integer.
% (Default value is 20 characters.)
%
% TEXTPROGRESSBAR(...,'updatestep',S) determines the minimum number of update
% steps S between consecutive bar re-renderings. The option controls how
% frequently the bar is rendered and in turn controls the computational
% overhead that the bar rendering incurs to the monitored task. It is
% especially useful when bar is used for loops with large number of rounds
% and short execution time per round.
% (Default value is S=10 steps.)
%
% TEXTPROGRESSBAR(...,'startmsg',str) determines the message string to be
% displayed before the progress bar.
% (Default is str='Completed '.)
%
% TEXTPROGRESSBAR(...,'endmsg',str) determines the message string to be
% displayed after progress bar when the task is completed.
% (Default is str=' Done.')
%
% TEXTPROGRESSBAR(...,'showremtime',b) logical parameter that controls
% whether an estimate of the remaining time is displayed.
% (Default is b=true.)
%
% TEXTPROGRESSBAR(...,'showbar',b) logical parameter that controls whether
% the progress bar is displayed. (Default is b=true.)
%
% TEXTPROGRESSBAR(...,'showpercentage',b) logical parameter that controls
% whether to display the percentage of completed items.
% (Default is true.)
%
% TEXTPROGRESSBAR(...,'showactualnum',b) logical parameter that controls
% whether to display the actual number of completed items.
% (Default is false.)
%
% TEXTPROGRESSBAR(...,'showfinaltime',b) logical parameter that controls
% whether to display the total run-time when completed.
% (Default is true.)
%
% TEXTPROGRESSBAR(...,'barsymbol',c) determines the symbol (character) to
% be used for the progress bar. c must be a single character.
% (Default is c='='.)
%
% TEXTPROGRESSBAR(...,'emptybarsymbol',c) determines the symbol (character)
% that is used to fill the un-completed part of the progress bar. c must be
% a single character.
% (Default is c=' '.)
%
% Example:
%
%   n = 150;
%   upd = textprogressbar(n);
%   for i = 1:n
%      pause(0.05);
%      upd(i);
%   end
%

% Default Parameter values:
defaultbarCharLen = 20;
defaultUpdStep = 10;
defaultstartMsg = 'Completed ';
defaultendMsg = ' Done.';
defaultShowremTime = true;
defaultShowBar = true;
defaultshowPercentage = true;
defaultshowActualNum = false;
defaultshowFinalTime = true;
defaultbarCharSymbol = '=';
defaultEmptybarCharSymbol = ' ';

% Auxiliary functions for checking parameter values:
ischarsymbol = @(c) (ischar(c) && length(c) == 1);
ispositiveint = @(x) (isnumeric(x) && mod(x, 1) == 0 && x > 0);

% Register input parameters:
p = inputParser;
addRequired(p,'n', ispositiveint);
addParameter(p, 'barlength', defaultbarCharLen, ispositiveint)
addParameter(p, 'updatestep', defaultUpdStep, ispositiveint)
addParameter(p, 'startmsg', defaultstartMsg, @ischar)
addParameter(p, 'endmsg', defaultendMsg, @ischar)
addParameter(p, 'showremtime', defaultShowremTime, @islogical)
addParameter(p, 'showbar', defaultShowBar, @islogical)
addParameter(p, 'showpercentage', defaultshowPercentage, @islogical)
addParameter(p, 'showactualnum', defaultshowActualNum, @islogical)
addParameter(p, 'showfinaltime', defaultshowFinalTime, @islogical)
addParameter(p, 'barsymbol', defaultbarCharSymbol, ischarsymbol)
addParameter(p, 'emptybarsymbol', defaultEmptybarCharSymbol, ischarsymbol)

% Parse input arguments:
parse(p, n, varargin{:});
n = p.Results.n;
barCharLen = p.Results.barlength;
updStep = p.Results.updatestep;
startMsg = p.Results.startmsg;
endMsg = p.Results.endmsg;
showremTime = p.Results.showremtime;
showBar = p.Results.showbar;
showPercentage = p.Results.showpercentage;
showActualNum = p.Results.showactualnum;
showFinalTime = p.Results.showfinaltime;
barCharSymbol = p.Results.barsymbol;
emptybarCharSymbol = p.Results.emptybarsymbol;

% Initialize progress bar:
bar = ['[', repmat(emptybarCharSymbol, 1, barCharLen), ']'];

nextRenderPoint = 0;
startTime = tic;

% Initalize block for actual number of completed items:

ind = 1;

% Start message block:
startMsgLen = length(startMsg);
startMsgStart = ind;
startMsgEnd = startMsgStart + startMsgLen - 1;
ind = ind + startMsgLen;

% Bar block:
barLen = length(bar);
barStart = 0;
barEnd = 0;
if showBar
    barStart = ind;
    barEnd = barStart + barLen - 1;
    ind = ind + barLen;
end

% Actual Num block:
actualNumDigitLen = numel(num2str(n));
actualNumFormat = sprintf(' %%%dd/%d', actualNumDigitLen, n);
actualNumStr = sprintf(actualNumFormat, 0);
actualNumLen = length(actualNumStr);
actualNumStart = 0;
actualNumEnd = 0;
if showActualNum
    actualNumStart = ind;
    actualNumEnd = actualNumStart + actualNumLen-1;
    ind = ind + actualNumLen;
end

% Percentage block:
percentageFormat = sprintf(' %%3d%%%%');
percentageStr = sprintf(percentageFormat, 0);
percentageLen = length(percentageStr);
percentageStart = 0;
percentageEnd = 0;
if showPercentage
    percentageStart = ind;
    percentageEnd = percentageStart + percentageLen-1;
    ind = ind + percentageLen;
end

% Remaining Time block:
remTimeStr = time2str(Inf);
remTimeLen = length(remTimeStr);
remTimeStart = 0;
remTimeEnd = 0;
if showremTime
    remTimeStart = ind;
    remTimeEnd = remTimeStart + remTimeLen - 1;
    ind = ind + remTimeLen;
end


% End msg block:
endMsgLen = length(endMsg);
if showBar
    endMsgStart = barEnd + 1; % Place end message right after bar;
else
    endMsgStart = startMsgEnd + 1;
end
endMsgEnd = endMsgStart + endMsgLen - 1;

ind = max([ind, endMsgEnd]);

% Determine size of buffer:
arrayLen = ind - 1;
array = repmat(' ', 1, arrayLen);

% Initial render:
array(startMsgStart:startMsgEnd) = sprintf('%s', startMsg);

delAll = repmat('\b', 1, arrayLen);

% Function to update the status of the progress bar:
    function update(i)
        
        if i < nextRenderPoint
            return;
        end
        if i > 0
            fprintf(delAll);
        end
        %pause(1)
        nextRenderPoint = min([nextRenderPoint + updStep, n]);
        
        if showremTime
            % Delete remaining time block:
            array(remTimeStart:remTimeEnd) = ' ';
        end
        
        if showPercentage
            % Delete percentage block:
            array(percentageStart:percentageEnd) = ' ';
        end
        
        if showActualNum
            % Delete actual num block:
            array(actualNumStart:actualNumEnd) = ' ';
        end
        
        if showBar
            % Update progress bar (only if needed):
            barsToPrint = floor( i / n * barCharLen );
            bar(2:1+barsToPrint) = barCharSymbol;
            array(barStart:barEnd) = bar;
        end
        
        % Check if done:
        if i >= n
            array(endMsgStart:endMsgEnd) = endMsg;
            array(endMsgEnd+1:end) = ' ';
            
            if showFinalTime
                finalTimeStr = ...
                    sprintf(' [%d seconds]', round(toc(startTime)));
                finalTimeLen = length(finalTimeStr);
                if endMsgEnd + finalTimeLen < arrayLen
                    array(endMsgEnd+1:endMsgEnd+finalTimeLen) = ...
                        finalTimeStr;
                else
                    array = [array(1:endMsgEnd), finalTimeStr];
                end
            end
            
            fprintf('%s', array);
            fprintf('\n');
            return;
        end
        
        if showActualNum
            % Delete actual num block:
            actualNumStr = sprintf(actualNumFormat, i);
            array(actualNumStart:actualNumEnd) = actualNumStr;
        end
        
        if showPercentage
            % Render percentage block:
            percentage = floor(i / n * 100);
            percentageStr = sprintf(percentageFormat, percentage);
            array(percentageStart:percentageEnd) = percentageStr;
        end
        
        % Print remaining time block:
        if showremTime
            t = toc(startTime);
            remTime = t/ i * (n-i);
            remTimeStr = time2str(remTime);
            array(remTimeStart:remTimeEnd) = remTimeStr;
        end
        fprintf('%s', array);
    end

% Do the first render:
update(0);

upd = @update;

end

function timestr = time2str(t)

if t == Inf
    timestr = sprintf(' --:--:--');
else
    [hh, mm, tt] = sec2hhmmss(t);
    timestr = sprintf(' %02d:%02d:%02d', hh, mm, tt);
end
end

function [hh, mm, ss] = sec2hhmmss(t)
hh = floor(t / 3600);
t = t - hh * 3600;
mm = floor(t / 60);
ss = round(t - mm * 60);
end

function [nn,err]=plot_norm(t,varargin)
nn=nargin;
err=varargin;
% if nn==2
%     err{1}='color';
%     a=rand;
%     pause(0.01)
%     b=rand;
%     pause(0.01)
%     c=rand;
%     err{2}=[a b c];
%
% % return
% plot(x,y/max(y(:)),err{1},err{2});%,err{3},err{4})
%
% else
%     plot(x,y/max(y(:)),err{1})%,err{2});%,err{3},err{4})

if nargin==1
    y=t;
    x=1:length(y);
else
    y=varargin{1};
    x=t;
end

plot(x,y/max(y(:)))%,err{2});%,err{3},err{4})

end

function [sp_img,T_medium]=gen_speckle_controlsize(grid_size,speckle_size_in_macropixels, varargin)
NA=1/speckle_size_in_macropixels;

x=1:(grid_size);
x=x-length(x)/2-1;
[X,Y]=meshgrid(x);
R=sqrt(X.^2+Y.^2); % will be used for limiting the NA of the imaging lens

T_medium=exp(1i*2*pi*rand(size(X))); % scattering-medium's transmission (speckle k-space representation)
T_medium=T_medium.*(R<=(grid_size/2*NA)); % limit lens' NA (so speckles would look nice on camera, and MTF max freq. resolved)
sp_img=abs(fftshift(fft2(fftshift(T_medium)))).^2;

if nargin>2
    sigma=varargin{1};
    g=exp(-(X.^2+Y.^2)/2/sigma^2);
    sp_img=sp_img.*g;
end
end
