addpath('C:\Users\Tobias\Desktop\Thomas\Boston DM-01_140') %pour utiliser les fonctions utilisées par le boston
close all
clear all
imaqreset
vid = videoinput('gige', 2, 'Mono8');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
src.TriggerActivation = 'RisingEdge';
vid.TriggerRepeat = Inf;
src.TriggerMode = 'On'; %--------->à voir si nécessaire
src.TriggerSource = 'Line1';
% création d'une ROI
% ROI=[290 190 200 200];
ROI=[0 0 780 580];

vid.ROIPosition = ROI;
%temps d'exposition
src.ExposureTimeAbs = 100250;
src.packetSize=1500;
% src.FrameRate=4;
bin_size=1;

%% Initialize the SLM
% Define two amplitude control arrays
num_actuators = 144;
n_SLM=sqrt(num_actuators);
% Open and initialize MultiDM driver USB connection
mapping_ID = 2; %IMPORTANT
[error_code, driver_info] = OPEN_multiDM(mapping_ID);

%% Définition des paramètres
n_step=160;
n_phase=3;
phase_ramp=((1:n_phase)-1)/n_phase*2*pi; % just a helpful vector for fitting a cosine(phi) function
p_ramp=((1:20)-1)/20*2*pi;
C=transpose(cos(phase_ramp));
S=transpose(sin(phase_ramp));
% u_min=42.6;
% u_max=51.5;
u_min=28;
u_max=38.2;

X=580;
Y=160;
Lmod=50;
Lvis=30;

Hada=hadamard(160);
Hada=double(Hada>0); % zeros and ones
Int = zeros(2*Lmod+1,2*Lmod+1,n_phase,n_step);
Imod=zeros(n_phase,1);
amplitudes = zeros(num_actuators,1)+u_min;
UPDATE_multiDM(driver_info, amplitudes);
% img1=double(getsnapshot(vid));
start(vid);
tic

%% Matrix

for k=1:n_step
    disp(k)
    pixchoice=Hada(k,1:(n_SLM^2));
    
    for n=1:n_phase
        flushdata(vid);
        
        SLM_phase=phase_ramp(n) + pi*pixchoice - pi;
        amplitudes=SLM_phase*(u_max-u_min)/2/pi+u_min;
        UPDATE_multiDM(driver_info, amplitudes);
        img=double(getdata(vid));
        img=double(getdata(vid));
        %Int(:,:,n,k)=binning(img,bin_size);
        Int(:,:,n,k)=img(Y+[-Lmod:Lmod],X+[-Lmod:Lmod]);
        N(n,k)=mean2(Int(:,:,n,k));
        Int(:,:,n,k)=Int(:,:,n,k)/N(n,k);
    end
    
    for i=1:(2*Lmod+1)%580%floor(ROI(4)/bin_size)
        for j=1:(2*Lmod+1)%780%floor(ROI(3)/bin_size)
            Imod(:)=Int(i,j,:,k);
            cos_amp(i,j,k)=Imod.'*C;
            sin_amp(i,j,k)=Imod.'*S;
        end
    end
    mod_hada(:,:,k)=cos_amp(:,:,k)+1i*sin_amp(:,:,k);
    %     figure(3), imagesc(abs(mod_hada(:,:,k)))

    Imod(:)=Int(28,28,:,k);
    figure(4), plot(phase_ramp,Imod), title([num2str(abs(mod_hada(i,j,k))),...
        '    ',num2str(angle(mod_hada(28,28,k)))])
    hold on
    plot(p_ramp,mean(Imod)+2*abs(mod_hada(28,28,k)/n_phase)*cos(p_ramp-angle(mod_hada(28,28,k))),'r')
    legend('measurements','cosine fit');
  xlabel(['k=' num2str(k)])
  title(num2str(max(img(:))))
    hold off
end


%% Kill all connections
% Disable and close MultiDM driver USB connection
error_code = CLOSE_multiDM(driver_info);
stop(vid);
delete(vid);
clear vid;

clear Int
c=clock;
save(['D:\Thomas_Data\Juillet 2014\',sprintf('Matrix_CAM_chickenbackplane_%i_%i_%i_%ih%i',c(1),c(2),c(3),c(4),c(5)),'.mat'],'-v7.3')


%% Analyse chicken backplane matrix
for i=1:160
figure(5), 
set(gcf,'position',[680         441        1221         537])
subplot(1,2,1),imagesc(angle(mod_hada(:,:,i))), axis image
subplot(1,2,2),imagesc(abs(mod_hada(:,:,i))), axis image

a=conv2(angle(mod_hada(:,:,i)),rot90(angle(mod_hada(:,:,i)),2),'same');
figure(6), imagesc(a)
    pause
end