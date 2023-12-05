close all
clear all
instrreset
addpath('C:\Users\Tobias\Desktop\Thomas\Boston DM-01_140') %SLM control functions
addpath('C:\Users\Tobias\Desktop\Thomas\SureliteII10') %Surelite control functions

%% Initialization
%SLM
num_actuators = 144;
mapping_ID = 2; %IMPORTANT
[error_code, driver_info] = OPEN_multiDM(mapping_ID);
%Scope
SamplingFreq= 1e9;     % Fréquence d'échantillonnage
Nacq= 10e3;  % Nombres d'échantillons à acquérir
tdiv = Nacq / (SamplingFreq * 10);
figure(666)
DSO = actxcontrol('LeCroy.ActiveDSOCtrl.1'); % Load ActiveDSO control
adr='IP:10.10.20.27';
invoke(DSO,'MakeConnection',adr) % Returns True on success, False on failure.Substitute your choice of IP address here
invoke(DSO,'WriteString','BUZZ BEEP',true);  %Beep if connection is successful
invoke(DSO,'WriteString','TDIV num2str(tdiv)');
%Laser
ComPort=4;
QS_delay=420;
Controller=sureliteII10OpenAndInitialise(ComPort,QS_delay); %set the right power

%Aixplorer
%????

%% Définition des paramètres
n_phase=8;
n_avg=10;
N=zeros(n_avg,1);
invoke(DSO,'WriteString','SEQ ON,n_avg,250',true); % Turn On Sequence Mode 'seq on ou off,nb de segments, longueur max d'un segment en unité de sample'
%le 250 est encore à éclaircir, rien d'autre ne marche
phase_ramp=((1:n_phase)-1)/n_phase*2*pi; % just a helpful vector for fitting a cosine(phi) function
C=transpose(cos(phase_ramp));
S=transpose(sin(phase_ramp));
u_min=42.2;
u_max=52.2;
Hada=double(hadamard(160)>0);
valeur = zeros(n_phase,num_actuators);
Delta_n = zeros(num_actuators,num_actuators);
amplitudes = zeros(num_actuators,1)+u_min;
%focus region
Focus_x=284:365;
Focus_y=55:674;
PD_sig_ind=1980:2040;

tic
for k=1:160
    disp(k);
    pixchoice=Hada(k,1:num_actuators)';
    for n=1:n_phase
        cur_phases=phase_ramp(n) + pi*pixchoice - pi;
        amplitudes=u_min + (u_max-u_min)*mod(cur_phases,2*pi)/2/pi;
        UPDATE_multiDM(driver_info, amplitudes);
        %shot laser, trig scope, record n_avg PA signals
        invoke(DSO,'WriteString','ARM;WAIT',true); %the scope waits for the next trigger event
        sureliteII10SendCommand(Controller, 'SH 1'); %open the shutter
        shotsCounter = sureliteII10Query(Controller, 'SC'); %attention à cette histoire de command echoed, à essayer
        shotsCounter2 = sureliteII10Query(Controller, 'SC');
        while (shotsCounter2-shotsCounter)<n_avg
            shotsCounter2 = sureliteII10Query(Controller, 'SC');
        end
        %scope en mode single+seq->donne les n_avg premiers shots (donc
        %c'est parfait si la sonde PA fait la même chose)
        PD_sig=invoke(DSO,'GetScaledWaveformWithTimes','C1',1000000,1); % Get Wavefrom Data - 1 Mpts is chosen as the arbitrary maximum file transfer file size (substitute your own max value to use)
        for u=1:10
            N(u)=sum(PD_sig(2,((u-1)*(Nacq+2)+1+4960):(u*(Nacq+2)-4940))); %on doit restreindre la zone d'int sinon trop de bruit
%             figure(12)
%             plot(PD_sig(2,((u-1)*(Nacq+2)+1+4960):(u*(Nacq+2)-4940)));
        end
        %Normaliser la uieme image par N(u)
        %ici??
        acquire_PA_images(n_avg); %on préferait avoir n_avg images pour normaliser chacune avant d'en prendre la moyenne
        PA_image_avg=blabla;
        valeur(n,k)=sum(sum(PA_image_avg(Focus_y,Focus_x)));
        %       valeur(n,k)=mean2(PA_image_avg(Focus_y,Focus_x));
        
        if (k==1)&&(n==1)
            figure(56)
            plot(avg_y_PA)
            hold on
            plot(avg_y_PD,'r')
            legend('PA','PD')
            hold off
            disp('press any key')
            figure(1000)
            plot(y_PA(:,n,k));
            hold on
            plot(sig_ind_for_valeur-PA_sig_ind(1)+1,avg_y_PA(sig_ind_for_valeur),'r')
            hold off
            pause
        end
    end
    %% fit the current iteration results to a cosine and find its phase:
    % (not sure a cosine is good/right for broadband input)------> In the
    % calib of the SLM, make sure that we scan 2*pi interval for the
    % "highest" frequency
    fit=1; %pour si on utilisait ça ou pas plus tard
    I_to_optimize=valeur(:,k);
    cos_amp=sum(I_to_optimize.*C);
    sin_amp=sum(I_to_optimize.*S);
    max_phase=angle(cos_amp+1i*sin_amp);
    Mod_depth(k)=cos_amp+i*sin_amp;
    figure(23)
    plot(phase_ramp,I_to_optimize,...
        phase_ramp,mean(I_to_optimize)+2*abs(Mod_depth(k)/n_phase)*cos(phase_ramp-max_phase))
    legend('measurements','cosine fit');
    title(['k=' num2str(k)])
end
temps_passe = toc
%% Kill all connections
%SLM
error_code = CLOSE_multiDM(driver_info);
%scope
invoke(DSO,'Disconnect');
close(666);
%laser
sureliteII10StopLaser(Controller);

%%
Mod_depth_real=inv(hadamard(160))*transpose(Mod_depth);

figure;
subplot(2,1,1); plot(abs(Mod_depth)); xlabel('Hadamard vector'); title('Modulation depth - Hadamard basis')
subplot(2,1,2); plot(1:160,abs(Mod_depth_real)); xlabel('real vector'); title('Modulation depth - "real" basis')
% plot(abs(inv(hadamard(160))*(hadamard(160)*[SLM_amp_phase; zeros(16,1)])))

figure;
% subplot(2,1,1);
imagesc(reshape(abs(Mod_depth_real(1:144)),12,12))
title('Modulation depth of each SLM pixel')
colorbar
colormap gray(256)

save('D:\Thomas_Data\Mars 2013\11-03-2013_data\TM_meas_10avg_16pts1.mat')

