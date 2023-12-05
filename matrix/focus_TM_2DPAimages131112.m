close all
% clear all
instrreset
addpath('C:\Users\Tobias\Desktop\Thomas\Boston DM-01_140') %SLM control functions
addpath('C:\Users\Tobias\Desktop\Thomas\SureliteII10') %Surelite control functions
addpath('C:\Users\Tobias\Desktop\Thomas\Juillet 2013\legHAL_r23754')
%addpath('C:\Users\Tobias\Desktop\Thomas\Juillet 2013\legHAL_r23754\scriptJerome') 

%   addPathLegHAL
%    Sequence=usse.usse;
%    Sequence= Sequence.selectProbe(); % comment if the probe did not change
%    
% filename='D:\Thomas_Data\Octobre 2013\23-10-2013_data\TM_meas_2DPA_2013_11_26_17h49.mat';
% load(filename);
clear DSO
close all;
%%

% tic
% [BP_images,SS]=Beamforming_2Dimages_130402(c,ImgInfo,NbElemts,Pitch,PAT,PA_sig);
% 
%  clear PA_sig;
% toc
% load('D:\Thomas_Data\Octobre 2013\23-10-2013_data\BPimages_feuille_2013_11_26_17h50.mat')
m_img=mean(BP_images(:,:,:),3);
figure;  imagesc(squeeze(m_img));
            colormap(gray);
            colorbar;  
title('PAT Image mean');
drawnow

% std_img=squeeze(std(BP_images(:,:,:),[],3));
     figure;  imagesc(std_img);
            colormap(gray);
            colorbar;  
title('PAT Image std');
drawnow



% save('D:\Thomas_Data\Juillet 2013\17-07-2013_data\BPimages1.mat','BP_images','SS', 'c' )
% pause
% return
%% focus region
Focus_y=127;
Focus_x=239;
Focusing.Focus_x=Focus_x;
Focusing.Focus_y=Focus_y;
Focusing.l=2;
Focusing.L=2;
hold on; plot(Focus_x,Focus_y,'or'); hold off

% figure, imagesc(BP_images(Focus_y+(-Focusing.l:Focusing.l),Focus_x+(-Focusing.L:Focusing.L),1,1));
% colormap(jet);
% colorbar;
% axis equal;
% axis tight;
% title('PAT Image');
% drawnow
%% Compute the phase masks for each wire
plot_cosine_fit=0;

Focusing.Mod_depth_real1=fn_analyze_TM_2Dimages(BP_images,Focus_x,Focus_y,plot_cosine_fit,Focusing.l,Focusing.L);
% Focusing.Mod_depth_real1=fn_analyze_TM_2Dimages_extinction(-BP_images,Focus_x,Focus_y,plot_cosine_fit,Focusing.l,Focusing.L);



%%=====================================================================
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
invoke(DSO,'WriteString',['TDIV ',num2str(tdiv)],true);
invoke(DSO,'WriteString','TRSE EDGE,SR,C2,HT,OFF',true); % Set trigger source to be another channel, for example channel 2

%% Laser
ComPort=7;
% disp('did you check the laser power???')
% pause
QS_delay=240;
Controller=sureliteII10OpenAndInitialise(ComPort,QS_delay); %set the right power


%% Définition des paramètres
n_avg=20;

amplitudes = zeros(num_actuators,1)+u_min;
invoke(DSO,'WriteString',['SEQ ON,',num2str(n_avg)],true); % Turn On Sequence Mode 'seq on ou off,nb de segments, longueur max d'un segment en unité de sample'



%% Aixplorer
ini_aixplorer_forSLM131023; % Aixplorer
startupcl;
[b ,a] =butter(3,[5e6 20e6]/(fs/2));

Sequence = Sequence.initializeRemote('IPaddress', AixplorerIP, 'InitTGC',PAT.TGC, 'InitRxFreq',PAT.RxFreq ); % Aixplorer
% Load sequence
Sequence = Sequence.loadSequence(); % Aixplorer


%% Focusing
% flat SLM
amplitudes = u_min + zeros(1,144);
UPDATE_multiDM(driver_info, amplitudes);
invoke(DSO,'WriteString','ARM;WAIT',true); %the scope waits for the next trigger event
clear buffer; % Aixplorer
Sequence = Sequence.startSequence('Wait', 0); % Aixplorer
sureliteII10SendCommand(Controller, 'SH 1'); %open the shutter
buffer = Sequence.getData('Realign',1); % Aixplorer
sureliteII10SendCommand(Controller, 'SH 0'); %open the shutter
Sequence = Sequence.stopSequence('Wait', 1);% focus on 1st target:
PD_sig=invoke(DSO,'GetScaledWaveformWithTimes','C3',1000000,1); % Get Wavefrom Data - 1 Mpts is chosen as the arbitrary maximum file transfer file size (substitute your own max value to use)
for u=1:n_avg
    N_norm(u)=sum(PD_sig(2,((u-1)*(Nacq+2)+1+4990):(u*(Nacq+2)-4910))); %on doit restreindre la zone d'int sinon trop de bruit
end

i_norm=find(N_norm>(max(N_norm)/2));


sigMat = single(buffer(1).data(:,(1:128)+round((NbElemts-128)/2),:))/4;      
        for i=1:size(sigMat,3)
            sigMat(:,:,i)=filtfilt(b,a,double(squeeze(sigMat(:,:,i))));
            sigMat(:,:,i) = single(sigMat(:,:,i)./N_norm(i));
         end

        LsigMat=size(sigMat,1);
        
        for i=1:size(sigMat,3)
           MatCorr=zeros(2*size(sigMat,1)-1,size(sigMat,2));
            for j=1:size(sigMat,2)
                MatCorr(:,j)=xcorr( sigMat_ref(:,j),squeeze(sigMat(:,j,i)),'coeff');
            end
            
            [H hi]=max(MatCorr(size(sigMat,1)+(-5:5),:),[],1);
            iid=find(H>(max(H)/2));
            decal=round(mean(hi(iid)-6));
          if decal>=0
                
                SigMat1 =  cat(1, zeros(decal,size(sigMat,2),'single'), squeeze(sigMat(:,:,i)));
                sigMat(:,:,i)=SigMat1(1:LsigMat,:);
          else
                SigMat1 =  cat(1,  squeeze(sigMat(:,:,i)), zeros(-decal,size(sigMat,2),'single'));
                sigMat(:,:,i)=SigMat1((-decal)+(1:LsigMat),:);
                            end
            clear SigMat1 decal
        end

        sigMat=mean(sigMat(:,:,i_norm),3);      
        
        sigMat = cat(1, zeros(ind_plus,size(sigMat,2),'single'), sigMat);


Focusing.sigMat_flat = sigMat;

sigMat=diff(sigMat);
       
        Sb = reshape( squeeze( sigMat), size( sigMat, 1 ), size( sigMat, 2 ) ) ;
                
        bp_image = backprojects3d_cl( Sb, n_image, P, c, 1, t(1:(end-1))', fs, image_width );
     
Focusing.bp_image_flat=bp_image;
% figure(11);imagesc(bp_image);
% colormap(jet);
% colorbar;
% axis equal;
% axis tight;
% title('Flat SLM');
% hold on; plot(Focus_x,Focus_y,'or'); hold off;
% drawnow
clear sigMat;

saveas(gcf,[sprintf('PAimg_flat131125',Focus_x,Focus_y),'.fig'])

%% focusing phase mask
cur_phases1=angle(Focusing.Mod_depth_real1);
amplitudes1 = u_min + (u_max-u_min)*mod(cur_phases1(1:144),2*pi)/2/pi; 

UPDATE_multiDM(driver_info, amplitudes1);
invoke(DSO,'WriteString','ARM;WAIT',true); %the scope waits for the next trigger event
clear buffer; % Aixplorer
Sequence = Sequence.startSequence('Wait', 0); % Aixplorer
sureliteII10SendCommand(Controller, 'SH 1'); %open the shutter
buffer = Sequence.getData('Realign',1); % Aixplorer
sureliteII10SendCommand(Controller, 'SH 0'); %open the shutter
Sequence = Sequence.stopSequence('Wait', 1);% focus on 1st target:
PD_sig=invoke(DSO,'GetScaledWaveformWithTimes','C3',1000000,1); % Get Wavefrom Data - 1 Mpts is chosen as the arbitrary maximum file transfer file size (substitute your own max value to use)
for u=1:n_avg
    N_norm(u)=sum(PD_sig(2,((u-1)*(Nacq+2)+1+4990):(u*(Nacq+2)-4910))); %on doit restreindre la zone d'int sinon trop de bruit
end
i_norm=find(N_norm>(max(N_norm)/2));

clear sigMat
sigMat = single(buffer(1).data(:,(1:128)+round((NbElemts-128)/2),:))/4;      
        for i=1:size(sigMat,3)
            sigMat(:,:,i)=filtfilt(b,a,double(squeeze(sigMat(:,:,i))));
            sigMat(:,:,i) = single(sigMat(:,:,i)./N_norm(i));
         end

        LsigMat=size(sigMat,1);
        
        for i=1:size(sigMat,3)
           MatCorr=zeros(2*size(sigMat,1)-1,size(sigMat,2));
            for j=1:size(sigMat,2)
                MatCorr(:,j)=xcorr( sigMat_ref(:,j),squeeze(sigMat(:,j,i)),'coeff');
            end
            
            [H hi]=max(MatCorr(size(sigMat,1)+(-5:5),:),[],1);
            iid=find(H>(max(H)/2));
            decal=round(mean(hi(iid)-6));
          if decal>=0
                
                SigMat1 =  cat(1, zeros(decal,size(sigMat,2),'single'), squeeze(sigMat(:,:,i)));
                sigMat(:,:,i)=SigMat1(1:LsigMat,:);
          else
                SigMat1 =  cat(1,  squeeze(sigMat(:,:,i)), zeros(-decal,size(sigMat,2),'single'));
                sigMat(:,:,i)=SigMat1((-decal)+(1:LsigMat),:);
                            end
            clear SigMat1 decal
        end

        sigMat=mean(sigMat(:,:,i_norm),3);      
        
        sigMat = cat(1, zeros(ind_plus,size(sigMat,2),'single'), sigMat);


Focusing.sigMat_focus = sigMat;

sigMat=diff(sigMat);
       
        Sb = reshape( squeeze( sigMat), size( sigMat, 1 ), size( sigMat, 2 ) ) ;
                
        bp_image = backprojects3d_cl( Sb, n_image, P, c, 1, t(1:(end-1))', fs, image_width );

Focusing.bp_image_focus=bp_image;

figure(22);
set(gcf,'position',[256          49        1383         948])
subplot(2,2,1)

imagesc(Focusing.bp_image_flat);
colormap(jet);
colorbar;
axis equal;
axis tight;
title('Flat SLM');
hold on; plot(Focus_x,Focus_y,'or'); hold off;
drawnow

subplot(2,2,2)
imagesc(bp_image);
colormap(jet);
colorbar;
axis equal;
axis tight;
title('Focusing SLM');
hold on; plot(Focus_x,Focus_y,'ok'); hold off
drawnow

subplot(2,2,3)
imagesc(abs(Focusing.bp_image_focus-Focusing.bp_image_flat));
colormap(jet);
colorbar;
axis equal;
axis tight;
title('difference with flat slm');
hold on; plot(Focus_x,Focus_y,'ok'); hold off
drawnow
subplot(2,2,4)
imagesc(abs(Focusing.bp_image_focus-m_img));
colormap(jet);
colorbar;
axis equal;
axis tight;
title('difference with mean image');
hold on; plot(Focus_x,Focus_y,'ok'); hold off
drawnow

saveas(gcf,[sprintf('PAimg_focus131125_x%i_y%i',Focus_x,Focus_y),'.fig'])
save([sprintf('PAimg131125_feuille_x%i_y%i',Focus_x,Focus_y),'.mat'],'Focusing')
clear sigMat;

error_code = CLOSE_multiDM(driver_info);
%scope
invoke(DSO,'Disconnect');
close(666);
%laser
sureliteII10StopLaser(Controller);


return
%% SCan focus with memory effect

mask=reshape(cur_phases1(1:144),12,12);
mask2=mask;
[X, Y]=meshgrid(1:length(mask(:,1)));
step_width=0.03;

% while not(s==27) %%esc to exit the scan mode
%     figure(22)
%     disp('which direction?')
% s=getkey;% gauche: 28% haut: 30% droite: 29% bas: 31
% if s==30
% mask2=mod(mask+step_width*Y,2*pi);
% end
% if s==31
% mask2=mod(mask+step_width*flipud(Y),2*pi);
% end
% if s==29
% mask2=mod(mask+step_width*X,2*pi);
% end
% if s==28
% mask2=mod(mask+step_width*fliplr(X),2*pi);
% end
mask2=mod(mask+step_width*fliplr(X-1),2*pi);%left
% mask2=mod(mask+step_width*X,2*pi);%right
% mask2=mod(mask+step_width*Y,2*pi);%up?
% mask2=mod(mask+step_width*flipud(Y),2*pi);%down?


amplitudes = u_min + (u_max-u_min)*reshape(mask2,144,1)/2/pi;

UPDATE_multiDM(driver_info, amplitudes);
invoke(DSO,'WriteString','ARM;WAIT',true); %the scope waits for the next trigger event
clear buffer; % Aixplorer
Sequence = Sequence.startSequence('Wait', 0); % Aixplorer
sureliteII10SendCommand(Controller, 'SH 1'); %open the shutter
buffer = Sequence.getData('Realign',1); % Aixplorer
sureliteII10SendCommand(Controller, 'SH 0'); %open the shutter
Sequence = Sequence.stopSequence('Wait', 1);% focus on 1st target:
PD_sig=invoke(DSO,'GetScaledWaveformWithTimes','C1',1000000,1); % Get Wavefrom Data - 1 Mpts is chosen as the arbitrary maximum file transfer file size (substitute your own max value to use)
for u=1:n_avg
    N_norm(u)=sum(PD_sig(2,((u-1)*(Nacq+2)+1+4990):(u*(Nacq+2)-4910))); %on doit restreindre la zone d'int sinon trop de bruit
end
sigMat = buffer(1).data(:,(1:128)+round((NbElemts-128)/2),:)/4;
clear sigMat2;
ii=1;
for i=1:size(sigMat,3)
    if N_norm(i)> (max(N_norm)*0.7)
        sigMat2(:,:,ii) = double(sigMat(:,:,i))./N_norm(i); %sigMAt en int16 par défaut, double pour ne pas arrondir
        ii=ii+1;
    end
end
sigMat=mean(sigMat2,3);
Focusing.sigMat_focus=sigMat;

bp_image = zeros(size(SS,1), size(SS,2));
sigMat=-diff(sigMat);
 [b a]=butter(3,20/30,'low');
sigMat= filtfilt(b,a,sigMat);

for j = 1:length(x_array)
    A0 = sigMat(:,j);
    s1 = A0(squeeze(SS(:,:,j)));
    bp_image = bp_image + s1;
end ;
Focusing.bp_image_focus_tilt=bp_image;
disp('focus is moving')

figure(32)
subplot(2,1,1)
imagesc(-Focusing.bp_image_focus_tilt);
colormap(jet);
colorbar;
axis equal;
axis tight;
title('Focusing SLM');
hold on; plot(Focus_x,Focus_y,'+k'); hold off
drawnow

subplot(2,1,2)
imagesc(abs(Focusing.bp_image_focus_tilt-Focusing.bp_image_flat));
colormap(jet);
colorbar;
axis equal;
axis tight;
title('difference');
hold on; plot(Focus_x,Focus_y,'+k'); hold off
drawnow
saveas(gcf,[sprintf('PAimg_focusmemeff11left_step%i_x%i_y%i',step_width,Focus_x,Focus_y),'.fig'])
save([sprintf('PAimg170713_focusmemeff11left_step%i_x%i_y%i',step_width,Focus_x,Focus_y),'.mat'],'Focusing')

figure(33)
imagesc(abs(Focusing.bp_image_focus_tilt-Focusing.bp_image_focus));
colormap(jet);
colorbar;
axis equal;
axis tight;
% end
%% Kill all connections
%SLM
error_code = CLOSE_multiDM(driver_info);
%scope
invoke(DSO,'Disconnect');
close(666);
%laser
sureliteII10StopLaser(Controller);


c=clock;
save(['D:\Thomas_Data\Octobre 2013\23-10-2013_data\',sprintf('focus_2DPA_%i_%i_%i_%ih%i',c(1),c(2),c(3),c(4),c(5)),'.mat'],'-v7.3')



