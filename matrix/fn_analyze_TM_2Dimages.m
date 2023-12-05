function [Mod_depth_real]=fn_analyze_TM_2Dimages(BP_images,Focus_x,Focus_y,plot_cosine_fit,l,L)
% opens a measurement file of Hadamard modulation and analyzes TM
% (modulation depth amplitude and phase)
% (plot_cosine_fit=1 will plot each of the fit traces)
%
% use e.g.:
% Mod_depth_real=fn_analyze_TM_matched_filter('\\rococo\groupe\Photoacoustique\21-01-2013_data\Mod_Had1.mat',30,40,160:300,400,0);

% plot_cosine_fit=0; % don't plot each of the 160 fit curve...
valeur = zeros(size(BP_images,3),size(BP_images,4));
n_phase=size(BP_images,3);
phase_ramp=((1:n_phase)-1)/n_phase*2*pi; % just a helpful vector for fitting a cosine(phi) function
C=transpose(cos(phase_ramp));
S=transpose(sin(phase_ramp));


for k=1:size(BP_images,4) % modulation measurement in Hadamard basis (1:160)
    for n=1:size(BP_images,3)
        valeur(n,k)=max(max(BP_images(Focus_y+(-l:l),Focus_x+(-L:L),n,k)));
    end
    
    % analyze modulation depth (amp and phase) for each Hadamard k-vector
    I_to_optimize=valeur(:,k);
    cos_amp=sum(I_to_optimize.*C);
    sin_amp=sum(I_to_optimize.*S);
    Mod_depth_Hadamard(k)=cos_amp+1i*sin_amp;
    
    %% Jerome 13/11/27
    
    Err_Hadamard(k)=sum(abs(I_to_optimize'-(mean(I_to_optimize)+2*abs(Mod_depth_Hadamard(k)/n_phase)*cos(phase_ramp-angle(Mod_depth_Hadamard(k))))));
    
    % plot of each cosine fit (not a must)
    if plot_cosine_fit
        figure(23)
        plot(phase_ramp,I_to_optimize,...
            phase_ramp,mean(I_to_optimize)+2*abs(Mod_depth_Hadamard(k)/n_phase)*cos(phase_ramp-angle(Mod_depth_Hadamard(k))))
        legend('measurements','cosine fit');
        title(['k=' num2str(k)])
        drawnow
    end
end

%% convert measured modulation from Hadamard to SLM pixel basis:
Mod_depth_real = inv(hadamard(160)) * transpose(Mod_depth_Hadamard);

%% plots
figure;
subplot(2,2,1);
plot(abs(Mod_depth_Hadamard));
% hold on; plot(Err_Hadamard,'k'); hold off;

xlabel('Hadamard vector');
title('Modulation depth - Hadamard basis')
subplot(2,2,3);
dead_pixels_indexes=[1 12 133 144:160];
plot(1:160,abs(Mod_depth_real),dead_pixels_indexes,abs(Mod_depth_real(dead_pixels_indexes)),'.  r');
xlabel('real vector');
title('Modulation depth - "real" basis  (LAST 16 VALUES SHOULD =0!)')


% plot(abs(inv(hadamard(160))*(hadamard(160)*[SLM_amp_phase; zeros(16,1)])))

subplot(1,2,2);
imagesc(reshape(abs(Mod_depth_real(1:156)),12,13))
axis image
title('Modulation depth of each SLM pixel')
colorbar
colormap gray(256)
