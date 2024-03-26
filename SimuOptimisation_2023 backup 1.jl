## import packages
using Images, Plots, Hadamard, FFTW, Statistics, ProgressBars, Makie

## add missing functions
whos(x) = (size(x), typeof(x), x);

function imagesc(img)
    sc = scaleminmax(minimum(img), maximum(img))
    Gray.(sc.(img))
end

nextpow2(x) = nextprod((2,), x);

squeeze(x) = reshape(x, filter(!isone, size(x)));

had = Hadamard.hadamard;

function meshgrid(x, y)
    X = x' .* ones(size(y))
    Y = ones(size(x))' .* y
    return X, Y
end

# function cimage(complex_img) # plot a complex image with abs as hue and phase as color (red=0, blue=pi)
#     img = sc(abs.(complex_img));
#     # img = img .- minimum(img);
#     # img = img ./ maximum(img);
#     # img = img .* 255;
#     phase = sc(angle.(complex_img));
#     # phase = phase .- minimum(phase);
#     # phase = phase ./ maximum(phase);
#     # phase = phase .* 255;
#     return colorview(RGB,Gray.(img),Gray.(phase),Gray.(phase))
# end

function gen_speckle_controlsize(grid_size, speckle_size_in_macropixels)
    NA = 1 / speckle_size_in_macropixels
    x = 1:grid_size
    x = x .- size(x, 1) / 2 .- 1
    X, Y = meshgrid(x, x)
    R = sqrt.(X .^ 2 + Y .^ 2) # will be used for limiting the NA of the imaging lens
    T_medium = exp.(1im * 2 * pi * rand(size(X, 1), size(X, 2))) # scattering-medium's transmission (speckle k-space representation)
    T_medium = T_medium .* (R .<= (grid_size / 2 * NA)) # limit lens' NA (so speckles would look nice on camera, and MTF max freq. resolved)
    sp_img = abs.(fftshift(fft(fftshift(T_medium)))) .^ 2 #translate this line in julia
    return sp_img, T_medium
end

# using multiple dispatch, i.e. different methods depeneding on the number of arguments
function gen_speckle_controlsize(grid_size, speckle_size_in_macropixels, sigma)
    sp_img, T_medium = gen_speckle_controlsize(grid_size, speckle_size_in_macropixels)
    x = 1:grid_size
    x = x .- (size(x, 1) / 2) .- 1
    X, Y = meshgrid(x, x)
    g = exp.(-(X .^ 2 + Y .^ 2) ./ (2 * sigma^2))
    sp_img = sp_img .* g
    return sp_img, T_medium
end

## set parameters
# sizes in pixels
grid1Size = 400;
grainSize = 5;
SLMsize = 20;

# if rem(grid1Size,grainSize*SLMsize) !=0
if grid1Size % (grainSize * SLMsize) != 0
    error("grid1Size/grainSize/SLMsize should be integer and is currently equal to  $(grid1Size/grainSize/SLMsize)")
end

Nact = SLMsize^2; # SLM has Nact pixels
n_step = 1 * Nact;
n_phase = 3;
NL_coef = 1;

diskFlag = true;
gaussFlag = false;
squareFlag = false;

padCoord = floor(Int, 0.5 * grid1Size * (1 - 1 / grainSize)):floor(Int, 0.5 * grid1Size * (1 + 1 / grainSize));
padCoord = padCoord[2:floor(Int, grid1Size / grainSize + 1)];

## define physical parameters
lambda = 1000e-9; # in m
pixelSize = lambda / 2 / grainSize; # in m, assuming diffraction limited grains (likely even smaller in 3D)

fUS = 100e6; # in Hz
NAUS = 5;
wUS = 1500 / fUS / NAUS; # in m, ~width of transducer profile
wUS_pix = wUS / pixelSize;


## define optimisation target

if diskFlag
    radius = 20
    N = 0.5 * grid1Size
    # [X,Y] = meshgrid (-N+1:N,-N+1:N);
    x = -N+1:N
    y = -N+1:N
    X = x' .* ones(grid1Size)
    Y = ones(grid1Size)' .* y
    # to check, use Gray.(sc.(X)), after defining sc = scaleminmax(-200,200) --> write an imagesc function
    # [theta,r] = cart2pol(X,Y);
    cart2pol(x, y) = [sqrt(x^2 + y^2), atan(y, x)]
    pol = reduce(hcat, cart2pol.(X, Y))
    theta = reshape(pol[2, :], grid1Size, grid1Size)
    r = reshape(pol[1, :], grid1Size, grid1Size)
    # area = pi*radius^2;
    idx = (r .<= radius)
    imagesc(idx)
end

if gaussFlag
    N = 0.5 * grid1Size
    # [X,Y] = meshgrid(-N+1:N,-N+1:N);
    x = -N+1:N
    y = -N+1:N
    X = x' .* ones(grid1Size)
    Y = ones(grid1Size)' .* y
    radius = wUS_pix
    idx = exp.(-X .^ 2 / (2 * radius^2) .- Y .^ 2 / (2 * radius^2))
end

if squareFlag
    l = 120
    x = floor(grid1Size * 0.5)
    y = floor(grid1Size * 0.5)
    RTO_x = (x-0.5*l):(x+0.5*l)
    RTO_y = (y-0.5*l):(y+0.5*l)
    # RTO_x=x;
    # RTO_y=y;
    idx = 0 .* idx;
    idx[RTO_y, RTO_x] = 1;
    figure(1), imagesc(idx) #axis image
end

## Optimisation

n_med = 1;
optVal = fill(0.0, (n_phase, n_step, n_med));
peakVal = fill(0.0, (n_phase, n_step, n_med));

phase_ramp = ((1:n_phase) .- 1) / n_phase * 2 * pi; # just a helpful vector for fitting a cosine(phi) function
C = cos.(phase_ramp);
S = sin.(phase_ramp);
Hada = had(nextpow2(SLMsize .^ 2));
Hada = Hada .> 0; # zeros and ones
global ss = fill(0.0, (grid1Size, grid1Size,n_phase));
global s = fill(0.0, (grid1Size, grid1Size, n_med)); #zeros(grid1Size,grid1Size,1); #stocke tous les speckles successifs
global sd = fill(0.0, (grid1Size, grid1Size));
# sfin = fill(0.0, (grid1Size, grid1Size));

for imed = 1:n_med
    println("imed = $imed")

    sp, phase_mask = gen_speckle_controlsize(grid1Size, grainSize)
    s[:, :, imed] = sp
    Delta_n = fill(0., (SLMsize, SLMsize, n_step))
    pad = fill(0.0 + 0.0 * 1im, (grid1Size, grid1Size))
    max_phase = fill(0.,n_step);

    # for k = 2:n_step
    for k in ProgressBar(2:n_step)
        # println("Step: $k")
        # pixchoice=rand(SLMsize)>0.5; # random pixel selection
        pixchoice = reshape(Hada[mod(k, nextpow2(SLMsize .^ 2))+1, 1:(SLMsize^2)], SLMsize, SLMsize)';

        for n = 1:n_phase
            Delta_n[:, :, k] = mod.(Delta_n[:, :, k-1] + pixchoice * (n - 1), n_phase); # ramp the phase of chosen pixels only
            SLM_phase = kron(exp.(1im * 2 * pi * Delta_n[:, :, k] ./ n_phase), fill(1.0, (floor(Int, grid1Size / grainSize / SLMsize), floor(Int, grid1Size / grainSize / SLMsize))));

            pad[padCoord, padCoord] = SLM_phase .* phase_mask[padCoord, padCoord];

            if k == 2
                # mask[:, :, n] = Delta_n[:, :, k]
                ss[:, :, n] = abs.(fftshift(fft(fftshift(pad)))) .^ 2;
            end
            global sd = abs.(fftshift(fft(fftshift(pad)))) .^ 2;
            # if k == 2
            #     println("I was here")
    
            # end
            optVal[n, k, imed] = mean((sd .* idx) .^ NL_coef);
            peakVal[n, k, imed] = maximum((sd .* idx) .^ NL_coef);

        end

        I_to_optimize = optVal[:, k, imed];
        cos_amp = sum(I_to_optimize .* C);
        sin_amp = sum(I_to_optimize .* S);
        max_phase[k] = angle(cos_amp + 1im * sin_amp);
        maxInd = (max_phase[k] / (2 * pi) * n_phase) + 1;

        Delta_n[:, :, k] = mod.(Delta_n[:, :, k-1] + pixchoice * (maxInd - 1), n_phase); # ramp the phase of chosen pixels only
        SLM_phase = kron(exp.(1im * 2 * pi * Delta_n[:, :, k] ./ n_phase), fill(1.0, (floor(Int, grid1Size / grainSize / SLMsize), floor(Int, grid1Size / grainSize / SLMsize))));

        pad[padCoord, padCoord] = SLM_phase .* phase_mask[padCoord, padCoord];
        global sd = abs.(fftshift(fft(pad))) .^ 2
        # println(mean(sd))
    end
    println(mean(sd))
end
println(mean(sd))

## Plotting
# p1 = heatmap(s[:,:,1],aspect_ratio=:equal) # axis image,colormap jet,set(gca,'visible','off')
p1 = plot(s[:, :, 1], st=:heatmap, aspect_ratio=:equal)
# heatmap(s[:, :, 1], aspect_ratio=:equal) 

# h = drawcircle('Center',.5*[size(sd,1), size(sd,2)],'Radius',radius,...
#     'linewidth',1,'StripeColor','white','InteractionsAllowed','none');
p2 = plot(sd ./ maximum(sd), st=:heatmap, aspect_ratio=:equal, clim=(0, 1))
# heatmap(sd./maximum(sd),aspect_ratio=:equal); #axis image ,colormap jet,set(gca,'visible','off')
# subplot(133), imagesc(sdd/max(sd(:)),[0 1]), axis image,set(gca,'visible','off')
# plot(max(optVal)/mean2(s1),'linewidth',3), axis tight, ylim([0 max(optVal(:))/mean2(s1)+10])
# colormap hot
plot(p1, p2, layout=(1, 2))

MoptVal = squeeze(maximum(optVal, dims=1));
MpeakVal = squeeze(maximum(peakVal, dims=1));

Ngrains = radius^2 / (grainSize / 2)^2; # there is probably a sqrt(2) because of definition of grain size
EF_th = round(pi / 4 * Nact / Ngrains, digits=2); # see Popoff 2011 NJP for pre-factor justification

p3 = plot([MoptVal / MoptVal[2], MpeakVal / MpeakVal[2]], label=["Mean value in ROI (optimized)" "Peak value in ROI"], linewidth=3)
title!("Theoretical enhancement factor =  $EF_th", titlefontsize=10)
xlabel!("Steps")
ylabel!("Enhancement factor")
plot(p1, p2, p3, layout=(3, 1), legend=:best)

# for k = 1:n_step/Nact
#    p3 = plot!([k*Nact k*Nact],ylim,"color",[0 0 0 .2])
# end

# plot Alexis w. Makie

f2 = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (800, 600))
 
               gridlayout_mips2 = f2[1, 1] = GridLayout()
 
               ax_p1 = gridlayout_mips2[1, 1] = WGLMakie.Axis(f2, title = "Speckle 1")
               ax_p2 = gridlayout_mips2[1, 2] = WGLMakie.Axis(f2, title = "Speckle 2")              
               ax_peak = gridlayout_mips2[2, 1:2] = WGLMakie.Axis(f2, title = "Steps/Enhancement factor")
              
               ## Plotting
               heatmap!(ax_p1,s[:, :, 1], aspect_ratio=:equal)
              
               # h = drawcircle('Center',.5*[size(sd,1), size(sd,2)],'Radius',radius,...
               #     'linewidth',1,'StripeColor','white','InteractionsAllowed','none');
               heatmap!(ax_p2,sd./maximum(sd),aspect_ratio=:equal)
                             
               MoptVal = squeeze(maximum(optVal, dims=1));
               MpeakVal = squeeze(maximum(peakVal, dims=1));
              
               Ngrains = radius^2 / (grainSize / 2)^2; # there is probably a sqrt(2) because of definition of grain size
               EF_th = round(pi / 4 * Nact / Ngrains, digits=2); # see Popoff 2011 NJP for pre-factor justification
 
               max_line = [MoptVal / MoptVal[2], MpeakVal / MpeakVal[2]]
               lines!(ax_peak,MoptVal / MoptVal[2]   , label="Mean value in ROI (optimized)")
               lines!(ax_peak,MpeakVal / MpeakVal[2] , label="Peak value in ROI")
               axislegend()
               # title!("Theoretical enhancement factor =  $EF_th", titlefontsize=10)
               # xlabel!("Steps")
               # ylabel!("Enhancement factor")
               # plot(p1, p2, p3, layout=(3, 1), legend=:best)
              
               # for k = 1:n_step/Nact
               #    p3 = plot!([k*Nact k*Nact],ylim,"color",[0 0 0 .2])
               # end
              
               f2