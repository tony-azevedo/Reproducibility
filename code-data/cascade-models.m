%----------------------------------------------------------------------------------
% TRANSDUCTION MODELS
%----------------------------------------------------------------------------------

%%
% define plot color sequence, axis fonts
PlotColors = 'bgrkymcbgrkymcbgrkymcbgrkymc';
set(0, 'DefaultAxesFontName','Helvetica')
set(0, 'DefaultAxesFontSize', 12)
colormap([0 0 0])
scrsz = get(0, 'ScreenSize');

%%
%----------------------------------------------------------------------------------
% THEORETICAL PREDICTIONS
%----------------------------------------------------------------------------------

% CV of response areas vs number of steps in rhodopsin inactivation and
% number of transducin activations

TAct = 1:1000;
RhSteps = [2 4 8 16 32 64];
clear CVArea;

for step = 1:length(RhSteps)
    CVArea(:, step) = sqrt(1/RhSteps(step) + 2./TAct);
end
figure(1);

semilogx(CVArea);
xlabel('transducin activations');
ylabel('CV_{area}');
legend('m=2', '4', '8', '16', '32', '64');

%%
cd ~/Dropbox/Reproducibility/code-data 

cd /Users/Shared/Dropbox/Reproducibility/code-data/

%%
%----------------------------------------------------------------------------------
% LIKELIHOOD PARAMETER FITTING
%----------------------------------------------------------------------------------
%%
cd ~/Dropbox/Reproducibility/code-data

%%
% scan across several parameters and assess likelihood of each given
% properties of experimental responses

NumSteps =  [8 10 12];              % number shutoff steps
RhTCon = [200 250 300 350 400];             % rhodopsin decay time constant in ms
RhDecayFact = [0.8 0.9 1];  % extend of rhodopsin decay for each step
NumResponses = 2000;       % with 5000 resps simulation std ~1%, about 10% std in likelihood

load GCAPSingles36
Resp = Target.AveSingle;
FreqCutoff = 10;                                  
NumLoops = 1;
Verbose = 1;
CovarFit = 0;                   % if 0, likelihood assumes parameters independent

for loop = 1:NumLoops
    
    clear TPeakRatio VarWidth Likelihood CV;

    for step = 1:length(NumSteps)

        for tcon = 1:length(RhTCon)

            for decay = 1:length(RhDecayFact)

                RespLength = length(Resp);
                tme = 1:RespLength;                          % time in ms

                % run stochastic simulation to get mean rhodopsin time course
                clear SimulateCondition;
                SimulateCondition.RhSteps = NumSteps(step);
                SimulateCondition.RhDecayFact = RhDecayFact(decay);
                SimulateCondition.SamplingInterval = 1;         % ms
                SimulateCondition.EpochPts = length(Resp);
                SimulateCondition.EpochNumbers(1:NumResponses) = 0;
                SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
                SimulateCondition.RhodopsinCompression = 0;
                SimulateCondition.ResponseCompression = 0;
                SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps(step) / RhTCon(tcon);
                SimulateCondition.NumPtsToPeak = 100;
                SimulateCondition.GainFact = 1;
                SimulateCondition.Verbose = 0;
                
                SimulateCondition = RhCascadeModelTA(SimulateCondition, Resp, NumResponses);
                
                % compute transduction cascade filter
                filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanStochTimeCourse);
                filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
                WTFilt = real(ifft(filtFFT));    

                % generate responses
                SimulateCondition = RhCascadeModelTA(SimulateCondition, WTFilt, NumResponses);
                
                % compare modeled and measured mean squared and variance
                CV(decay, tcon, step) = std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData'));

                SimulateCondition.AverageResponse = mean(SimulateCondition.EpochData);
                SimulateCondition.VarianceResponse = var(SimulateCondition.EpochData);

                [MaxVal, MeanTPeak] = max(SimulateCondition.AverageResponse);
                [MaxVal, VarTPeak] = max(SimulateCondition.VarianceResponse);
                TPeakRatio(decay, tcon, step) = VarTPeak / MeanTPeak;
                Indices = find(SimulateCondition.VarianceResponse > 0.5*MaxVal);
                FirstTime = Indices(min(find(Indices < VarTPeak)));
                SecondTime = Indices(max(find(Indices > VarTPeak)));
                VarWidth(decay, tcon, step) = (SecondTime - FirstTime) / MeanTPeak;
                
                % calculate likelihood of model parameters
                if (CovarFit)       % take into account covariance of measured parameters
                    ParamsVec(1) = CV(decay, tcon, step) - Target.MeanCVArea; 
                    ParamsVec(2) = TPeakRatio(decay, tcon, step) - Target.MeanTPeakRatio;
                    ParamsVec(3) = VarWidth(decay, tcon, step) - Target.MeanVarWidth;
                    ParamsProjs = ParamsVec * Target.ParamsEigVec;
                    Likelihood(decay, tcon, step) = normpdf(abs(ParamsProjs(1)/Target.ParamsEigSEM(1)));
                    Likelihood(decay, tcon, step) = Likelihood(decay, tcon, step) * normpdf(abs(ParamsProjs(2)/Target.ParamsEigSEM(2)));;
                    Likelihood(decay, tcon, step) = Likelihood(decay, tcon, step) * normpdf(abs(ParamsProjs(3)/Target.ParamsEigSEM(3)));
                else    % assume parameters independent
                    Likelihood(decay, tcon, step) = normpdf((abs(CV(decay, tcon, step) - Target.MeanCVArea) / Target.SEMCVArea));
                    Likelihood(decay, tcon, step) = Likelihood(decay, tcon, step) * normpdf((abs(TPeakRatio(decay, tcon, step) - Target.MeanTPeakRatio) / Target.SEMTPeakRatio));
                    Likelihood(decay, tcon, step) = Likelihood(decay, tcon, step) * normpdf((abs(VarWidth(decay, tcon, step) - Target.MeanVarWidth) / Target.SEMVarWidth));
                end
                
                if (Verbose)
                    fprintf('%d %d %d\n', step, tcon, decay); 
                    fprintf('\tCV = %d TPeakRatio = %d VarWidth = %d Likelihood = %d\n', CV(decay, tcon, step), TPeakRatio(decay, tcon, step), VarWidth(decay, tcon, step), Likelihood(decay, tcon, step));
                end
                
            end

        end

    end
    [MaxVal, MaxLoc] = max(Likelihood(:));
    MaxLoc = MaxLoc-1;
    Index1 = floor(MaxLoc / (length(RhDecayFact) * length(RhTCon)));
    Index2 = floor((MaxLoc - Index1 * (length(RhDecayFact) * length(RhTCon))) / length(RhDecayFact));
    Index3 = MaxLoc - Index1 * (length(RhDecayFact) * length(RhTCon)) - Index2 * length(RhDecayFact);

    fprintf(1, '%d %d %d\n', Index1+1, Index2+1, Index3+1);
    fprintf(1, '%d %d %d\n', NumSteps(Index1+1), RhTCon(Index2+1), RhDecayFact(Index3+1));

end

ContourLevels = [-4:0.5:0] + log10(max(Likelihood(:)));
figure(2);clf

for step = 1:length(NumSteps)
    subplot(1, length(NumSteps), step);
    contour(RhTCon, RhDecayFact, log10(Likelihood(:, :, step)), sort(ContourLevels))
    temp = strcat(num2str(NumSteps(step)), ' steps');
    text(RhTCon(2), RhDecayFact(2), temp);
end
subplot(1, length(NumSteps), 1)
ylabel('Rh decay factor');
subplot(1, length(NumSteps), ceil(length(NumSteps)/2))
xlabel('Rh* duration');


%%
% scan across several parameters and assess likelihood of each given
% properties of experimental responses - transducin compression

NumSteps = [12];
RhTCon = [100 200 300];
TranTCon = 100;
RhDecayFact = 0.9;
NumResponses = 5000;       % with 5000 resps simulation std ~1%, about 10% std in likelihood
load WTSingles36
Resp = Target.AveSingle;
FreqCutoff = 10;                                
NumLoops = 1;
Verbose = 1;
TransRate = [0.001 0.01 0.1 1];
TransCompression = [0];
TransRateCompressFlag = 0;
CovarFit = 0;

for loop = 1:NumLoops
    
    clear TPeakRatio VarWidth Likelihood CV TransCount;

    for step = 1:length(NumSteps)

        for tcon = 1:length(RhTCon)

            for comp = 1:length(TransRate)

                RespLength = length(Resp);
                tme = 1:RespLength;                          % time in ms

                % run stochastic simulation to get mean rhodopsin time
                % course
                clear SimulateCondition;
                SimulateCondition.RhSteps = NumSteps(step);
                SimulateCondition.RhDecayFact = RhDecayFact;
                SimulateCondition.SamplingInterval = 1;         % ms
                SimulateCondition.EpochPts = length(Resp);
                SimulateCondition.EpochNumbers(1:NumResponses) = 0;
                SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
                SimulateCondition.RhodopsinCompression = 0;
                SimulateCondition.ResponseCompression = 0;
                SimulateCondition.DeterministicRhModel = 0;
                SimulateCondition.DeterministicTrModel = 0;
                SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps(step) / RhTCon(tcon);
                SimulateCondition.NumPtsToPeak = 100;
                SimulateCondition.TransCompression = 0;
                SimulateCondition.TransRateCompression = TransCompression;
                SimulateCondition.BaseTransRate = TransRate(comp);                  % cannot exceed 1 - decrease sampling interval to get higher rate
                SimulateCondition.GainFact = 1;
                SimulateCondition.Verbose = 0;
                SimulateCondition.TransDecayTimeConst = TranTCon / SimulateCondition.SamplingInterval;
                SimulateCondition.NumPtsToPeak = 100;
                SimulateCondition.ConvFlag = 0;
                
                % multi-step shutoff
                tic
                SimulateCondition = RhPlusTransducinCascadeModel(SimulateCondition, Resp, NumResponses);
                toc
                
                % compute transduction cascade filter
                filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanTrTimeCourse);
                filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
                WTFilt = real(ifft(filtFFT));    
                 
                % multi-step shutoff
                tic
                SimulateCondition = RhPlusTransducinCascadeModel(SimulateCondition, WTFilt, NumResponses);
                toc
                
                % compare modeled and measured mean squared and variance
                CV(comp, tcon, step) = std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData'));

                SimulateCondition.AverageResponse = mean(SimulateCondition.EpochData);
                SimulateCondition.VarianceResponse = var(SimulateCondition.EpochData);

                [MaxVal, MeanTPeak] = max(SimulateCondition.AverageResponse);
                [MaxVal, VarTPeak] = max(SimulateCondition.VarianceResponse);
                TPeakRatio(comp, tcon, step) = VarTPeak / MeanTPeak;
                Indices = find(SimulateCondition.VarianceResponse > 0.5*MaxVal);
                FirstTime = Indices(min(find(Indices < VarTPeak)));
                SecondTime = Indices(max(find(Indices > VarTPeak)));
                VarWidth(comp, tcon, step) = (SecondTime - FirstTime) / MeanTPeak;
                TransCount(comp, tcon, step) = mean(SimulateCondition.TransCount);
                if (CovarFit)       % take into account covariance of measured parameters
                    ParamsVec(1) = CV(comp, tcon, step) - Target.MeanCVArea; 
                    ParamsVec(2) = TPeakRatio(comp, tcon, step) - Target.MeanTPeakRatio;
                    ParamsVec(3) = VarWidth(comp, tcon, step) - Target.MeanVarWidth;
                    ParamsProjs = ParamsVec * Target.ParamsEigVec;
                    Likelihood(comp, tcon, step) = normpdf(abs(ParamsProjs(1)/Target.ParamsEigSEM(1)));
                    Likelihood(comp, tcon, step) = Likelihood(comp, tcon, step) * normpdf(abs(ParamsProjs(2)/Target.ParamsEigSEM(2)));;
                    Likelihood(comp, tcon, step) = Likelihood(comp, tcon, step) * normpdf(abs(ParamsProjs(3)/Target.ParamsEigSEM(3)));
                else    % assume parameters independent
                    Likelihood(comp, tcon, step) = normpdf((abs(CV(comp, tcon, step) - Target.MeanCVArea) / Target.SEMCVArea));
                    Likelihood(comp, tcon, step) = Likelihood(comp, tcon, step) * normpdf((abs(TPeakRatio(comp, tcon, step) - Target.MeanTPeakRatio) / Target.SEMTPeakRatio));
                    Likelihood(comp, tcon, step) = Likelihood(comp, tcon, step) * normpdf((abs(VarWidth(comp, tcon, step) - Target.MeanVarWidth) / Target.SEMVarWidth));
                end

                if (Verbose)
                    fprintf('%d %d %d %d\n', step, tcon, comp, mean(SimulateCondition.TransCount)); 
                    fprintf('\tCV = %d TPeakRatio = %d VarWidth = %d Likelihood = %d\n', CV(comp, tcon, step), TPeakRatio(comp, tcon, step), VarWidth(comp, tcon, step), Likelihood(comp, tcon, step));
                end

            end

        end
    end

    [MaxVal, MaxLoc] = max(Likelihood(:));
    MaxLoc = MaxLoc-1;
    Index1 = floor(MaxLoc / (length(TransRate) * length(RhTCon)));
    Index2 = floor((MaxLoc - Index1 * (length(TransRate) * length(RhTCon))) / length(TransRate));
    Index3 = MaxLoc - Index1 * (length(TransRate) * length(RhTCon)) - Index2 * length(TransRate);
    
    fprintf(1, '%d %d %d\n', Index1+1, Index2+1, Index3+1);
    fprintf(1, '%d %d %d\n', NumSteps(Index1+1), RhTCon(Index2+1), TransCompression(Index3+1));
end


ContourLevels = [-4:0.5:0] + log10(max(Likelihood(:)));
figure(1);clf

for step = 1:length(NumSteps)
    subplot(1, length(NumSteps), step);
    contour(RhTCon, TransRate, log10(Likelihood(:, :, step)), sort(ContourLevels))
    temp = strcat(num2str(NumSteps(step)), ' steps');
    text(RhTCon(2), TransRate(2), temp);
end
subplot(1, length(NumSteps), 1)
ylabel('Tr compress factor');
subplot(1, length(NumSteps), ceil(length(NumSteps)/2))
xlabel('Rh* duration');

%%
% variability in rhodopsin vs transducin

NumSteps = [6 8 10 12 14];
RhTCon = 250;
TranTCon = 100;
RhDecayFact = 1;
NumResponses = 1000;       % with 5000 resps simulation std ~1%, about 10% std in likelihood
Resp = WTAveNormSingle36C;
Target = WTSingles36;
FreqCutoff = 10;                                
NumLoops = 1;
Verbose = 1;
TransRate = [0.1 0.2 0.4 0.6 0.8];
TransCompression = 0;

for loop = 1:NumLoops
    
    clear TPeakRatio VarWidth Likelihood CV;

    for step = 1:length(NumSteps)

        for rate = 1:length(TransRate)

            RespLength = length(Resp);
            tme = 1:RespLength;                          % time in ms

            % run stochastic simulation to get mean rhodopsin time course
            clear SimulateCondition;
            SimulateCondition.NumSteps = NumSteps(step);
            SimulateCondition.RhDecayFact = RhDecayFact;
            SimulateCondition.SamplingInterval = 1;         % ms
            SimulateCondition.EpochPts = length(Resp);
            SimulateCondition.EpochNumbers(1:NumResponses) = 0;
            SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
            SimulateCondition.RhodopsinCompression = 0;
            SimulateCondition.ResponseCompression = 0;
            SimulateCondition.DeterministicRhModel = 0;
            SimulateCondition.DeterministicTrModel = 0;
            SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps(step) / RhTCon;
            SimulateCondition.NumPtsToPeak = 100;
            SimulateCondition.BaseTransRate = TransRate(rate);                  % cannot exceed 1 - decrease sampling interval to get higher rate
            SimulateCondition.TransCompression = TransCompression;
            SimulateCondition.TransDecayTimeConst = TranTCon / SimulateCondition.SamplingInterval;
            SimulateCondition.NumPtsToPeak = 100;

            % multi-step shutoff
            SimulateCondition = RhPlusTransducinCascadeModel(SimulateCondition, Resp, NumResponses);

            % compute transduction cascade filter
            filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanTrTimeCourse);
            filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
            WTFilt = real(ifft(filtFFT));    

            % multi-step shutoff
            SimulateCondition = RhPlusTransducinCascadeModel(SimulateCondition, WTFilt, NumResponses);

            % compare modeled and measured mean squared and variance
            CV(rate, step) = std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData'));

            SimulateCondition.AverageResponse = mean(SimulateCondition.EpochData);
            SimulateCondition.VarianceResponse = var(SimulateCondition.EpochData);

            [MaxVal, MeanTPeak] = max(SimulateCondition.AverageResponse);
            [MaxVal, VarTPeak] = max(SimulateCondition.VarianceResponse);
            TPeakRatio(rate, step) = VarTPeak / MeanTPeak;
            Indices = find(SimulateCondition.VarianceResponse > 0.5*MaxVal);
            FirstTime = Indices(min(find(Indices < VarTPeak)));
            SecondTime = Indices(max(find(Indices > VarTPeak)));
            VarWidth(rate, step) = SecondTime - FirstTime;
            Likelihood(rate, step) = normpdf((abs(CV(rate, step) - Target.MeanCVArea) / Target.SEMCVArea));
            Likelihood(rate, step) = Likelihood(rate, step) * normpdf((abs(TPeakRatio(rate, step) - Target.MeanTPeakRatio) / Target.SEMTPeakRatio));
            Likelihood(rate, step) = Likelihood(rate, step) * normpdf((abs(VarWidth(rate, step) - Target.MeanVarWidth) / Target.SEMVarWidth));

            if (Verbose)
                fprintf('%d %d\n', rate, comp); 
                fprintf('\tCV = %d TPeakRatio = %d VarWidth = %d Likelihood = %d\n', CV(rate, step), TPeakRatio(rate, step), VarWidth(rate, step), Likelihood(rate, step));
            end

        end
    end

%     [MaxVal, MaxLoc] = max(Likelihood(:));
%     MaxLoc = MaxLoc-1;
%     Index1 = floor(MaxLoc / (length(TransCompression) * length(TranTCon)));
%     Index2 = floor((MaxLoc - Index1 * (length(TransCompression) * length(TranTCon))) / length(TransCompression));
%     Index3 = MaxLoc - Index1 * (length(TransCompression) * length(TranTCon)) - Index2 * length(TransCompression);
%     
%     fprintf(1, '%d %d %d\n', Index1+1, Index2+1, Index3+1);
%     fprintf(1, '%d %d %d\n', NumSteps(Index1+1), TranTCon(Index2+1), TransCompression(Index3+1));
end

%%
%----------------------------------------------------------------------------------
% LIKELIHOOD ANALYSIS USING PCA
%----------------------------------------------------------------------------------

%%
% scan across several parameters and assess likelihood of each given
% properties of experimental responses

NumSteps =  [8 10 12];
RhTCon = [400 450 500 550];
RhDecayFact = [0.6 0.8 1];
NumResponses = 2000;       % with 5000 resps simulation std ~1%, about 10% std in likelihood
load GCAPSingles36
Resp = Target.AveSingle;
NormResp = Target.AveNormSingle;
FreqCutoff = 10;                                
NumLoops = 1;
Verbose = 1;
NumEigVec = 3;

clear EigVec Likelihood;

for vec = 1:NumEigVec
    EigVec(1:length(NormResp), vec) = Target.EigVec(1:length(NormResp), vec);
end
CompVarMean = mean(Target.VarSinglesProjections - Target.VarFailuresProjections);
CompVarSEM = std(Target.VarSinglesProjections - Target.VarFailuresProjections)/sqrt(size(Target.VarSinglesProjections, 1));

for loop = 1:NumLoops
            
    RespLength = length(Resp);
    tme = 1:RespLength;                          % time in ms

    for step = 1:length(NumSteps)
        for tcon = 1:length(RhTCon)
            for decay = 1:length(RhDecayFact)

                % run stochastic simulation to get mean rhodopsin time course
                clear SimulateCondition;
                SimulateCondition.NumSteps = NumSteps(step);
                SimulateCondition.RhDecayFact = RhDecayFact(decay);
                SimulateCondition.SamplingInterval = 1;         % ms
                SimulateCondition.EpochPts = length(Resp);
                SimulateCondition.EpochNumbers(1:NumResponses) = 0;
                SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
                SimulateCondition.RhodopsinCompression = 0;
                SimulateCondition.ResponseCompression = 0;
                SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps(step) / RhTCon(tcon);
                SimulateCondition.NumPtsToPeak = 100;

                SimulateCondition = RhCascadeModel(SimulateCondition, Resp, NumResponses);

                % compute transduction cascade filter
                filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanRhTimeCourse);
                filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
                WTFilt = real(ifft(filtFFT));    

                % generate responses
                SimulateCondition = RhCascadeModel(SimulateCondition, WTFilt, NumResponses);

                % determine projections along each PC
                Singles = SimulateCondition.NormSingles(1:length(NormResp), :);
                SimSinglesProjs = Singles' * EigVec;
                VarSimSinglesProjs = var(SimSinglesProjs);
                
                % compare modeled and measured mean squared and variance
                CV(decay, tcon, step) = std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData'));

                Likelihood(decay, tcon, step) = normpdf((abs(CV(decay, tcon, step) - Target.MeanCVArea) / Target.SEMCVArea));

                for vec = 1:NumEigVec
                         Likelihood(decay, tcon, step) = Likelihood(decay, tcon, step) * normpdf((VarSimSinglesProjs(vec) - CompVarMean(vec)) / CompVarSEM(vec));
                end

                if (Verbose)
                    fprintf('%d %d %d CV = %d Likelihood = %d\n', step, tcon, decay, CV(decay, tcon, step), Likelihood(decay, tcon, step));
                end
    
            end
        end
    end
    [MaxVal, MaxLoc] = max(Likelihood(:));
    MaxLoc = MaxLoc-1;
    Index1 = floor(MaxLoc / (length(RhDecayFact) * length(RhTCon)));
    Index2 = floor((MaxLoc - Index1 * (length(RhDecayFact) * length(RhTCon))) / length(RhDecayFact));
    Index3 = MaxLoc - Index1 * (length(RhDecayFact) * length(RhTCon)) - Index2 * length(RhDecayFact);

    fprintf(1, '%d %d %d\n', Index1+1, Index2+1, Index3+1);
    fprintf(1, '%d %d %d\n', NumSteps(Index1+1), RhTCon(Index2+1), RhDecayFact(Index3+1));
end

ContourLevels = [-4:0.5:0] + log10(max(Likelihood(:)));
figure(2);clf

for step = 1:length(NumSteps)
    subplot(1, length(NumSteps), step);
    contour(RhTCon, RhDecayFact, log10(Likelihood(:, :, step)), sort(ContourLevels))
    temp = strcat(num2str(NumSteps(step)), ' steps');
    text(RhTCon(2), RhDecayFact(2), temp);
end
subplot(1, length(NumSteps), 1)
ylabel('Rh decay factor');
subplot(1, length(NumSteps), ceil(length(NumSteps)/2))
xlabel('Rh* duration');



%%
%**************************************************************************
% CONSTRAINTS ON LINEARITY OF TRANSDUCTION PROCESS
%**************************************************************************
%---------------------------------------------------------------------
% Check linearity of transduction process by looking at consistency of 
% average single photon responses across genetic manipulations of R*
% inactivation.  First compute linear filter connecting R*(t) and I(t) for
% one genetic background.  Then use this filter to predict response to step
% in R* activity (i.e. in absence of inactivation).  
%
% Created 7/10 FMR
%
% Dependencies: requires average single photon responses from wild-type,
% GRK1+/- and CSM rods.  Sampling time step is 1 ms.
%---------------------------------------------------------------------

% assuming linear transduction cascade
RespLength = length(AveSingles.RKAveSingle);

Resp = AveSingles.RKAveSingle(1:RespLength);   % base response to compute R*->I filter
tme = 1:RespLength;                          % time in ms
tcon1 = 640;                                      % R*(t) decay time constant for "base" response - units ms
tcon2 = 300;                                      % R*(t) decay time constant for first predicted response - units ms
FreqCutoff = 20;                                

% compute R*->I filter for "base" response
Rh = exp(-tme/tcon1);

filtFFT = fft(Resp) ./ fft(Rh);
filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
filt = real(ifft(filtFFT));
pred = real(ifft(fft(Rh) .* fft(filt)));        % sanity check - essentially checks frequency cutoff

% predict CSM (no inactivation) response as integral of filter (i.e. 
% R*(t) = 1 for all t)
pred2 = cumsum(filt);

% predict response with other R* decay time constant
RespLength = length(AveSingles.WTAveSingle);
Resp = AveSingles.RKAveSingle(1:RespLength);   % base response to compute R*->I filter
tme2 = 1:RespLength;                          % time in ms

% compute R*->I filter for "base" response
Rh = exp(-tme2/tcon1);
filtFFT = fft(Resp) ./ fft(Rh);
filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
filt = real(ifft(filtFFT));

Rh = exp(-tme2/tcon2);
pred1 = real(ifft(fft(Rh) .* fft(filt)));

% plot
figure(1);
plot(tme2, filt);
figure(2); clf
pred1 = pred1 + 0.1;
plot(tme, AveSingles.RKAveSingle(1:length(tme)), 'k', tme2, AveSingles.WTAveSingle(1:length(tme2)), 'b', tme, AveSingles.CSMAveSingle(1:length(tme)), 'r', tme, pred, 'k', tme2, pred1, 'b', tme, pred2, 'r');
set(2, 'DefaultAxesFontSize', 11)
legend('GRK+/-', 'wt', 'CSM', 'GRK+/- pred', 'wt pred', 'CSM pred', 'Location', 'East');
legend boxoff
set(2, 'DefaultAxesFontSize', 14)
xlabel('ms');
ylabel('pA');
xlim([0 2000]);


%%
%  Incorporate simple thresholding nonlinearity into transduction cascade
%  Nonlinearity is applied to PDE activity - i.e. all values of P* above
%  threshold are set to threshold.  

RespLength = length(RKAveSingle);
Resp = RKAveSingle(1:RespLength);   % base response to compute R*->I filter
tme = 1:RespLength;                          % time in ms

FreqCutoff = 10;
Thresh = 1.0;                             % threshold for nonlinearity - applied to PDE activity.  Units are relative
                                                 % to max PDE activity for "base" response
tcon1 = 640;
tcon = 300;                                 % T* decay time constant
tcon2 = 300;

Rh = exp(-tme/tcon1);
temp = exp(-tme/tcon);              % T*->PDE filter

P = real(ifft(fft(Rh) .* fft(temp)));   % PDE activity from convolution of R*(t) and T*->P filter
maxP = max(P);                             
Indices = find(P > maxP * Thresh);  % apply threshold
P(Indices) = maxP * Thresh;

filtFFT = fft(Resp) ./ fft(P);
filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
filt = real(ifft(filtFFT));
pred = real(ifft(fft(P) .* fft(filt)));

P = cumsum(temp);
P(length(tme)) = 0;
Indices = find(P > maxP * Thresh);
P(Indices) = maxP * Thresh;
pred2 = conv(filt, P, 'full');

RespLength = length(WTAveNormSingle);
Resp = RKAveSingle(1:RespLength);   % base response to compute R*->I filter
tme2 = 1:RespLength;                          % time in ms

Rh = exp(-tme2/tcon1);
temp = exp(-tme2/tcon);              % T*->PDE filter

P = real(ifft(fft(Rh) .* fft(temp)));   % PDE activity from convolution of R*(t) and T*->P filter
Indices = find(P > maxP * Thresh);  % apply threshold
P(Indices) = maxP * Thresh;

filtFFT = fft(Resp) ./ fft(P);
filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
filt = real(ifft(filtFFT));

Rh = exp(-tme2/tcon2);

P = real(ifft(fft(Rh) .* fft(temp)));
Indices = find(P > maxP * Thresh);
P(Indices) = maxP * Thresh;
pred1 = real(ifft(fft(P) .* fft(filt)));

figure(1); clf
plot(tme2, filt);
figure(2); clf
pred1 = pred1 + 0.1;
plot(tme, RKAveSingle(1:length(tme)), 'k', tme2, WTAveNormSingle(1:length(tme2)), 'b', tme, CSMAveSingle(1:length(tme)), 'r', tme, pred, 'k', tme2, pred1, 'b', tme, pred2(1:length(tme)), 'r');
set(2, 'DefaultAxesFontSize', 11)
legend('GRK+/-', 'wt', 'CSM', 'GRK+/- pred', 'wt pred', 'CSM pred', 'Location', 'East');
legend boxoff
set(2, 'DefaultAxesFontSize', 14)
xlabel('ms');
ylabel('pA');
xlim([0 2000]);

%%
% nonlinearity with transducin compression

Resp = AveSingles.RKAveSingle;   % base response to compute R*->I filter
temp = [1:6*length(Resp)]*0;
Resp = [Resp temp];
RespLength = length(Resp);
tme = 1:RespLength;                          % time in ms
tcon1 = 600;                                      % R*(t) decay time constant for "base" response - units ms
FreqCutoff = 4000;                                
NumSteps = 10;
RhDecayFact = 1;
NumResponses = 100;
TransCompression = [0 1e-4 1e-3 1e-2 1e-1 1];
TransRateCompressFlag = 0;
TranTCon = 50;
TransRate = 1;

clear Likelihood;

for tcomp = 1:length(TransCompression)
    % run stochastic simulation to get mean rhodopsin time course
    clear SimulateCondition;
    SimulateCondition.RhSteps = NumSteps;
    SimulateCondition.RhDecayFact = RhDecayFact;
    SimulateCondition.SamplingInterval = 1;         % ms
    SimulateCondition.EpochPts = length(Resp);
    SimulateCondition.EpochNumbers(1:NumResponses) = 0;
    SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
    SimulateCondition.RhodopsinCompression = 0;
    SimulateCondition.ResponseCompression = 0;
    SimulateCondition.NumPtsToPeak = 100;
    SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / tcon1;
    SimulateCondition.GainFact = 1;
    SimulateCondition.Verbose = 0;
    SimulateCondition.ConvFlag = 1;
    
    if (TransRateCompressFlag)
        SimulateCondition.TransCompression = 0;
        SimulateCondition.TransRateCompression = TransCompression(tcomp);
    else
        SimulateCondition.TransCompression = TransCompression(tcomp);
        SimulateCondition.TransRateCompression = 0;
    end
    SimulateCondition.DeterministicRhModel = 1;
    SimulateCondition.DeterministicTrModel = 1;
    SimulateCondition.BaseTransRate = TransRate;                  % cannot exceed 1 - decrease sampling interval to get higher rate
    SimulateCondition.TransDecayTimeConst = TranTCon / SimulateCondition.SamplingInterval;

    if (0)
        SimulateCondition = RhCascadeModelTA(SimulateCondition, Resp, NumResponses);
        filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanStochTimeCourse);
        filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
        WTFilt = real(ifft(filtFFT));
        pred = real(ifft(fft(WTFilt) .* fft(SimulateCondition.MeanStochTimeCourse)));
        pred2 = cumsum(WTFilt);
    else          
        SimulateCondition = RhPlusTransducinCascadeModel(SimulateCondition, Resp, NumResponses);
        filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanTrTimeCourse);
        filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
        WTFilt = real(ifft(filtFFT));
        pred = real(ifft(fft(WTFilt) .* fft(SimulateCondition.MeanTrTimeCourse)));
        SimulateCondition.NumSteps = 4;
        SimulateCondition.RhDecayFact = 0;
        SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval / 2000;
        SimulateCondition = RhPlusTransducinCascadeModel(SimulateCondition, Resp, NumResponses);
        temp = conv(SimulateCondition.MeanTrTimeCourse, WTFilt, 'full');
        pred2 = temp(1:SimulateCondition.EpochPts);
    end

    % plot
    figure(1);
    plot(tme, WTFilt);
    figure(2); clf
    plot(tme, Resp, 'k', tme, AveSingles.CSMAveSingle(1:length(tme)), 'r', tme, pred, 'k', tme, pred2, 'r');
    set(2, 'DefaultAxesFontSize', 11)
    legend('GRK+/-', 'CSM', 'GRK+/- pred', 'CSM pred', 'Location', 'East');
    legend boxoff
    set(2, 'DefaultAxesFontSize', 14)
    xlabel('ms');
    ylabel('pA');
    xlim([0 2000]);
    Likelihood(tcomp) = normpdf((max(pred2) - AveSingles.CSMPeakAmp) / AveSingles.CSMPeakAmpSEM);
    fprintf(1, '%d\n', Likelihood(tcomp));
    pause(1);
end 

figure(1);
loglog(TransCompression, Likelihood / max(Likelihood));
xlabel('compression');
ylabel('likelihood');

%%
% estimate trasduction cascade in two cases assuming same Rh kinetics
RespLength = length(WTAveNormSingle);

Resp = WTAveNormSingle(1:RespLength);   % base response to compute R*->I filter
tme = 1:RespLength;                          % time in ms
tcon = 320;                                      % R*(t) decay time constant - units ms
FreqCutoff = 10;                                

% compute R*->I filter for first response
Rh = exp(-tme/tcon);

filtFFT = fft(Resp) ./ fft(Rh);
filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
filt = real(ifft(filtFFT));
pred = real(ifft(fft(Rh) .* fft(filt)));        % sanity check - essentially checks frequency cutoff

% compute R*->I filter for second response
RespLength = length(GCAPAveNormSingle36C);
Resp = GCAPAveNormSingle36C(1:RespLength);   % base response to compute R*->I filter
tme = 1:RespLength;                          % time in ms
Rh = exp(-tme/tcon);
filtFFT = fft(Resp) ./ fft(Rh);
filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
filt2 = real(ifft(filtFFT));
pred2 = real(ifft(fft(Rh) .* fft(filt2)));        % sanity check - essentially checks frequency cutoff

% plot
figure(1);clf
plot(tme(1:length(filt)), filt, tme, filt2, 'r');
figure(2);clf
plot(tme(1:length(filt)), pred, tme, pred2, 'r');



%% 
%----------------------------------------------------------------------------------
% ARCHIVE
%----------------------------------------------------------------------------------

%%
%----------------------------------------------------------------------------------
% MODELS FOR TIME-DEPENDENT VARIANCE BASED ON NONLINEAR TRANSDUCTION CASCADES
%----------------------------------------------------------------------------------

%---------------------------------------------------------------------
% saturation of transducin activity

NumResponses = 	500;
NumSteps = 10;
RhTCon = 50;
TranTCon = 300;
RhDecayFact = 1;
TransRate = 0.3;
TransCompression = [0];
clear MeanTAct CVArea

if (1)
    ExpMean = WTAveNormSingle36C;
    ExpVar = WTVarNormSingle36C;
else
    ExpMean = GCAPAveNormSingle36C;
    ExpVar = GCAPVarNormSingle36C;
end

RespLength = length(ExpMean);
Resp = ExpMean(1:RespLength);   % base response to compute R*->I filter
tme = 1:RespLength;                          % time in ms
FreqCutoff = 10;                                

for cnt = 1:length(TransCompression)

    % general setup stuff
    clear SimulateCondition;
    SimulateCondition.SamplingInterval = 1;         % ms
    SimulateCondition.EpochPts = length(WTFilt);
    SimulateCondition.EpochNumbers(1:NumResponses) = 0;
    SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
    SimulateCondition.RhDecayFact = RhDecayFact;
    SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / RhTCon;
    SimulateCondition.NumSteps = NumSteps;
    SimulateCondition.ResponseCompression = 0;
    SimulateCondition.RhodopsinCompression = 0;
    SimulateCondition.DeterministicRhModel = 1;
    SimulateCondition.DeterministicTrModel = 0;
    SimulateCondition.BaseTransRate = TransRate;                  % cannot exceed 1 - decrease sampling interval to get higher rate
    SimulateCondition.TransCompression = TransCompression(cnt);
    SimulateCondition.TransDecayTimeConst = TranTCon / SimulateCondition.SamplingInterval;

    % multi-step shutoff
    SimulateCondition = RhPlusTransducinCascadeModel(SimulateCondition, Resp, NumResponses);

    % compute transduction cascade filter
    filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanTrTimeCourse);
    filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
    WTFilt = real(ifft(filtFFT));    

    % multi-step shutoff
    SimulateCondition = RhPlusTransducinCascadeModel(SimulateCondition, WTFilt, NumResponses);

    SimulateCondition.AverageResponse = mean(SimulateCondition.EpochData);
    SimulateCondition.VarianceResponse = var(SimulateCondition.EpochData);

    fprintf('CV = %d\n', std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData')));
    fprintf('mean number T activated = %d\n', mean(SimulateCondition.TransCount));
    CVArea(cnt) = std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData'));
    MeanTAct(cnt) = mean(SimulateCondition.TransCount);

    [MaxVal, MeanTPeak] = max(SimulateCondition.AverageResponse);
    avesqsim = (SimulateCondition.AverageResponse / MaxVal).^2;
    avesqsim = avesqsim / 10;
    varsim = SimulateCondition.VarianceResponse / MaxVal^2;
    tme2 = 1:length(varsim);
    tme2 = tme2*SimulateCondition.SamplingInterval;

    figure(2);clf
    plot(tme2, avesqsim, 'k', tme2, varsim, 'k--');
    axis([0 1500 -.01 .2]);
    xlabel('time (ms)');
    ylabel('normalize amp');
    pause(1);
    
end

figure(3);
semilogx(TransCompression, CVArea, 'o');
xlabel('ave # T activated');
ylabel('CV Area');
axis tight

%%
% compare predictions of models as abruptness of rhodopsin inactivation
% varied
% use PCA and every response
% has some problems 

RhDecayFact = [1];
NumResponses = 	500;       % with 5000 resps simulation std ~1%, about 10% std in likelihood
NumSteps = [6 8 10];
RhTCon = [300];
Resp = WTAveNormSingle30C;
Target = WTSingles;
FreqCutoff = 20;                                
Bootstraps = 200;
Verbose = 1;
NumLoops = 1;

clear TPeakRatio VarWidth Likelihood CV TempExpVar TempExpMean;

% step 1: identify space to describe responses and measure distribution of
% experimental responses in this space

% PCA of measured responses
resp = Target.Responses' - repmat(mean(Target.Responses'), size(Target.Responses, 2), 1);
%CovMeasuredSingles = cov(resp) - cov(Target.Failures');
CovMeasuredSingles = cov(resp);
opts.disp = 0;
[EigVec, EigVal] = eigs(CovMeasuredSingles, 3, 'LM', opts);
noise = Target.Failures';
SinglesProjections = resp * EigVec;
NoiseProjections = noise * EigVec;
CV2 = sqrt(var(sum(Target.Responses)) - var(sum(Target.Failures))) / mean(sum(Target.Responses));
ExpVar = var(Target.Responses') - var(Target.Failures');
ExpMean = mean(Target.Responses');

for loop = 1:NumLoops
    
    for perm = 1:Bootstraps
        Indices = randperm(size(resp, 1));
        TempExpMean(perm, :) = mean(SinglesProjections(Indices(1:floor(length(Indices)/4)), :));
        TempExpVar(perm, :) = var(SinglesProjections(Indices(1:floor(length(Indices)/4)), :)) - var(NoiseProjections(Indices(1:floor(length(Indices)/4)), :));
    end
    MeanExpVar = mean(TempExpVar);
    SDExpVar = std(TempExpVar);

    for step = 1:length(NumSteps)

        for tcon = 1:length(RhTCon)

            for decay = 1:length(RhDecayFact)

                RespLength = length(Resp);
                tme = 1:RespLength;                          % time in ms

                % run stochastic simulation to get mean rhodopsin time course
                clear SimulateCondition;
                SimulateCondition.NumSteps = NumSteps(step);
                SimulateCondition.RhDecayFact = RhDecayFact(decay);
                SimulateCondition.SamplingInterval = 1;         % ms
                SimulateCondition.EpochPts = length(Resp);
                SimulateCondition.EpochNumbers(1:NumResponses) = 0;
                SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
                SimulateCondition.RhodopsinCompression = 0;
                SimulateCondition.ResponseCompression = 0;
                SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps(step) / RhTCon(tcon);
                SimulateCondition.NumPtsToPeak = 100;
                SimulateCondition = RhCascadeModel(SimulateCondition, Resp, NumResponses);

                % compute transduction cascade filter
                filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanRhTimeCourse);
                filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
                WTFilt = real(ifft(filtFFT));    

                % generate responses
                clear SimulateCondition;
                SimulateCondition.SamplingInterval = 1;         % ms
                SimulateCondition.EpochPts = length(WTFilt);
                SimulateCondition.EpochNumbers(1:NumResponses) = 0;
                SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
                SimulateCondition.RhDecayFact = RhDecayFact(decay);
                SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps(step) / RhTCon(tcon);
                SimulateCondition.NumSteps = NumSteps(step);
                SimulateCondition.ResponseCompression = 0;
                SimulateCondition.RhodopsinCompression = 0;
                SimulateCondition.NumPtsToPeak = 100;

                % multi-step shutoff
                SimulateCondition = RhCascadeModel(SimulateCondition, WTFilt, NumResponses);

                % compare modeled and measured mean squared and variance
                clear SimProjections SimProb respLength
                respLength = min([size(resp, 2) size(SimulateCondition.NormSingles, 1)]);
                resp = SimulateCondition.NormSingles' - repmat(mean(SimulateCondition.NormSingles'), NumResponses, 1);
                SimProjections = resp(:, 1:respLength) * EigVec(1:respLength, :);
                SimVar = var(SimProjections);
                SimProb = abs(SimVar - MeanExpVar) ./ SDExpVar;

                MeanSingle = mean(SimulateCondition.NormSingles');
                SimEigVec = EigVec(1:respLength, :);
                resp = SimProjections * SimEigVec' + repmat(MeanSingle(1:respLength), NumResponses, 1);
                figure(2); clf
                plot(mean(resp).^2/10);
                hold on;
                plot(var(resp));
                plot(mean(SimulateCondition.NormSingles').^2/10, 'r');
                plot(var(SimulateCondition.NormSingles'), 'r');
                plot(ExpVar, 'k');
                plot(ExpMean.^2/10, 'k');
                
                CV = std(sum(SimulateCondition.NormSingles)) / mean(sum(SimulateCondition.NormSingles));
                
                Likelihood(decay, tcon, step) = prod(normpdf(SimProb));
                if (Verbose)
                    fprintf(1, '\tCV = %d (%d) Likelihood = %d\n', CV, CV2, Likelihood(decay, tcon, step));
                    SimProb
                end
                pause(0.1)

            end

        end

    end

    [MaxVal, MaxLoc] = max(Likelihood(:));
    MaxLoc = MaxLoc-1;
    Index1 = floor(MaxLoc / (length(RhDecayFact) * length(RhTCon)));
    Index2 = floor((MaxLoc - Index1 * (length(RhDecayFact) * length(RhTCon))) / length(RhDecayFact));
    Index3 = MaxLoc - Index1 * (length(RhDecayFact) * length(RhTCon)) - Index2 * length(RhDecayFact);

    fprintf(1, '%d %d %d\n', Index1+1, Index2+1, Index3+1);
    fprintf(1, '%d %d %d\n', NumSteps(Index1+1), RhTCon(Index2+1), RhDecayFact(Index3+1));

end

%%
%----------------------------------------------------------------------------------
% MODELS FOR TIME-DEPENDENT VARIANCE BASED ON LINEAR TRANSDUCTION CASCADES
%----------------------------------------------------------------------------------

%%
% predict time-dependent variance in two cases assuming same Rh kinetics
% but different (linear) transduction cascade kinetics
%
% do this across a range of inactivation kinetics for rhodopsin to see
% which provides best explanation of data

RhShutoffFact = [225]      % effective time constants for rhodopsin shutoff
NumResponses = 	1000;
NumSteps = 7;
RhDecayFact = 1;

clear WTTPeakRatio GCAPTPeakRatio;

for cnt = 1:length(RhShutoffFact)

    tcon = RhShutoffFact(cnt);                                      % R*(t) decay time constant - units ms
    FreqCutoff = 10;                                

    % compute cascade filter for first response
    RespLength = length(WTAveNormSingle36C);
    Resp = WTAveNormSingle36C(1:RespLength);   % base response to compute R*->I filter
    tme = 1:RespLength;                          % time in ms

    % run stochastic simulation to get mean rhodopsin time course
    clear SimulateCondition;
    SimulateCondition.NumSteps = NumSteps;
    SimulateCondition.RhDecayFact = RhDecayFact;
    SimulateCondition.SamplingInterval = 1;         % ms
    SimulateCondition.EpochPts = length(Resp);
    SimulateCondition.EpochNumbers(1:NumResponses) = 0;
    SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
    SimulateCondition.RhodopsinCompression = 0;
    SimulateCondition.ResponseCompression = 0;
    SimulateCondition.NumPtsToPeak = 100;
    SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / tcon;
    SimulateCondition = RhCascadeModel(SimulateCondition, Resp, NumResponses);

    filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanRhTimeCourse);
    filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
    WTFilt = real(ifft(filtFFT));

    % compute cascade filter for second response
    RespLength = length(GCAPAveNormSingle36C);
    Resp = GCAPAveNormSingle36C(1:RespLength);   % base response to compute R*->I filter
    tme = 1:RespLength;                          % time in ms

    % run stochastic simulation to get mean rhodopsin time course
    clear SimulateCondition;
    SimulateCondition.NumSteps = NumSteps;
    SimulateCondition.RhDecayFact = RhDecayFact;
    SimulateCondition.SamplingInterval = 1;         % ms
    SimulateCondition.EpochPts = length(Resp);
    SimulateCondition.EpochNumbers(1:NumResponses) = 0;
    SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
    SimulateCondition.RhodopsinCompression = 0;
    SimulateCondition.ResponseCompression = 0;
    SimulateCondition.NumPtsToPeak = 100;
    SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / tcon;
    SimulateCondition = RhCascadeModel(SimulateCondition, Resp, NumResponses);

    % compute R*->I filter for second response
    filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanRhTimeCourse);
    filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
    GCAPFilt = real(ifft(filtFFT));

    % plot
    figure(1);clf
    plot(tme(1:length(WTFilt)), WTFilt, 'k', tme, GCAPFilt, 'r');

   %----------------------------------------------------------------------------------
   % compare wt and GCAP knockout responses using linear filters calculated
    % assuming identical average time course of rhodopsin activity and linear
    % transduction cascade (i.e. fit to mean response).

    % general setup stuff
    clear SimulateCondition;
    SimulateCondition.NumSteps = NumSteps;
    SimulateCondition.RhDecayFact = RhDecayFact;
    SimulateCondition.SamplingInterval = 1;         % ms
    SimulateCondition.EpochPts = length(WTFilt);
    SimulateCondition.EpochNumbers(1:NumResponses) = 0;
    SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
    SimulateCondition.RhodopsinCompression = 0;
    SimulateCondition.ResponseCompression = 0;
    SimulateCondition.NumPtsToPeak = 100;

    % multi-step shutoff
    ModelType = 'MultiStep';
    SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / tcon;
    SimulateCondition = RhCascadeModel(SimulateCondition, WTFilt, NumResponses);

    % general setup stuff
    clear GCAPSimulateCondition;
    GCAPSimulateCondition.SamplingInterval = 1;         % ms
    GCAPSimulateCondition.EpochPts = length(GCAPFilt);
    GCAPSimulateCondition.EpochNumbers(1:NumResponses) = 0;
    GCAPSimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
    GCAPSimulateCondition.NumSteps = NumSteps;
    GCAPSimulateCondition.RhDecayFact = RhDecayFact;
    GCAPSimulateCondition.SamplingInterval = 1;         % ms
    GCAPSimulateCondition.EpochPts = length(GCAPFilt);
    GCAPSimulateCondition.EpochNumbers(1:NumResponses) = 0;
    GCAPSimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
    GCAPSimulateCondition.RhodopsinCompression = 0;
    GCAPSimulateCondition.ResponseCompression = 0;
    GCAPSimulateCondition.NumPtsToPeak = 100;
    
    GCAPSimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / tcon;

    % multi-step shutoff
    GCAPSimulateCondition = RhCascadeModel(GCAPSimulateCondition, GCAPFilt, NumResponses);

    % compare modeled and measured mean squared and variance
    fprintf('\nWT CV = %d\n', std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData')));
    fprintf('GCAP CV = %d\n', std(sum(GCAPSimulateCondition.EpochData')) / mean(sum(GCAPSimulateCondition.EpochData')));

    SimulateCondition.AverageResponse = mean(SimulateCondition.EpochData);
    SimulateCondition.VarianceResponse = var(SimulateCondition.EpochData);
    GCAPSimulateCondition.AverageResponse = mean(GCAPSimulateCondition.EpochData);
    GCAPSimulateCondition.VarianceResponse = var(GCAPSimulateCondition.EpochData);

    [MaxVal, WTMeanTPeak] = max(SimulateCondition.AverageResponse);
    avesqsim = (SimulateCondition.AverageResponse / MaxVal).^2;
    avesqsim = avesqsim / 10;
    varsim = SimulateCondition.VarianceResponse / MaxVal^2;
    tme2 = 1:length(varsim);
    tme2 = tme2*SimulateCondition.SamplingInterval;

    figure(2);clf
    plot(tme2, avesqsim, 'r', tme2, varsim, 'r--');
    axis([0 1500 -.01 .2]);
    xlabel('normalized time');
    ylabel('normalize amp');
    hold on
    [MaxVal, MeanTPeak] = max(WTAveNormSingle36C);
    plot(tme2, WTAveNormSingle36C.^2 / MaxVal^2 / 10, 'k--', tme2, WTVarNormSingle36C / MaxVal^2, 'k--');

    [MaxVal, GCAPMeanTPeak] = max(GCAPSimulateCondition.AverageResponse);
    avesqsim = (GCAPSimulateCondition.AverageResponse / MaxVal).^2;
    avesqsim = avesqsim / 10;
    varsim = GCAPSimulateCondition.VarianceResponse / MaxVal^2;
    tme2 = 1:length(varsim);
    tme2 = tme2*GCAPSimulateCondition.SamplingInterval;

    figure(3); clf
    plot(tme2, avesqsim, 'r', tme2, varsim, 'r');
    hold on
    axis([0 1500 -.01 .2]);
    set(2, 'DefaultAxesFontSize', 14)
    xlabel('time (ms)');
    ylabel('normalize amp');
    [MaxVal, MeanTPeak] = max(GCAPAveNormSingle36C);
    plot(tme2, GCAPAveNormSingle36C.^2 / MaxVal^2 / 10, 'k--', tme2, GCAPVarNormSingle36C / MaxVal^2, 'k--');
    set(2, 'DefaultAxesFontSize', 11)
    legend('I wt', 'var wt', 'exp I wt', 'exp var wt', 'I GCAP', 'var GCAP', 'exp I GCAP', 'exp var GCAP');
    legend boxoff

    [MaxVal, WTVarTPeak] = max(SimulateCondition.VarianceResponse);
    [MaxVal, GCAPVarTPeak] = max(GCAPSimulateCondition.VarianceResponse);
    WTTPeakRatio(cnt) = WTVarTPeak / WTMeanTPeak;
    GCAPTPeakRatio(cnt) = GCAPVarTPeak / GCAPMeanTPeak;
    pause(1);
    
end

[MaxVal, WTVarTPeak] = max(WTVarNormSingle36C);
[MaxVal, WTMeanTPeak] = max(WTAveNormSingle36C);
[MaxVal, GCAPVarTPeak] = max(GCAPVarNormSingle36C);
[MaxVal, GCAPMeanTPeak] = max(GCAPAveNormSingle36C);
ExpWTTPeakRatio = WTVarTPeak / WTMeanTPeak;
ExpGCAPTPeakRatio = GCAPVarTPeak / GCAPMeanTPeak;

figure(4);
plot(RhShutoffFact, WTTPeakRatio, 'ok', RhShutoffFact, GCAPTPeakRatio, 'or');
xlabel('\tau_{Rh}')
ylabel('t_{peak, \sigma^2} / t_{peak, m}')
fprintf(1, '%d\n%d\n', ExpWTTPeakRatio, ExpGCAPTPeakRatio);
legend('WT', 'GCAP-/-')
legend boxoff

%----------------------------------------------------------------------------------
%%
% compare predictions of models as abruptness of rhodopsin inactivation
% varied

RhDecayFact = [1];
NumResponses = 	100;
NumSteps = 9;
RhTCon = 350;
Resp = WTAveNormSingle30C;

clear TPeakRatio VarWidth;

for cnt = 1:length(RhDecayFact)
    
    RespLength = length(Resp);
    tme = 1:RespLength;                          % time in ms
    FreqCutoff = 10;                                

    % run stochastic simulation to get mean rhodopsin time course
    clear SimulateCondition;
    SimulateCondition.NumSteps = NumSteps;
    SimulateCondition.RhDecayFact = RhDecayFact(cnt);
    SimulateCondition.SamplingInterval = 1;         % ms
    SimulateCondition.EpochPts = length(Resp);
    SimulateCondition.EpochNumbers(1:NumResponses) = 0;
    SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
    SimulateCondition.RhodopsinCompression = 0;
    SimulateCondition.ResponseCompression = 0;
    SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / RhTCon;
    SimulateCondition.NumPtsToPeak = 100;
    SimulateCondition = RhCascadeModel(SimulateCondition, Resp, NumResponses);
    
    % compute transduction cascade filter
    filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanRhTimeCourse);
    filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
    WTFilt = real(ifft(filtFFT));    
    
    % generate responses
    clear SimulateCondition;
    SimulateCondition.SamplingInterval = 1;         % ms
    SimulateCondition.EpochPts = length(WTFilt);
    SimulateCondition.EpochNumbers(1:NumResponses) = 0;
    SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
    SimulateCondition.RhDecayFact = RhDecayFact(cnt);
    SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / RhTCon;
    SimulateCondition.NumSteps = NumSteps;
    SimulateCondition.ResponseCompression = 0;
    SimulateCondition.RhodopsinCompression = 0;
    SimulateCondition.NumPtsToPeak = 100;

    % multi-step shutoff
    ModelType = 'MultiStep';
    SimulateCondition = RhCascadeModel(SimulateCondition, WTFilt, NumResponses);

    % compare modeled and measured mean squared and variance
    fprintf('\nWT CV = %d\n', std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData')));

    SimulateCondition.AverageResponse = mean(SimulateCondition.EpochData);
    SimulateCondition.VarianceResponse = var(SimulateCondition.EpochData);

    [MaxVal, MeanTPeak] = max(SimulateCondition.AverageResponse);
    avesqsim = (SimulateCondition.AverageResponse / MaxVal).^2;
    avesqsim = avesqsim / 10;
    varsim = SimulateCondition.VarianceResponse / MaxVal^2;
    tme2 = 1:length(varsim);
    tme2 = tme2*SimulateCondition.SamplingInterval;

    figure(2);clf
    plot(tme2, avesqsim, 'k', tme2, varsim, 'k');
    axis([0 1500 -.01 .2]);
    xlabel('time (ms)');
    ylabel('normalize amp');
    hold on
    [MaxVal, MeanTPeak] = max(WTAveNormSingle36C);
    plot(tme2, WTAveNormSingle36C.^2 / MaxVal^2 / 10, 'k--', tme2, WTVarNormSingle36C / MaxVal^2, 'k--');
    set(2, 'DefaultAxesFontSize', 11)
    legend('I wt', 'var wt', 'exp I wt', 'exp var wt')
    legend boxoff
    pause(1)

    [MaxVal, VarTPeak] = max(SimulateCondition.VarianceResponse);
    TPeakRatio(cnt) = VarTPeak / MeanTPeak;
    Indices = find(SimulateCondition.VarianceResponse > 0.5*MaxVal);
    FirstTime = Indices(min(find(Indices < VarTPeak)));
    SecondTime = Indices(max(find(Indices > VarTPeak)));
    VarWidth(cnt) = SecondTime - FirstTime;
    
end

figure(3);clf
plot(RhDecayFact, TPeakRatio, 'ko');
xlabel('Rh* fractional decay')
ylabel('t_{peak, \sigma^2} / t_{peak, m}')
figure(4);clf
plot(RhDecayFact, VarWidth, 'ko');
xlabel('Rh* fractional decay')
ylabel('variance half width')

%----------------------------------------------------------------------------------
%%
% model including stochastic transducin activation

NumResponses = 	500;
NumSteps = 8;
RhTCon = 250;
TranTCon = 300;
RhDecayFact = 1;
TransRate = [0.08 0.1 0.2 0.4];

if (1)
    ExpMean = WTAveNormSingle36C;
    ExpVar = WTVarNormSingle36C;
else
    ExpMean = GCAPAveNormSingle36C;
    ExpVar = GCAPVarNormSingle36C;
end

RespLength = length(ExpMean);
Resp = ExpMean(1:RespLength);   % base response to compute R*->I filter
tme = 1:RespLength;                          % time in ms
FreqCutoff = 10;                                
clear MeanTAct CVArea

for cnt = 1:length(TransRate)

    % general setup stuff
    clear SimulateCondition;
    SimulateCondition.SamplingInterval = 1;         % ms
    SimulateCondition.EpochPts = length(ExpMean);
    SimulateCondition.EpochNumbers(1:NumResponses) = 0;
    SimulateCondition.ExcludeEpochs(1:NumResponses) = 0;
    SimulateCondition.RhDecayFact = RhDecayFact;
    SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / RhTCon;
    SimulateCondition.NumSteps = NumSteps;
    SimulateCondition.ResponseCompression = 0;
    SimulateCondition.RhodopsinCompression = 0;
    SimulateCondition.BaseTransRate = TransRate(cnt) * SimulateCondition.SamplingInterval;          % cannot exceed 1 - decrease sampling interval to get higher rate
    SimulateCondition.TransCompression = 0;
    SimulateCondition.TransDecayTimeConst = TranTCon / SimulateCondition.SamplingInterval;
    SimulateCondition.DeterministicRhModel = 0;
    
    % multi-step shutoff
    SimulateCondition = RhPlusTransducinCascadeModel(SimulateCondition, Resp, NumResponses);

    % compute transduction cascade filter
    filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanTrTimeCourse);
    filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
    WTFilt = real(ifft(filtFFT));    

    % multi-step shutoff
    SimulateCondition = RhPlusTransducinCascadeModel(SimulateCondition, WTFilt, NumResponses);

    SimulateCondition.AverageResponse = mean(SimulateCondition.EpochData);
    SimulateCondition.VarianceResponse = var(SimulateCondition.EpochData);

    fprintf('CV = %d\n', std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData')));
    fprintf('mean number T activated = %d\n', mean(SimulateCondition.TransCount));
    CVArea(cnt) = std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData'));
    MeanTAct(cnt) = mean(SimulateCondition.TransCount);

    [MaxVal, MeanTPeak] = max(SimulateCondition.AverageResponse);
    avesqsim = (SimulateCondition.AverageResponse / MaxVal).^2;
    avesqsim = avesqsim / 10;
    varsim = SimulateCondition.VarianceResponse / MaxVal^2;
    tme2 = 1:length(varsim);
    tme2 = tme2*SimulateCondition.SamplingInterval;

    figure(2);clf
    plot(tme2, avesqsim, 'k', tme2, varsim, 'k');
    axis([0 1500 -.01 .2]);
    xlabel('time (ms)');
    ylabel('normalize amp');
    hold on
    plot(tme2, ExpMean.^2 / 10, 'k--', tme2, ExpVar, 'k--');
    set(2, 'DefaultAxesFontSize', 11)
    legend('I wt', 'var wt', 'exp I wt', 'exp var wt')
    legend boxoff
    pause(1);

end

figure(3);clf
semilogx(MeanTAct, CVArea, 'o');
xlabel('ave # T activated');
ylabel('CV Area');
axis tight
