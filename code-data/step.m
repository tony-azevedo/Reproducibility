%%
% assuming linear transduction cascade
Resp = AveSingles.RKAveSingle;   % base response to compute R*->I filter
temp = [1:6*length(Resp)]*0;
Resp = [Resp temp];
RespLength = length(Resp);
tme = 1:RespLength;                          % time in ms
tcon1 = 580;                                      % R*(t) decay time constant for "base" response - units ms
FreqCutoff = 100;                                
NumSteps = 10;
RhDecayFact = 1;
NumResponses = 20;
TransCompression = [0 1e-3 2e-3 4e-3 6e-3 8e-3 1e-2];
TranTCon = 80;

for tcomp = 1:length(TransCompression)
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
    SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / tcon1;
    SimulateCondition.TransCompression = TransCompression(tcomp);
    SimulateCondition.DeterministicRhModel = 1;
    SimulateCondition.DeterministicTrModel = 1;
    SimulateCondition.BaseTransRate = TransRate;                  % cannot exceed 1 - decrease sampling interval to get higher rate
    SimulateCondition.TransDecayTimeConst = TranTCon / SimulateCondition.SamplingInterval;

    if (0)
        SimulateCondition = RhCascadeModel(SimulateCondition, Resp, NumResponses);
        filtFFT = fft(Resp) ./ fft(SimulateCondition.MeanRhTimeCourse);
        filtFFT(FreqCutoff+1:length(filtFFT)-FreqCutoff) = 0;
        WTFilt = real(ifft(filtFFT));
        pred = real(ifft(fft(WTFilt) .* fft(SimulateCondition.MeanRhTimeCourse)));
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
        pred2 = real(ifft(fft(WTFilt) .* fft(SimulateCondition.MeanTrTimeCourse)));
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

%%