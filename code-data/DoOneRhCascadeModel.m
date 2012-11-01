function SimReturns =  DoOneRhCascadeModel(NumSteps, RhDecayFact, RhTCon, SimParams)        

    RespLength = length(SimParams.Resp);
    tme = 1:RespLength;                          % time in ms

    % run stochastic simulation to get mean rhodopsin time course
    clear SimulateCondition;
    SimulateCondition.NumSteps = NumSteps;
    SimulateCondition.RhDecayFact = RhDecayFact;
    SimulateCondition.SamplingInterval = 1;         % ms
    SimulateCondition.EpochPts = RespLength;
    SimulateCondition.EpochNumbers(1:SimParams.NumResponses) = 0;
    SimulateCondition.ExcludeEpochs(1:SimParams.NumResponses) = 0;
    SimulateCondition.RhodopsinCompression = 0;
    SimulateCondition.ResponseCompression = 0;
    SimulateCondition.InitialShutoffRate = SimulateCondition.SamplingInterval * NumSteps / RhTCon;
    SimulateCondition.NumPtsToPeak = 100;

    SimulateCondition = RhCascadeModel(SimulateCondition, SimParams.Resp, SimParams.NumResponses);

    % compute transduction cascade filter
    filtFFT = fft(SimParams.Resp) ./ fft(SimulateCondition.MeanRhTimeCourse);
    filtFFT(SimParams.FreqCutoff+1:length(filtFFT)-SimParams.FreqCutoff) = 0;
    WTFilt = real(ifft(filtFFT));    

    % generate responses
    SimulateCondition = RhCascadeModel(SimulateCondition, WTFilt, SimParams.NumResponses);

    % compare modeled and measured mean squared and variance
    SimReturns.CV = std(sum(SimulateCondition.EpochData')) / mean(sum(SimulateCondition.EpochData'));

    SimulateCondition.AverageResponse = mean(SimulateCondition.EpochData);
    SimulateCondition.VarianceResponse = var(SimulateCondition.EpochData);

    [MaxVal, MeanTPeak] = max(SimulateCondition.AverageResponse);
    [MaxVal, VarTPeak] = max(SimulateCondition.VarianceResponse);
    SimReturns.TPeakRatio = VarTPeak / MeanTPeak;
    Indices = find(SimulateCondition.VarianceResponse > 0.5*MaxVal);
    FirstTime = Indices(min(find(Indices < VarTPeak)));
    SecondTime = Indices(max(find(Indices > VarTPeak)));
    SimReturns.VarWidth = SecondTime - FirstTime;

end