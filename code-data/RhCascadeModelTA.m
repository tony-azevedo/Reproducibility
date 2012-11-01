function ReturnedCondition = RhCascadeModelTA(Condition, Filter, NumResponses)
% Calculate series of dim flash responses for multi-step shutoff
% model for rhodopsin activity.  First generate time course of 
% rhodopsin activity assuming rhodopsin shutoff is described as a
% series of independent first order reactions.  Each shutoff 
% reaction is assumed to produce an equal decrease in rhodopsin's 
% catalytic activity.  The rate constant for each shutoff reaction
% increases as well, such that each reaction controls an equal amount
% of rhodopsin's activity
%
% Created 7/01 FMR

%tic
ReturnedCondition = Condition;
ReturnedCondition.EpochData = zeros(NumResponses,Condition.EpochPts);

numexamples = floor(log2(NumResponses));
randexamples = randi(NumResponses,numexamples,1);
randexamples = unique(randexamples);
ReturnedCondition.RhExamples = zeros(length(randexamples),Condition.EpochPts);
ReturnedCondition.IExamples = zeros(length(randexamples),Condition.EpochPts);

% generate series of responses
for resp = 1:NumResponses

	% initial settings for rhodopsin activity
	RhTimeCourse(1:Condition.EpochPts) = 0;
    CurrentStep = 0;
    CurrentCatalyticActivity = 1;
    ShutoffRate = Condition.InitialShutoffRate;
    
    % for each time point decide whether shutoff reaction has occurred.
    % if it has, update catalytic activity and shutoff rate.
    for cnt=1:Condition.EpochPts
        RhTimeCourse(cnt) = CurrentCatalyticActivity;
        % generate random number between 0 and 1 (uniform dist) and compare to shutoff rate
        if (rand(1) < ShutoffRate)
            ShutoffRate = Condition.InitialShutoffRate * (1 - CurrentStep * Condition.RhDecayFact / Condition.RhSteps);
            CurrentStep = CurrentStep + 1;
            CurrentCatalyticActivity = 1 - CurrentStep * Condition.RhDecayFact / Condition.RhSteps;
        end
        % are we done yet?
        if (CurrentStep == Condition.RhSteps)
            break;
        end
    end

    % compress rhodopsin activity if desired  
    if (Condition.RhodopsinCompression > 0)
        RhTimeCourse = RhTimeCourse ./ (1 + Condition.RhodopsinCompression * cumsum(RhTimeCourse));
    end
        
    if (resp == 1)
        Condition.MeanRhTimeCourse = RhTimeCourse;
    else            
        Condition.MeanRhTimeCourse = Condition.MeanRhTimeCourse + RhTimeCourse;
        ReturnedCondition.EpochData(resp, :) = RhTimeCourse;
        if sum(resp == randexamples)
            ReturnedCondition.RhExamples(find(resp==randexamples),:) = RhTimeCourse;
        end
    end
end

ReturnedCondition.MeanStochTimeCourse = Condition.MeanRhTimeCourse / NumResponses;

% convolve rhodopsin time course with linear filter to generate modeled
% response
if isempty(Filter)
    % Ifft the mean stochastic aspect of the SimFunc
    filtFFT = fft(ReturnedCondition.Resp) ./ fft(ReturnedCondition.MeanStochTimeCourse);
    filtFFT(ReturnedCondition.FreqCutoff+1:length(filtFFT)-ReturnedCondition.FreqCutoff) = 0;
    ReturnedCondition.Filter = real(ifft(filtFFT));
else
    ReturnedCondition.Filter = Filter;
end

for resp = 1:NumResponses
    ReturnedCondition.EpochData(resp, :) = real(ifft(ReturnedCondition.EpochData(resp, :) .* filtFFT));
    % compress if desired
    if (Condition.ResponseCompression > 0)
        ReturnedCondition.EpochData(resp, :) = (1-exp(-ReturnedCondition.EpochData(resp, :)*Condition.ResponseCompression))/Condition.ResponseCompression;
    end
    if sum(resp == randexamples)
        ReturnedCondition.IExamples(find(resp==randexamples),:) = ReturnedCondition.EpochData(resp, :);
    end
end

[Peakamp, PeakTime] = max(mean(ReturnedCondition.EpochData));
ReturnedCondition.time = (1:ReturnedCondition.EpochPts)*ReturnedCondition.NumPtsToPeak/PeakTime;
    
% Push resampling back to cascade model to avoid using it for feature
% comparison (CV, ttpk, varwidth)
% NormSingles = resample(ReturnedCondition.EpochData', Condition.NumPtsToPeak, PeakTime) / Peakamp;
% ReturnedCondition.NormSingles = NormSingles;
% ReturnedCondition.RhExamples = resample(ReturnedCondition.RhExamples', Condition.NumPtsToPeak, PeakTime);
% ReturnedCondition.IExamples = resample(ReturnedCondition.IExamples', Condition.NumPtsToPeak, PeakTime) / Peakamp;

