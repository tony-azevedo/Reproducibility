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

function ReturnedCondition = RhCascadeModel(Condition, Filter, NumResponses)

% put specified filter into frequency domain
filt = fft(Filter);
ReturnedCondition = Condition;
ReturnedCondition.EpochData= ones(NumResponses, Condition.EpochPts);

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
			ShutoffRate = Condition.InitialShutoffRate * (1 - CurrentStep * Condition.RhDecayFact / Condition.NumSteps);
			CurrentCatalyticActivity = 1 - CurrentStep * Condition.RhDecayFact / Condition.NumSteps;
			CurrentStep = CurrentStep + 1;
		end
		% are we done yet?
		if (CurrentStep == Condition.NumSteps)
			break;
		end
    end	
    
    % compress rhodopsin activity if desired
    if (Condition.RhodopsinCompression > 0)
        RhTimeCourse = RhTimeCourse ./ (1 + Condition.RhodopsinCompression * cumsum(RhTimeCourse));
    end
    
    if (resp == 1)
        ReturnedCondition.MeanRhTimeCourse = RhTimeCourse;
    else
        ReturnedCondition.MeanRhTimeCourse = ReturnedCondition.MeanRhTimeCourse + RhTimeCourse;
    end
    % convolve rhodopsin time course with linear filter to generate modeled response
	RhTimeCourse = fft(RhTimeCourse);
	ReturnedCondition.EpochData(resp, :) = real(ifft(RhTimeCourse .* filt));
    
    % compress if desired
	if (Condition.ResponseCompression > 0)
    	ReturnedCondition.EpochData(resp, :) = (1-exp(-ReturnedCondition.EpochData(resp, :)*Condition.ResponseCompression))/Condition.ResponseCompression;
    end
    
end

[Peakamp, PeakTime] = max(mean(ReturnedCondition.EpochData));
ReturnedCondition.MeanRhTimeCourse = ReturnedCondition.MeanRhTimeCourse / NumResponses;
NormSingles = resample(ReturnedCondition.EpochData', Condition.NumPtsToPeak, PeakTime) / Peakamp;
ReturnedCondition.NormSingles = NormSingles;
