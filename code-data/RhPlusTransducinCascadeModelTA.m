function ReturnedCondition = RhPlusTransducinCascadeModelTA(Condition, Filter, NumResponses)
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

% put specified filter into frequency domain
ReturnedCondition = Condition;
if ~isempty(Filter)
    filt = fft(Filter);
    ReturnedCondition.Filter = Filter;
else
    % Run the simulation, then calculate the filter, then filter epoch data
end

if (Condition.DeterministicTrModel)
    % filter for transducin activity
    tme = 1:Condition.EpochPts;
    TrFilt = exp(-tme / Condition.TransDecayTimeConst);
    TrFilt = TrFilt / sum(TrFilt);
    TrFilt = fft(TrFilt);
end

if (Condition.DeterministicRhModel)
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
        if (resp == 1)
            MeanRhTimeCourse = RhTimeCourse;
        else
            MeanRhTimeCourse = MeanRhTimeCourse + RhTimeCourse;
        end

    end
    MeanRhTimeCourse = MeanRhTimeCourse / NumResponses;
end

if (Condition.DeterministicRhModel && Condition.DeterministicTrModel)
   
    temp = conv(TrFilt, MeanRhTimeCourse, 'full');
    TransducinActivity = temp(1:Condition.EpochPts);
    
    % compress transducin activity if desired
    if (Condition.TransCompression > 0)
        TransducinActivity = TransducinActivity ./ (1 + Condition.TransCompression * cumsum(TransducinActivity));
    end
    if (Condition.TransRateCompression > 0)
        TransducinActivity = TransducinActivity ./ (1 + Condition.TransRateCompression * TransducinActivity);
    end
    
    ReturnedCondition.MeanStochTimeCourse = TransducinActivity;
            
	% convolve transducin time course with linear filter to generate modeled response
	TransducinActivity = fft(TransducinActivity);
    ReturnedCondition.EpochData = real(ifft(TransducinActivity .* filt));
else

    numexamples = floor(log2(NumResponses));
    randexamples = randi(NumResponses,numexamples,1);
    randexamples = unique(randexamples);
    ReturnedCondition.TExamples = zeros(length(randexamples),Condition.EpochPts);
    ReturnedCondition.IExamples = zeros(length(randexamples),Condition.EpochPts);
    
    % Pre definining EpochData matrix speeds up by 31X for 1000 Epochs
    ReturnedCondition.EpochData = zeros(NumResponses,Condition.EpochPts);
    ReturnedCondition.TransCount = zeros(1,NumResponses);

    %tic
    % generate series of responses
    for resp = 1:NumResponses
        
        % initial settings for rhodopsin activity
        RhTimeCourse(1:Condition.EpochPts) = 0;
        TransducinActivity(1:Condition.EpochPts) = 0;
        CurrentStep = 0;
        CurrentCatalyticActivity = 1;
        ShutoffRate = Condition.InitialShutoffRate;
        TransCount = 0;
        
        % for each time point decide whether shutoff reaction has occurred.
        % if it has, update catalytic activity and shutoff rate.
        for cnt=1:Condition.EpochPts
            if (Condition.DeterministicRhModel)
                CurrentCatalyticActivity = MeanRhTimeCourse(cnt);
            end
            % did we activate transducin in this time step?
            RhTimeCourse(cnt) = CurrentCatalyticActivity;
            if (Condition.DeterministicTrModel == 0)
                if (rand(1) < (CurrentCatalyticActivity * Condition.BaseTransRate))
                    TransducinLifetime = round(exprnd(Condition.TransDecayTimeConst));
                    if ( (cnt + TransducinLifetime) < Condition.EpochPts)
                        TransducinActivity(cnt:cnt+TransducinLifetime) = TransducinActivity(cnt:cnt+TransducinLifetime) + 1;
                    else
                        TransducinActivity(cnt:Condition.EpochPts) = TransducinActivity(cnt:Condition.EpochPts) + 1;
                    end
                    TransCount = TransCount + 1;
                end
            end
            % generate random number between 0 and 1 (uniform dist) and compare to shutoff rate
            if (Condition.DeterministicRhModel == 0)
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
        end
        
        % generate transductin activity if not stochastic
        if (Condition.DeterministicTrModel)
            RhTimeCourse = fft(RhTimeCourse);
            TransducinActivity = real(ifft(RhTimeCourse .* TrFilt));
            %           temp = conv(TrFilt, RhTimeCourse, 'full');
            %          TransducinActivity = temp(1:Condition.EpochPts);
        end

        % compress transducin activity if desired
        if (Condition.TransCompression > 0)
            TransducinActivity = TransducinActivity ./ (1 + Condition.TransCompression * cumsum(TransducinActivity));
        end
        if (Condition.TransRateCompression > 0)
            TransducinActivity = TransducinActivity ./ (1 + Condition.TransRateCompression * TransducinActivity);
        end
    
        if (resp == 1)
            ReturnedCondition.MeanTrTimeCourse = TransducinActivity;
        else
            ReturnedCondition.MeanTrTimeCourse = ReturnedCondition.MeanTrTimeCourse + TransducinActivity;
        end
        if sum(resp == randexamples)
            ReturnedCondition.TExamples(find(resp==randexamples),:) = TransducinActivity;
        end
        
        
        ReturnedCondition.EpochData(resp, :) = TransducinActivity;
        
        ReturnedCondition.TransCount(resp) = TransCount;
    end
    ReturnedCondition.MeanStochTimeCourse = ReturnedCondition.MeanTrTimeCourse / NumResponses;
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
end


% Push resampling back to cascade model to avoid using it for feature
% comparison (CV, ttpk, varwidth)
% NormSingles = resample(ReturnedCondition.EpochData', Condition.NumPtsToPeak, PeakTime) / Peakamp;
% ReturnedCondition.NormSingles = NormSingles;
% ReturnedCondition.RhExamples = resample(ReturnedCondition.RhExamples', Condition.NumPtsToPeak, PeakTime);
% ReturnedCondition.IExamples = resample(ReturnedCondition.IExamples', Condition.NumPtsToPeak, PeakTime) / Peakamp;
