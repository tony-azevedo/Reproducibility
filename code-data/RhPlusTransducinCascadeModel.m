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

function ReturnedCondition = RhPlusTransducinCascadeModel(Condition, Filter, NumResponses)

% put specified filter into frequency domain
filt = fft(Filter);
ReturnedCondition = Condition;

% filter for transducin activity
TrFilt(1:Condition.EpochPts) = 0;
tme = 1:Condition.EpochPts;
TrFilt = exp(-tme / Condition.TransDecayTimeConst);
TrFilt = TrFilt / sum(TrFilt);

if (Condition.ConvFlag == 0)
    TrFilt = fft(TrFilt);
end

ReturnedCondition.EpochData = ones(NumResponses, Condition.EpochPts);

if (Condition.DeterministicRhModel)
    for resp = 1:NumResponses
        RhTimeCourse = RhTrajectoryGenerator(Condition);
        if (resp == 1)
            MeanRhTimeCourse = RhTimeCourse;
        else
            MeanRhTimeCourse = MeanRhTimeCourse + RhTimeCourse;
        end

    end
    MeanRhTimeCourse = MeanRhTimeCourse / NumResponses;
end

if (Condition.DeterministicRhModel && Condition.DeterministicTrModel)

    if (Condition.ConvFlag)
        temp = conv(TrFilt, MeanRhTimeCourse, 'full');
        TransducinActivity = temp(1:Condition.EpochPts);
    else
        RhTimeCourse = fft(RhTimeCourse);
        TransducinActivity = real(ifft(RhTimeCourse .* TrFilt));
    end

    % compress transducin activity if desired
    if (Condition.TransCompression > 0)
        TransducinActivity = TransducinActivity ./ (1 + Condition.TransCompression * cumsum(TransducinActivity));
    end
    if (Condition.TransRateCompression > 0)
        TransducinActivity = TransducinActivity ./ (1 + Condition.TransRateCompression * TransducinActivity);
    end
    
    ReturnedCondition.MeanTrTimeCourse = TransducinActivity;
            
	% convolve transducin time course with linear filter to generate modeled response
	if (Condition.ConvFlag)
        temp = conv(TransducinActivity, Filter, 'full');
        ReturnedCondition.EpochData = temp(1:Condition.EpochPts);
    else
        TransducinActivity = fft(TransducinActivity);
        ReturnedCondition.EpochData = real(ifft(TransducinActivity .* filt));
    end
    
else
    
    % generate series of responses
    for resp = 1:NumResponses

        RhTimeCourse = RhTrajectoryGenerator(Condition);
        TransducinActivity(1:Condition.EpochPts) = 0;
        CurrentStep = 0;
        CurrentCatalyticActivity = 1;
        TransCount = 0;

        % for each time point decide whether shutoff reaction has occurred.
        % if it has, update catalytic activity and shutoff rate.
        if (Condition.DeterministicTrModel == 0)
            if (Condition.DeterministicRhModel)
                TransGenerator = rand(1, Condition.EpochPts) ./ (MeanRhTimeCourse * Condition.BaseTransRate);           
            else
                TransGenerator = rand(1, Condition.EpochPts) ./ (RhTimeCourse * Condition.BaseTransRate);           
            end
            TransIndices = find(TransGenerator < 1);
            
            for cnt=1:length(TransIndices)

                TransducinLifetime = round(exprnd(Condition.TransDecayTimeConst));
                if ( (TransIndices(cnt) + TransducinLifetime) < Condition.EpochPts)
                    TransducinActivity(TransIndices(cnt):TransIndices(cnt)+TransducinLifetime) = TransducinActivity(TransIndices(cnt):TransIndices(cnt)+TransducinLifetime) + 1;
                else
                    TransducinActivity(TransIndices(cnt):Condition.EpochPts) = TransducinActivity(TransIndices(cnt):Condition.EpochPts) + 1;
                end
                TransCount = TransCount + 1;
            end
        end
        
        % generate transductin activity if not stochastic
        if (Condition.DeterministicTrModel)
            if (Condition.ConvFlag)
                temp = conv(TrFilt, MeanRhTimeCourse, 'full');
                TransducinActivity = temp(1:Condition.EpochPts);
            else
                RhTimeCourse = fft(RhTimeCourse);
                TransducinActivity = real(ifft(RhTimeCourse .* TrFilt));
            end
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

        % convolve transducin time course with linear filter to generate modeled response
        % convolve transducin time course with linear filter to generate modeled response
        if (Condition.ConvFlag)
            temp = conv(TransducinActivity, Filter, 'full');
            ReturnedCondition.EpochData = temp(1:Condition.EpochPts);
        else
            TransducinActivity = fft(TransducinActivity);
            ReturnedCondition.EpochData(resp, :) = real(ifft(TransducinActivity .* filt));
        end

        % compress if desired
        if (Condition.ResponseCompression > 0)
            ReturnedCondition.EpochData(resp, :) = (1-exp(-ReturnedCondition.EpochData(resp, :)*Condition.ResponseCompression))/Condition.ResponseCompression;
        end

        ReturnedCondition.TransCount(resp) = TransCount;
    end
    ReturnedCondition.MeanTrTimeCourse = ReturnedCondition.MeanTrTimeCourse / NumResponses;
end

% [Peakamp, PeakTime] = max(mean(ReturnedCondition.EpochData));
% NormSingles = resample(ReturnedCondition.EpochData', Condition.NumPtsToPeak, PeakTime) / Peakamp;
% ReturnedCondition.NormSingles = NormSingles;



