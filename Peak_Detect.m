function DetectedPeaks = Peak_Detect(SigToInspect, StartIndex, FinishIndex)
% Peak_Detect(Signal to Inspect, Starting Index, Finishing Index)
% Function inspects the signal for peaks & valleys between indexes provided
% and stores the sample numbers in a vector which is passed back.

DetectedPeaks = [];

for N = StartIndex : FinishIndex
    % if the value at the current index is less then the last index save
    % the previous index to peaks
    if ( abs(SigToInspect(N-1)) < abs(SigToInspect(N)) ) && ( abs(SigToInspect(N+1)) < abs(SigToInspect(N)) ); 
        DetectedPeaks = [DetectedPeaks N];      
    end            
        
end


end