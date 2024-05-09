function [cfar_targets] = cfar_ca_1D(data, numGuardCells, numRefCells, thresholdFactor)
    N = length(data);
    cfar_targets = zeros(N, 1);
    
    numTotalCells = numGuardCells + numRefCells;
    
    for i = numTotalCells + 1 : N - numTotalCells
        % Guard cells and target cell are excluded from the noise estimate
        leadingEdge = data(i - numTotalCells : i - numGuardCells - 1);
        trailingEdge = data(i + numGuardCells + 1 : i + numTotalCells);
        
        noiseLevel = mean([leadingEdge; trailingEdge]);
        threshold = noiseLevel * thresholdFactor;
        
        if data(i) > threshold
            cfar_targets(i) = 1;
        else
            cfar_targets(i) = 0;
        end
    end
    
end