numTethers = length(thr_struct);

attchPtMat = zeros(3,numTethers);
for ii = 1:numTethers
    if attach_point == 1
        attchPtMat(:,ii) = thr_struct(ii).R1_g;
    elseif attach_point == 2
        attchPtMat(:,ii) = thr_struct(ii).Rn_cm;
    end
    
end