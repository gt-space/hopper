altList = -50:5:0;
z = z_trajectory(1355:end);
pos1 = zeros(1,length(altList));
vel1 = zeros(1,length(altList));
for i = 1:length(altList)
    testVar = altList(i);
    [~, idx] = min(abs(z - testVar));
    closestValue = z(idx);
    pos1(i) = closestValue; 

    if closestValue < testVar
        idx1 = idx + 1354;
    else
        idx1 = idx + 1353;
    end
    interpVal = interp1([z_trajectory(idx1) z_trajectory(idx1+1)], [vz_trajectory(idx1), vz_trajectory(idx1+1)], testVar);
    %interpVal = interp1([z_trajectory(idx1) z_trajectory(idx1+1)], [vz_trajectory(idx1), vz_trajectory(idx+1)] testVar);
    vel1(i) = interpVal;

end
