%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidar_model_base.m
% 
% Function to output the distance measurement from each LIDAR given the
% true distance between the sensor and ground.
% Depends on filter_sim_initializer being ran first, which has bias/noise
% params.
% 
% Author: David Meisinger
% Date: 10/24/26
% 
% Inputs:
%   dist_true: true distances to measure (m) (4x1)
%   measurementNoise: optional additive Sensors measurement noise (26x1)
% 
% Outputs:
%   dist_meas: measured distances with noise (m) (4x1)
%   covariance: variance of each beam for the filter to use (4x1)
%   validity: one logical flag, true if at least one beam has a return
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist_meas, covariance, validity] = lidar_model_base(dist_true, measurementNoise) %#ok<INUSD>
%#genned but reviewed
% One true range or four beam ranges (m). Distance and covariance are 4-by-1.
% The histogram supplies the noise. The optional second input is ignored.

    useBrick = true;  % false selects the white wall
    %checks to make sure that imput is valid
    assert(isnumeric(dist_true) && isreal(dist_true), 'Ranges must be real numbers.');

    %rescales to deal with each individual lidar measurement seperately
    if isscalar(dist_true)
        d = repmat(double(dist_true),4,1);
    else
        assert(isvector(dist_true) && numel(dist_true)==4, 'Pass one range or four beam ranges.');
        d = reshape(double(dist_true),4,1);
    end
    %gets the data set we need
    persistent H
    if isempty(H), H = fixedHistograms(); end
    if useBrick
        runs = H.brick;
    else
        runs = H.wall;
    end

    %How far away the recorded measurement is
    targets = [runs.nominal_m];
    %check if the true distance is within parameters Add a 1 meter grace
    %period on top of the last measured distance. Thus only goes up to 4, 5
    %for the brick and wall measurements respectively.
    beamValid = isfinite(d) & d >= 0.2 & d <= targets(end) + 1;
    validity = any(beamValid);
    % A zero configured noise level requests the true range for that distance band.
    params = nav_sensor_params();
    noiseOff = (d <= 3 & params.lidar.noise_std_1 == 0) | ...
               (d > 3 & params.lidar.noise_std_2 == 0);
    useHistogram = beamValid & ~noiseOff;
    
    [~, runIndex] = min(abs(d - targets), [], 2);
    recordedStd = [runs.std_dev_m].';
    varying = useHistogram & recordedStd(runIndex) > 0;
    draw = zeros(4,1);
    draw(varying) = rand(nnz(varying),1);

    dist_meas = zeros(4,1);
    dist_meas(beamValid) = d(beamValid);
    covariance = zeros(4,1);
    
    for k = 1:4
        if useHistogram(k)
            h = runs(runIndex(k));
            sample_m = h.mean_m;  % The recorded value when the spread is zero.
            if varying(k)
                bin = 1 + sum(h.cdf <= draw(k));
                sample_m = h.reading_m(bin);
                covariance(k) = h.std_dev_m^2;
            end
            %adds error
            dist_meas(k) = d(k) + (sample_m - h.nominal_m);
            %computes covarience which is a vector containing the standard
            %deviation of the chosen distribution squared. Will inflate
            %to account for bias. Added additional
            %checks in case more tuning is needed. This is only for the
            %brick measurements a better solution should be developed
            %eventually. uses the squared measured bias as the inflation as
            %opposed to tunable parameters for consistency.
            Rinflation = (h.mean_m - h.nominal_m)^2;  % m²
            covariance(k) = h.std_dev_m^2 + Rinflation;
        end
    end
end
%% This is just a datastructure that stores histograms
function H = fixedHistograms()
    % Each entry gives the filename, nominal range (m), histogram limits
    % (mm), and [bucket number, recorded range in mm, observation count].
    % The filenames are labels. All data are embedded here.
    H.brick = [ ...
        makeHist("qtr_meter_brick_run.csv",.25,[220 230],[1 220 5022; 20 230 6146]), ...
        makeHist("half_meter_brick_run.csv",.5,[470 480],[1 470 6560; 20 480 5007]), ...
        makeHist("1m_brick_run.csv",1,[970 980],[1 970 1282; 20 980 10985]), ...
        makeHist("2m_brick_run.csv",2,[1960 2000], ...
            [1 1960 4; 6 1970 854; 11 1980 7427; 16 1990 3718; 20 2000 82]), ...
        makeHist("3m_brick_run.csv",3,[2940 3020], ...
            [1 2940 6; 3 2950 125; 6 2960 1118; 8 2970 3734; 11 2980 4654; ...
             13 2990 1666; 16 3000 209; 18 3010 21; 20 3020 2])];

    H.wall = [ ...
        makeHist("qtr_meter_run.csv",.25,[210 220],[1 210 18; 20 220 11298]), ...
        makeHist("half_meter_run.csv",.5,[490 510],[1 490 12317]), ...
        makeHist("1m_run.csv",1,[990 1000],[1 990 4773; 20 1000 6348]), ...
        makeHist("2m_run.csv",2,[1980 1990],[1 1980 4850; 20 1990 6416]), ...
        makeHist("3m_run.csv",3,[2990 3010],[1 2990 237; 11 3000 10171; 20 3010 607]), ...
        makeHist("4m_run.csv",4,[3990 4030], ...
            [1 3990 47; 6 4000 1696; 11 4010 7575; 16 4020 1835; 20 4030 13])];
end

function h = makeHist(file, nominal_m, limits_mm, occupied)
    counts = zeros(20,1);
    readings = zeros(20,1);
    indices = occupied(:,1);
    counts(indices) = occupied(:,3);
    readings(indices) = occupied(:,2)/1000;

    n = sum(counts);
    p = counts/n;
    mu = sum(p .* readings);
    edges = linspace(limits_mm(1),limits_mm(2),21)'/1000;

    h = struct('file',file,'nominal_m',nominal_m,'n',n,'num_buckets',20, ...
        'buckets_m',[edges(1:end-1),edges(2:end)],'probabilities',p, ...
        'reading_m',readings,'cdf',cumsum(counts)/n,'mean_m',mu, ...
        'std_dev_m',sqrt(sum(p .* (readings-mu).^2)));
end
