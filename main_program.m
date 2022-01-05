%% Clear Workspace
close all
clear all
clc
%% Add Size to Box
inputCells = inputdlg({'Frame X Size', 'Frame Y Size', 'Frame Z Size'});
frame2X = str2double(inputCells{1});
frame2Y = str2double(inputCells{2});
frame2Z = str2double(inputCells{3});
inputCells = inputdlg({'Box X Size', 'Box Y Size', 'Box Z Size'});
box2X = str2double(inputCells{1});
box2Y = str2double(inputCells{2});
box2Z = str2double(inputCells{3});
inputCells = inputdlg('Input points');
nPoints = str2double(inputCells{1});
%
xLimits = [-frame2X/2 frame2X/2];
yLimits = [-frame2Y/2 frame2Y/2];
zLimits = [-frame2Z/2 frame2Z/2];
% open new blank figure with defined limits
figure
points = zeros(nPoints + 1, 3);
for i=1:nPoints
    %view([0,0,1]);
    xlabel('X');
    ylabel('Y');
    axis('equal');
    title('Set XY');
    xlim(xLimits);
    ylim(yLimits);
    coords = ginput(1);
    points(i, [1, 2]) = coords;
    if i == 1
        points(nPoints + 1, [1, 2]) = coords;
    end
    plot(points(:, 1), points(:, 2));
end 

for i=1:nPoints
    plot(points(:, 1), points(:, 3));
    hold on
    plot(points(i, 1), points(i, 3), 'r.', 'MarkerSize', 20);
    hold off
    %view([0,1,0]);
    xlabel('X');
    ylabel('Z');
    axis('equal');
    title('Set Z');
    xlim(xLimits);
    ylim(zLimits);
    coords = ginput(1);
    points(i, 3) = coords(2);
    if i == 1
        points(nPoints + 1, 3) = coords(2);
    end
end
plot3(points(:, 1), points(:, 2), points(:, 3));
spline_result = cscvn(points');
title('Curve Smoothing Result');
hold on
fnplt(spline_result, 'r', 2);
hold off

spline_start = spline_result.breaks(1);
spline_end = spline_result.breaks(end);
%% Fourier Reconstruction
linspacePts = 8192;
val_domain = linspace(spline_start, spline_end, linspacePts + 1);
val_domain = val_domain(1:end-1);
centerPosition = fnval(spline_result, val_domain);
trajX = centerPosition(1, :);
trajY = centerPosition(2, :);
trajZ = centerPosition(3, :);
% Fourier Reconstruction
[reconPX, reconVX, reconAX] = fourierRecon(trajX, linspacePts, nPoints);
[reconPY, reconVY, reconAY] = fourierRecon(trajY, linspacePts, nPoints);
[reconPZ, reconVZ, reconAZ] = fourierRecon(trajZ, linspacePts, nPoints);
centerPosition = [reconPX; reconPY; reconPZ];
figure
subplot(311);
hold on
plot(reconPX);
plot(reconPY);
plot(reconPZ);
hold off
title('Center Position');
legend('X', 'Y', 'Z');
xlabel('Frame');

subplot(312);
hold on
plot(reconVX);
plot(reconVY);
plot(reconVZ);
hold off
title('Center Velocity');
legend('X', 'Y', 'Z');
xlabel('Frame');

subplot(313);
hold on
plot(reconAX);
plot(reconAY);
plot(reconAZ);
hold off
title('Center Acceleration');
legend('X', 'Y', 'Z');
xlabel('Frame');

figure
ax1 = subplot(131);
plot3(reconPX, reconPY, reconPZ);
title('3D Position');
xlabel('X');
ylabel('Y');
zlabel('Z');

ax2 = subplot(132);
plot3(reconVX, reconVY, reconVZ);
title('3D Velocity');
xlabel('X');
ylabel('Y');
zlabel('Z');

ax3 = subplot(133);
plot3(reconAX, reconAY, reconAZ);
title('3D Acceleration');
xlabel('X');
ylabel('Y');
zlabel('Z');

Link = linkprop([ax1, ax2, ax3],{'View'});
setappdata(gcf, 'StoreTheLink', Link);
%% Cable Lengths
innerCorners = cell([1, 8]);
outerCorners = cell([1, 8]);
directionVectors = cell([1, 8]);
cableLengths = cell([1, 8]);
legendsC = cell([1, 8]);
boxCoords = cell([1, 8]);
figure
hold on
for i=1:8
    positionM = double(dec2bin(i - 1, 3) - '0') * 2 - 1; % -1, -1, -1 thru 1, 1, 1
    boxCoords{i} = (positionM .* [box2X/2 box2Y/2 box2Z/2])';
    innerCorners{i} = centerPosition + boxCoords{i};
    outerCorners{i} = (positionM .* [frame2X/2 frame2Y/2 frame2Z/2])';
    directionVectors{i} = outerCorners{i} - innerCorners{i};
    cableLengths{i} = vecnorm(directionVectors{i});
    plot(cableLengths{i});
    legendsC{i} = ['Cable ' num2str(i)];
end
hold off
legend(legendsC);
title('Cable Lengths');
xlabel('Frame');
%% Add External Force
inputCells = inputdlg({'External Force X', 'External Force Y', 'External Force Z', 'External Torque X', 'External Torque Y', 'External Torque Z'});
ExtForceTorque = str2double(inputCells)';
forceScrew = cell([1, 8]);
cableForces = zeros(8, linspacePts);
residues = zeros(6, linspacePts);
screwMatrices = cell([1, linspacePts]);
%figure
%hold on
for j=1:linspacePts
    screwMatrix = zeros(6, 8);
    for i=1:8
        normVector = directionVectors{i}(:, j) / norm(directionVectors{i}(:, j));
        forceScrewRF = cross(boxCoords{i}, normVector); % use floating "center" of the box as reference
        forceScrewF = [normVector ; forceScrewRF]';
        screwMatrix(:, i) = forceScrewF;
    end
    screwMatrices{j} = screwMatrix;
    %cableForce = lsqnonneg(screwMatrix,ExtForceTorque'); %;
    cableForce = pinv(screwMatrix)*ExtForceTorque';% + (eye(8) - (pinv(screwMatrix)*screwMatrix)) *([1;1;1;1;1;1;1;1]*2*sum(abs(ExtForceTorque))); %;
    %cableForce = screwMatrix \ ExtForceTorque'; %;
    %residues(:, j) = (screwMatrix * cableForce - ExtForceTorque')';
    cableForces(:, j) = cableForce';%(cableForce + null(screwMatrix,'r') * SolK)';
end

correctionTerm = -min(cableForces, 2);

for j=1:linspacePts
    screwMatrix = screwMatrices{j};
    cableForce = pinv(screwMatrix)*ExtForceTorque' + (eye(8) - (pinv(screwMatrix)*screwMatrix)) *(correctionTerm*sum(abs(ExtForceTorque))); %;
    residues(:, j) = (screwMatrix * cableForce - ExtForceTorque')';
    cableForces(:, j) = cableForce';
end
figure
plot(cableForces');
%hold off
title('Cable Forces');
legend(legendsC);
xlabel('Frame');
figure
plot(residues')
legend('F_X', 'F_Y', 'F_Z', 'T_X', 'T_Y', 'T_Z');
title('Cable Forces Error [ScrewMatrix][CableForces] - [ExternalScrewLoad]');
xlabel('Frame');
residueRMSE = sqrt(mean(residues' .^ 2, 'all'))
%% Reconstruction using fourier basis
function [reconP, reconV, reconA] = fourierRecon(traj, linspacePts, nPoints)

fftP = fft(traj - mean(traj));
magP = abs(fftP(1:linspacePts/2));
angleP = angle(fftP(1:linspacePts/2));
[~, inds] = maxk(magP, nPoints + 1);
basePi = linspace(0, 2*pi, linspacePts+1);
basePi = basePi(1:end-1);
freqs = (inds - 1)' * basePi;
reconP = sum(magP(inds) * cos(freqs' + angleP(inds))', 1) .* (2 / linspacePts) + mean(traj);
reconV = sum((magP(inds).*(inds-1)) * cos(freqs' + angleP(inds) + pi/2)', 1) .* (4*pi / (linspacePts ^ 2));
reconA = sum((magP(inds).*((inds-1).^2)) * cos(freqs' + angleP(inds) + pi)', 1) .* (8*pi / (linspacePts ^ 3));

end
