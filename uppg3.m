
xStart = 0;
xEnd = pi;
xPart = 1000;

tStart = 0;
tEnd = 4;
tPart = 1000;

n = 10;

xVec = linspace(xStart, xEnd, xPart);
tVec = linspace(tStart, tEnd, tPart);

fourierEnkelSol = fourierEnkel(xVec, tVec', n);
[X, T] = meshgrid(xVec, tVec);
figure;
surf(X, T, fourierEnkelSol, 'EdgeColor', 'none');
xlabel('X');
ylabel('T');
zlabel('u');
title('FourierEnkel');
colorbar;

n = 100;
fourierDisSol = fourierDiskont(xVec, tVec', n);
figure;
surf(X, T, fourierDisSol, 'EdgeColor', 'none');
xlabel('X');
ylabel('T');
zlabel('u');
title('FourierDiskontinuerlig');
colorbar;


deltaX = pi/10; 
deltaT = 1/10^2;
lambda = deltaT/deltaX^2;

diskT = 0:deltaT:4;
diskX = 0:deltaX:pi;

finitenkelSol = zeros(length(diskX), length(diskT));
enkelStartVektor = diskX.*(pi-diskX);
finitenkelSol(1:length(diskX), 1) = enkelStartVektor;

for t = diskT(2:end)
    for x = diskX(2:end-1)
        addOn = lambda*(finitenkelSol(round(x/deltaX+2), round(t/deltaT)) - 2*finitenkelSol(round(x/deltaX+1), round(t/deltaT)) + finitenkelSol(round(x/deltaX), round(t/deltaT)));
        finitenkelSol(round(x/deltaX+1), round(t/deltaT+1)) = finitenkelSol(round(x/deltaX+1), round(t/deltaT)) + addOn;
    end
end




figure;
surf(diskX, diskT, finitenkelSol', 'EdgeColor', 'none');
xlabel('X');
ylabel('T');
zlabel('u');
title('FinitEnkel');
colorbar;


finitDiskontSol = zeros(length(diskX), length(diskT));
diskontStartVektor = diskX .* (diskX < pi/2);
finitDiskontSol(1:length(diskX), 1) = diskontStartVektor;

for t = diskT(2:end)
    for x = diskX(2:end-1)
        addOn = lambda*(finitDiskontSol(round(x/deltaX+2), round(t/deltaT)) - 2*finitDiskontSol(round(x/deltaX+1), round(t/deltaT)) + finitDiskontSol(round(x/deltaX), round(t/deltaT)));
        finitDiskontSol(round(x/deltaX+1), round(t/deltaT+1)) = finitDiskontSol(round(x/deltaX+1), round(t/deltaT)) + addOn;
    end
end

figure;
surf(diskX, diskT, finitDiskontSol', 'EdgeColor', 'none');
xlabel('X');
ylabel('T');
zlabel('u');
title('FinitDiskont');
colorbar;



%% Fourier Konvergens

xStart = 0;
xEnd = pi;
xPart = 1000;

tStart = 0;
tEnd = 4;
tPart = 1000;

Nlist = [10, 20, 50, 100, 200, 500, 2000];

xVec = linspace(xStart, xEnd, xPart);
tVec = linspace(tStart, tEnd, tPart);
A = zeros(xPart, tPart, length(Nlist));

for i = 1:length(Nlist)
    A(:,:,i) = fourierEnkel(xVec, tVec', Nlist(i));
end

errors = zeros(length(Nlist)-1);
for i = 1:length(errors)
    errors(i) = norm(A(:,:,i)-A(:,:,length(Nlist)));
end

figure;
loglog(Nlist(1:end-1), errors);
hold on
x = 1:Nlist(end-1);
loglog(x, 1./x.^(2.5));


B = zeros(xPart, tPart, length(Nlist));

for i = 1:length(Nlist)
    B(:,:,i) = fourierDiskont(xVec, tVec', Nlist(i));
end

errors = zeros(length(Nlist)-1);
for i = 1:length(errors)
    errors(i) = norm(B(:,:,i)-B(:,:,length(Nlist)));
end

figure;
loglog(Nlist(1:end-1), errors);
hold on
loglog(x, 1./sqrt(x));

%% Finit differens konvergens

xParts = [10, 30, 50, 100, 200];
deltaXs = pi./xParts;
lambda = 1/pi^2;
deltaTs = lambda.*deltaXs.^2;

diskT = 0:deltaTs(1):4;
diskX = 0:deltaXs(1):pi;
finitbasSol = zeros(length(diskX), length(diskT));
enkelStartVektor = diskX .* (diskX < pi/2);
finitbasSol(1:length(diskX), 1) = enkelStartVektor;
deltaX = deltaXs(1);
deltaT = deltaTs(1);

for t = diskT(2:end)
    for x = diskX(2:end-1)
        addOn = lambda*(finitbasSol(round(x/deltaX+2), round(t/deltaT)) - 2*finitbasSol(round(x/deltaX+1), round(t/deltaT)) + finitbasSol(round(x/deltaX), round(t/deltaT)));
        finitbasSol(round(x/deltaX+1), round(t/deltaT+1)) = finitbasSol(round(x/deltaX+1), round(t/deltaT)) + addOn;
    end
end



diskT = 0:deltaTs(2):4;
diskX = 0:deltaXs(2):pi;
finit2Sol = zeros(length(diskX), length(diskT));
enkelStartVektor = diskX .* (diskX < pi/2);
finit2Sol(1:length(diskX), 1) = enkelStartVektor;
deltaX = deltaXs(2);
deltaT = deltaTs(2);

for t = diskT(2:end)
    for x = diskX(2:end-1)
        addOn = lambda*(finit2Sol(round(x/deltaX+2), round(t/deltaT)) - 2*finit2Sol(round(x/deltaX+1), round(t/deltaT)) + finit2Sol(round(x/deltaX), round(t/deltaT)));
        finit2Sol(round(x/deltaX+1), round(t/deltaT+1)) = finit2Sol(round(x/deltaX+1), round(t/deltaT)) + addOn;
    end
end



diskT = 0:deltaTs(3):4;
diskX = 0:deltaXs(3):pi;
finit3Sol = zeros(length(diskX), length(diskT));
enkelStartVektor = diskX .* (diskX < pi/2);
finit3Sol(1:length(diskX), 1) = enkelStartVektor;
deltaX = deltaXs(3);
deltaT = deltaTs(3);

for t = diskT(2:end)
    for x = diskX(2:end-1)
        addOn = lambda*(finit3Sol(round(x/deltaX+2), round(t/deltaT)) - 2*finit3Sol(round(x/deltaX+1), round(t/deltaT)) + finit3Sol(round(x/deltaX), round(t/deltaT)));
        finit3Sol(round(x/deltaX+1), round(t/deltaT+1)) = finit3Sol(round(x/deltaX+1), round(t/deltaT)) + addOn;
    end
end




diskT = 0:deltaTs(4):4;
diskX = 0:deltaXs(4):pi;
finit4Sol = zeros(length(diskX), length(diskT));
enkelStartVektor = diskX .* (diskX < pi/2);
finit4Sol(1:length(diskX), 1) = enkelStartVektor;
deltaX = deltaXs(4);
deltaT = deltaTs(4);

for t = diskT(2:end)
    for x = diskX(2:end-1)
        addOn = lambda*(finit4Sol(round(x/deltaX+2), round(t/deltaT)) - 2*finit4Sol(round(x/deltaX+1), round(t/deltaT)) + finit4Sol(round(x/deltaX), round(t/deltaT)));
        finit4Sol(round(x/deltaX+1), round(t/deltaT+1)) = finit4Sol(round(x/deltaX+1), round(t/deltaT)) + addOn;
    end
end



diskT = 0:deltaTs(5):4;
diskX = 0:deltaXs(5):pi;
finitFinalSol = zeros(length(diskX), length(diskT));
enkelStartVektor = diskX .* (diskX < pi/2);
finitFinalSol(1:length(diskX), 1) = enkelStartVektor;
deltaX = deltaXs(5);
deltaT = deltaTs(5);

for t = diskT(2:end)
    for x = diskX(2:end-1)
        addOn = lambda*(finitFinalSol(round(x/deltaX+2), round(t/deltaT)) - 2*finitFinalSol(round(x/deltaX+1), round(t/deltaT)) + finitFinalSol(round(x/deltaX), round(t/deltaT)));
        finitFinalSol(round(x/deltaX+1), round(t/deltaT+1)) = finitFinalSol(round(x/deltaX+1), round(t/deltaT)) + addOn;
    end
end






[xbas, ybas] = meshgrid(1:size(finitbasSol, 2), 1:size(finitbasSol, 1));
[x2, y2] = meshgrid(1:size(finit2Sol, 2), 1:size(finit2Sol, 1));
[x3, y3] = meshgrid(1:size(finit3Sol, 2), 1:size(finit3Sol, 1));
[x4, y4] = meshgrid(1:size(finit4Sol, 2), 1:size(finit4Sol, 1));
[xFinal, yFinal] = meshgrid(1:size(finitFinalSol, 2), 1:size(finitFinalSol, 1));



finit2interp = interp2(x2, y2, finit2Sol, xbas, ybas, 'linear');
finit3interp = interp2(x3, y3, finit3Sol, xbas, ybas, 'linear');
finit4interp = interp2(x4, y4, finit4Sol, xbas, ybas, 'linear');
finitFinalinterp = interp2(xFinal, yFinal, finitFinalSol, xbas, ybas, 'linear');

error = zeros(4);
error(1) = abs(finitFinalinterp(6, 101)-finitbasSol(6, 101));
error(2) = abs(finitFinalinterp(6, 101)-finit2interp(6, 101));
error(3) = abs(finitFinalinterp(6, 101)-finit3interp(6, 101));
error(4) = abs(finitFinalinterp(6, 101)-finit4interp(6, 101));

loglog(deltaTs(1:4), error);
hold on
xcurve = 0.0001:0.0001:0.01;
loglog(xcurve, sqrt(xcurve));



%%
function u = fourierEnkel(x, t, n)
    sum = 0;
    for i=1:n
        sum = sum + (8 ./ (pi*(2*i-1)^3)) .* exp(-t .* (2*i-1)^2) .* sin(x .* (2*i-1));
    end
    u = sum;
end

function u = fourierDiskont(x, t, n)
    sum = 0;
    for i =1:n
        a_1 = 2*sin(i*pi/2)/(pi*i^2);
        a_2 = -cos(i*pi/2)/i;
        sum = sum + (a_1+a_2).* exp(-t .* i^2) .* sin(i * x);
    end
    u = sum;


end

