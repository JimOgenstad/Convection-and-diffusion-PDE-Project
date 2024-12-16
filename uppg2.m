
M = 1000;
punkter = rand(M, 2);

v1AvPunkter = punkter(1:M, 2);
v2AvPunkter = 1 - punkter(1:M, 1);

fAvPunkter = fEnkel(punkter, M);

%Minsta kvadratanpassning. Ax = b -> A^TAx = A^Tb
N = 8;
nLista = -N:N;

minstaKvaA = zeros(M, (2*N+1)^2);


for k=1:M
    for n1 = nLista
        for n2 = nLista
            
            minstaKvaA(k, (n1 + N) * (2 * N + 1) + (n2 + N + 1)) = exp(1i*2*pi*(n1*punkter(k, 1)+n2*punkter(k, 2)));
        end
    end
end

ATASymetrisk = minstaKvaA'*minstaKvaA;
ATv1 = minstaKvaA'*v1AvPunkter;
ATv2 = minstaKvaA'*v2AvPunkter;
ATfunk = minstaKvaA'*fAvPunkter;

v1approx = ATASymetrisk\ATv1;
v2approx = ATASymetrisk\ATv2;
fapprox = ATASymetrisk\ATfunk;


%Med KxK stycken punkter på [0, 1]x[0, 1] vill vi nu vandra till randen med
%euler och beräkna f stegvis för att kunna bestämma denna i den
%ursprunkliga punkten




K = 50;
uSol = zeros(K, K);
delta = 0.01;
for i=1:K
    for j=1:K
        disp([i, j])
        x = i/K;
        y = j/K;
        uSum = 0;
        while x>0 && y>0
           
            xSuggest = x - delta*approxV1func(v1approx, x, y, nLista, N);
            ySuggest = y - delta*approxV2func(v2approx, x, y, nLista, N);
            
           
            if xSuggest>0 && ySuggest>0
                uSum = uSum + delta * approxF(fapprox, x, y, nLista, N);
            else
                
                break;
            end
            
            y = ySuggest;
            x = xSuggest;

        end   
            
        %Vi tar ett ministeg för att landa precis på randen. Vad ska
        %stegländgen vara?
        deltaVal1 = x/approxV1func(v1approx, x, y, nLista, N);
        deltaVal2 = y/approxV2func(v2approx, x, y, nLista, N);

        minDelta = min(deltaVal1, deltaVal2);


        uSum = uSum + minDelta*approxF(v1approx, x, y, nLista, N);
        uSol(j, i) = uSum;

        %Detta spelar mycket lite roll eftersom F kommer vara noll här
        %ändå.
        


    end
end



%plots
figure;
[X, Y] = meshgrid((1:K)/K, (1:K)/K); 
contour(X, Y, uSol, 20);
xlabel('X');
ylabel('Y');
title('Konturplot av uSol');
colorbar;


figure;
surf(X, Y, uSol, 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
zlabel('uSol');
title('Meshplot av uSol');
colorbar;







%%
function approxV1func = approxV1func(v1Bas, x, y, nList, N)
    summa = 0;
    for n1 = nList
        for n2 = nList
            summa = summa + v1Bas((n1 + N) * (2 * N + 1) + (n2 + N + 1))*exp(1i*2*pi*(n1*x+n2*y));
        end
    end

    approxV1func = real(summa);
end

function approxV2func = approxV2func(v2Bas, x, y, nList, N)
    summa = 0;
    for n1 = nList
        for n2 = nList
            summa = summa + v2Bas((n1 + N) * (2 * N + 1) + (n2 + N + 1))*exp(1i*2*pi*(n1*x+n2*y));
        end
    end

    approxV2func = real(summa);
end

function approxF = approxF(FBas, x, y, nList, N)
    summa = 0;
    for n1 = nList
        for n2 = nList
            summa = summa + FBas((n1 + N) * (2 * N + 1) + (n2 + N + 1))*exp(1i*2*pi*(n1*x+n2*y));
        end
    end

    approxF = real(summa);
end

function enkelProduktion = fEnkel(punkter, storlek)
    enkelProduktion = zeros(storlek, 1);

    for i=1:storlek
        if norm([0.5, 0.5] - punkter(i, 1:2)) < 0.1
            enkelProduktion(i) = 1;
        else
            enkelProduktion(i) = 0;
        end

    end
    
end