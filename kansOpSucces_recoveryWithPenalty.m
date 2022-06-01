drempel = 10e-5;
rijen = 7;
kolommen = 6;
rang = 2;
aantalKeer = 10;
start= (7+6-2)*2;
aantalMetingen = 2*start;
y = zeros(aantalKeer*(aantalMetingen-start+1));
k=0;
m=0;
lambda = 10;
for i= start :aantalMetingen
    for j=1:aantalKeer
        k=k+1;
        y(k) = recoverMetPenalty(rijen,kolommen,rang,i,lambda);
    end
end

%%gemiddelde
g=0;
k=0;
succesRecovery = zeros((aantalMetingen-start+1),1);
for j=start:aantalMetingen
    g=g+1;
    succes =0;
    for i = 1:aantalKeer
        k=k+1;
        if (y(k)<=drempel)
            succes = succes + 1;
        end
    end
    succesRecovery(g) = succes/aantalKeer;
end
x = start : 1 : aantalMetingen;
fig1 = figure(1)
plot(x,succesRecovery,'--ks')
title('kans op succes (afwijking <= 1e-5) matrixrecovery 5x4-matrix van rang 2 met fminunc en gradient door Matlab en penalty op matrix ((2norm A²) - (2normB)²)² ')
xlabel(' metingen')
ylabel('succesvol hersteld [%]')
grid on