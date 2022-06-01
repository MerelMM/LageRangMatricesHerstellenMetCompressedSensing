conditieV = 1e-9;
drempel = 10e-5;
rijen = 7;
kolommen = 6;
rang = 2;
oversamplingsNoemer = (kolommen+rijen-rang)*rang;
oversamplingsTellerStart = oversamplingsNoemer *1;
oversamplingsTellerStop = oversamplingsNoemer * 2; %42 metingen voor 5x4 rang 2
aantalKeer = 10;
y = zeros(aantalKeer*(oversamplingsTellerStop-oversamplingsTellerStart+1),1);
k=0;
m=0;
iteraties= zeros(aantalKeer*(oversamplingsTellerStop-oversamplingsTellerStart+1),1);

for i= oversamplingsTellerStart:oversamplingsTellerStop
    for j=1:aantalKeer
        k=k+1;
        [y(k), iteraties(k)] = gradientEnPenalty(rijen,kolommen,rang,i,10)
    end
end

%%gemiddelde
g=0;
k=0;
succesRecovery = zeros((oversamplingsTellerStop-oversamplingsTellerStart+1),1);
gemiddeldIter = zeros((oversamplingsTellerStop-oversamplingsTellerStart+1),1);
for j=oversamplingsTellerStart:oversamplingsTellerStop
    g=g+1;
    succes =0;
    gemiddeldI = 0;
    for i = 1:aantalKeer
        k=k+1;
        gemiddeldI = gemiddeldI + iteraties(k);
        if (y(k)<=drempel)
            succes = succes + 1;
        end
    end
    gemiddeldI
    succesRecovery(g) = succes/aantalKeer;
    gemiddeldIter(g) = gemiddeldI / aantalKeer
end

x= oversamplingsTellerStart/oversamplingsNoemer:(1/oversamplingsNoemer):oversamplingsTellerStop/oversamplingsNoemer;
fig2= figure(2)
plot(x,gemiddeldIter, '--ks')
title('gemiddelde aantal iteraties')
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]')
ylabel ('aantal iteraties')
grid on

fig1 = figure(1)
plot(x,succesRecovery,'--ks')
title('kans Op Succes met kolom A * 1')
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]')
ylabel('succesvol hersteld [%]')
grid on