drempel = 10e-5;
rijen = 7;
kolommen = 6;
rang = 2;
start = (rijen+kolommen-rang)*rang;
aantalMetingen = 2*start;
aantalKeer = 20;
y = zeros(aantalKeer*(aantalMetingen-start+1),1);
k=0;
m=0;
iter = zeros(aantalKeer*(aantalMetingen-start+1),1);
slechteAB = zeros(aantalKeer*(aantalMetingen-start+1),1);
slechteAB5 = zeros(aantalKeer*(aantalMetingen-start+1),1);
slechteAB10 = zeros(aantalKeer*(aantalMetingen-start+1),1);

for i= start:aantalMetingen
    for j=1:aantalKeer
        k=k+1;
        [y(k),iter(k),slechteAB(k),slechteAB5(k),slechteAB10(k)] = recover(rijen,kolommen,rang,i);
    end
end

%%gemiddelde
g=0;
k=0;
succesRecovery = zeros((aantalMetingen-start+1),1);
plotIteraties = zeros((aantalMetingen-start+1),1);
plotSlechteAB = zeros((aantalMetingen-start+1),1);
plotSlechteAB5 = zeros((aantalMetingen-start+1),1);
plotSlechteAB10 = zeros((aantalMetingen-start+1),1);
for j=start:aantalMetingen
    g=g+1;
    succes =0;
    gemiddeldIteraties = 0;
    succesSlechteAB=0;
    succesSlechteAB5=0;
    succesSlechteAB10=0;

    for i = 1:aantalKeer
        k=k+1;
        gemiddeldIteraties = gemiddeldIteraties + iter(k);
        if (y(k)<=drempel)
            succes = succes + 1;
        end
        if (slechteAB(k)<=drempel)
            succesSlechteAB = succesSlechteAB + 1;
        end
        if (slechteAB5(k)<=drempel)
            succesSlechteAB5 = succesSlechteAB5 + 1;
        end
        if (slechteAB10(k)<=drempel)
            succesSlechteAB10 = succesSlechteAB10 + 1;
        end
    end
    
    succesRecovery(g) = succes/aantalKeer;
    plotIteraties(g) = gemiddeldIteraties/aantalKeer;
    plotSlechteAB(g) = succesSlechteAB/aantalKeer;
    plotSlechteAB5(g) = succesSlechteAB5/aantalKeer;
    plotSlechteAB10(g) = succesSlechteAB10/aantalKeer;
end
x=start/start:1/start:aantalMetingen/start;
fig1 = figure(1)
plot(x,succesRecovery,'--ks')
title('\alpha = 1',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig2 = figure(2)
plot(x,plotSlechteAB,'--ks')
title('\alpha = 1e1',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig3 = figure(3)
plot(x,plotSlechteAB5,'--ks')
title('\alpha = 1e2',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

fig4 = figure(4)
plot(x,plotSlechteAB10,'--ks')
title('\alpha = 1e3',fontsize=16)
xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
ylabel('succesvol hersteld [%]',fontsize=16)
grid on

%fig5 = figure(5)
%plot(x,plotIteraties,'--ks')
%title('gemiddeld aantal iteraties',fontsize=16)
%xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
%ylabel('succesvol hersteld [%]',fontsize=16)
%grid on

