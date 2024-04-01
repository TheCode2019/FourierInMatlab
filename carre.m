function main()
    clear all
    close all
    clc

    Xmax = 499;
    Xmin = 0;
    H = input('Combien d''harmoniques voulez-vous?');
    a = input('Choisissez votre amplitude');
    f = input('Entrez la frequence du signal');
    fe= input('Indiquez la frequence d''echantillonnage');
    fin=input('Donnez la fin de votre vecteur temps');
    t = 0:1/fe:fin;

    [TabConcat, TabCarre, TabReconstr, reconstruction_carre] = calculateSignals(H, a, f, t);

    plotResults(TabConcat, TabCarre, TabReconstr, reconstruction_carre, Xmin, Xmax);
end

function [TabConcat, TabCarre, TabReconstr, reconstruction_carre] = calculateSignals(H, a, f, t)
    TabConcat = zeros(2*H, 1);
    TabCarre = zeros(2*H, length(t));

    for i = 1:2*H
        CoeffCarre = calculateXCarre(a, i);
        TabConcat(i,1) = CoeffCarre;

        HarmoCarre = cos(2*pi*i*f*t);
        TabHarmos(i,:) = HarmoCarre;
        TabCarre(i,:) = CoeffCarre*HarmoCarre;
    end

    TabReconstr = cumsum(TabCarre);
    reconstruction_carre = TabReconstr(2*H,:);
end

function CoeffCarre = calculateXCarre(a, i)
    CoeffCarre = (2*a)/(pi*i)*sin(pi*i/2);
end

function plotResults(TabConcat, TabCarre, TabReconstr, reconstruction_carre, Xmin, Xmax)
    TabAbsolu = abs(TabConcat(find(abs(TabConcat)>=3.898171832519376e-17)));
    
    figure(1)
    hold on
    stem(abs(TabAbsolu))
    plot(abs(TabAbsolu),'r--')
    title('Les harmoniques du signal Carré')
    xlabel('Coefficent du signal carré')
    ylabel('Amplitude')

    figure(2)
    plot(TabReconstr')
    title('Représentation 2D des harmoniques du signal Carré')
    xlabel('Temps(ms)')
    ylabel('Amplitude')
    axis([Xmin Xmax min(reconstruction_carre)-0.1 max(reconstruction_carre)+0.1])

    figure(3)
    plot(reconstruction_carre)
    title('Signal carre reconstruit grâce à la série de Fourier réelle ')
    xlabel('Temps (ms)')
    ylabel('Amplitude')
    axis([Xmin Xmax min(reconstruction_carre)-0.1 max(reconstruction_carre)+0.1])

    figure(4)
    mesh(TabReconstr)
    title('Visualisation 3D de la décomposition en série de Fourier réelle du signal Carré')
    xlabel('Temps (ms)')
    ylabel('Le rang des coefficients')
    zlabel('Amplitude')

    figure(5)
    stem(abs(fft(reconstruction_carre)))
    axis([Xmin 60 0 160])
    title('Estimation de la TF du signal carre grâce à l''algorithme FFT')
    xlabel('Coefficient du signal carre')
    ylabel('Amplitude')

    figure(6)
    plot(TabCarre')
    title('Représentation des harmoniques pondérés')
    xlabel('Fréquence')
    ylabel('Amplitude')
    axis([Xmin Xmax+1 -0.7 0.7])

    figure(7)
    mesh(TabCarre)
    title('Rerésentation 3D de la décroissance en 1/n du signal Carré')
    xlabel('Le nombre d''harmoniques')
end
