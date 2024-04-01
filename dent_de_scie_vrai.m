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

    [TabConcat, TabDDS, TabReconstr, reconstruction_DDS] = calculateSignals(H, a, f, t);

    plotResults(TabConcat, TabDDS, TabReconstr, reconstruction_DDS, Xmin, Xmax);
end

function [TabConcat, TabDDS, TabReconstr, reconstruction_DDS] = calculateSignals(H, a, f, t)
    TabConcat = zeros(2*H, 1);
    TabDDS = zeros(2*H, length(t));

    for i = 1:2*H
        CoeffDDS = calculateXDDS(a, i);
        TabConcat(i,1) = CoeffDDS;

        HarmoDDS = sin(2*pi*i*f*t);
        TabDDS(i,:) = CoeffDDS*HarmoDDS;
    end

    TabReconstr = cumsum(TabDDS);
    reconstruction_DDS = TabReconstr(2*H,:);
end

function CoeffDDS = calculateXDDS(a, i)
    CoeffDDS = ((-a)/(2*pi*i))*cos(i*pi);
end

function plotResults(TabConcat, TabDDS, TabReconstr, reconstruction_DDS, Xmin, Xmax)
    TabAbsolu = abs(TabConcat(find(abs(TabConcat)>=3.898171832519376e-17)));
    
    figure(1)
    hold on
    stem(abs(TabAbsolu))
    plot(abs(TabAbsolu),'r--')
    title('Les harmoniques du signal dent de scie')
    xlabel('Coefficent du signal dent de scie')
    ylabel('Amplitude')

    figure(2)
    plot(TabReconstr')
    title('Représentation 2D des harmoniques du signal dent de scie')
    xlabel('Temps(ms)')
    ylabel('Amplitude')
    axis([Xmin Xmax min(reconstruction_DDS)-0.1 max(reconstruction_DDS)+0.1])

    figure(3)
    plot(reconstruction_DDS)
    title('Signal dent de scie reconstruit grâce à la série de Fourier réelle ')
    xlabel('Temps (ms)')
    ylabel('Amplitude')
    axis([Xmin Xmax min(reconstruction_DDS)-0.1 max(reconstruction_DDS)+0.1])

    figure(4)
    mesh(TabReconstr)
    title('Visualisation 3D de la décomposition en série de Fourier réelle du signal dent de scie')
    xlabel('Temps (ms)')
    ylabel('Le rang des coefficients')
    zlabel('Amplitude')

    figure(5)
    stem(abs(fft(reconstruction_DDS)))
    axis([Xmin 60 0 160])
    title('Estimation de la TF du signal dent de scie grâce à l''algorithme FFT')
    xlabel('Coefficient du signal triangle')
    ylabel('Amplitude')

    figure(6)
    plot(TabDDS')
    title('Représentation des harmoniques pondérés')
    xlabel('Fréquence')
    ylabel('Amplitude')
    axis([Xmin Xmax+1 -0.7 0.7])

    figure(7)
    mesh(TabDDS)
    title('Rerésentation 3D de la décroissance en 1/n du signal dent de scie')
    xlabel('Le nombre d''harmoniques')
end
