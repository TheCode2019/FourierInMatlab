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

    [TabConcat, TabTriangle, TabReconstr, reconstruction_triangle] = calculateSignals(H, a, f, t);

    plotResults(TabConcat, TabTriangle, TabReconstr, reconstruction_triangle, Xmin, Xmax);
end

function [TabConcat, TabTriangle, TabReconstr, reconstruction_triangle] = calculateSignals(H, a, f, t)
    TabConcat = zeros(2*H, 1);
    TabTriangle = zeros(2*H, length(t));

    for i = 1:2*H
        CoeffTriangle = calculateXTriangle(a, i);
        TabConcat(i,1) = CoeffTriangle;

        HarmoTriangle = sin(2*pi*i*f*t);
        TabHarmos(i,:) = HarmoTriangle;
        TabTriangle(i,:) = CoeffTriangle*HarmoTriangle;
    end

    TabReconstr = cumsum(TabTriangle);
    reconstruction_triangle = TabReconstr(2*H,:);
end

function CoeffTriangle = calculateXTriangle(a, i)
    CoeffTriangle = (-a)/(((i))^2)*sin((i*pi)/2);
end

function plotResults(TabConcat, TabTriangle, TabReconstr, reconstruction_triangle, Xmin, Xmax)
    TabAbsolu = abs(TabConcat(find(abs(TabConcat)>=3.898171832519376e-17)));
    
    figure(1)
    hold on
    stem(abs(TabAbsolu))
    plot(abs(TabAbsolu),'r--')
    title('Les harmoniques du signal triangle')
    xlabel('Coefficent du signal triangle')
    ylabel('Amplitude')

    figure(2)
    plot(TabReconstr')
    title('Représentation 2D des harmoniques du signal triangle')
    xlabel('Temps(ms)')
    ylabel('Amplitude')
    axis([Xmin Xmax min(reconstruction_triangle)-0.1 max(reconstruction_triangle)+0.1])

    figure(3)
    plot(reconstruction_triangle)
    title('Signal triangle reconstruit grâce à la série de Fourier réelle ')
    xlabel('Temps (ms)')
    ylabel('Amplitude')
    axis([Xmin Xmax min(reconstruction_triangle)-0.1 max(reconstruction_triangle)+0.1])

    figure(4)
    mesh(TabReconstr)
    title('Visualisation 3D de la décomposition en série de Fourier réelle du signal triangle')
    xlabel('Temps (ms)')
    ylabel('Le rang des coefficients')
    zlabel('Amplitude')

    figure(5)
    stem(abs(fft(reconstruction_triangle)))
    axis([Xmin 60 0 160])
    title('Estimation de la TF du signal triangle grâce à l''algorithme FFT')
    xlabel('Coefficient du signal triangle')
    ylabel('Amplitude')

    figure(6)
    plot(TabTriangle')
    title('Représentation des harmoniques pondérés')
    xlabel('Fréquence')
    ylabel('Amplitude')
    axis([Xmin Xmax+1 -0.7 0.7])

    figure(7)
    mesh(TabTriangle)
    title('Rerésentation 3D de la décroissance en 1/n² du signal triangle')
    xlabel('Le nombre d''harmoniques')
end
