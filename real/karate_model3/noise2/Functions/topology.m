function topology()
load('karate.txt');
dlmwrite('Data/connectivity.dat', karate, 'delimiter', '\t', 'precision', 4);
end