function simulate3(MODEL,N,NI,S,M)

resolution=1;

models={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
if any(ismember(models,MODEL))==0
    disp('ERROR: MODEL must be a valid string: kuramoto1, kuramoto2, michaelis_menten, roessler');
else
    system('mkdir Data');
    disp('Creating network structure...')
    topology(N,'homogeneous','directed',NI);
    disp('Simulating time series...')
    Y=[];
    switch MODEL
        case 'kuramoto1'
            w=-2 + (4).*rand(N,1);
            dlmwrite('Data/frequencies.dat', w, 'delimiter', '\t', 'precision', 4);
            for s=1:S
                init=-3.14 +(3.14+3.14)*rand(N,1);
                tspan=0:resolution:M;
                [~,y] = ode45('kuramoto1',tspan,init);
                if rand<0.5
                    snr=10*rand;
                    px_dBW=10*rand;
                    y=awgn(y,snr,px_dBW);
%                     y=awgn(y,snr,'measured');
                end
                Y=[Y;y];
            end
        case 'kuramoto2'
            w=-2 + (4).*rand(N,1);
            dlmwrite('Data/frequencies.dat', w, 'delimiter', '\t', 'precision', 4);
            for s=1:S
                init=-3.14 +(3.14+3.14)*rand(N,1);
                tspan=0:resolution:M;
                [~,y] = ode45('kuramoto2',tspan,init);
                if rand<0.5
                    snr=10*rand;
                    px_dBW=10*rand;
                    y=awgn(y,snr,px_dBW);
%                     y=awgn(y,snr,'measured');
                end
                Y=[Y;y];
            end
        case 'michaelis_menten'
            for s=1:S
                init=1+1*rand(N,1);
                tspan=0:resolution:M;
                [~,y] = ode45('michaelis_menten',tspan,init);
                if rand<0.5
                    snr=10*rand;
                    px_dBW=10*rand;
                    y=awgn(y,snr,px_dBW);
%                     y=awgn(y,snr,'measured');
                end
                Y=[Y;y];
            end
        case 'roessler'
            for s=1:S
                init=-5 +(5+5)*rand(3*N,1);
                tspan=0:resolution:M;
                [~,y] = ode45('roessler',tspan,init);
                if rand<0.5
                    snr=10*rand;
                    px_dBW=10*rand;
                    y=awgn(y,snr,px_dBW);
%                     y=awgn(y,snr,'measured');
                end
                Y=[Y;y];
            end     
    end
    
    ts_param=[S,M];
    dlmwrite('Data/data.dat', Y, 'delimiter', '\t');
    dlmwrite('Data/ts_param.dat', ts_param, 'delimiter', '\t');
    clear;
    disp('Simulation finished!');
end
end