%% set parameters and initialize values
spot_interest_rate = 0.05,
mean_reversion = 0.1;
equilibrium = 0.05;
volatility = 0.01;
swap_rate = 0.05;


row_count = 40;

%allocate
time_interval = zeros(row_count,1);
B = zeros(row_count,1);
A = zeros(row_count,1);



IR = zeros(row_count+1,1);
fixed_leg = zeros(row_count+1,1);
floating_leg = zeros(row_count+1,1);
swap_mtm = zeros(row_count+1,1);


IR(1) = spot_interest_rate;

for i = 1 : row_count
    time_interval(i) = i*0.25;
    B(i) = (1-exp(-mean_reversion*time_interval(i)))/mean_reversion;
    A(i) = exp(((B(i)-time_interval(i))*(mean_reversion^2*equilibrium-...
           volatility^2/2))/mean_reversion^2-(volatility^2*B(i)^2/(4*...
           mean_reversion)));
end
%% calculation
for j = 1
        discount_factor = zeros((41-j),1);
        for k = 1 : (41-j)
            discount_factor(k) = A(k)*exp(-B(k)*IR(1));
        end
        fixed_leg(j) = swap_rate*sum(discount_factor)*0.25;
        floating_leg(j) = 1-discount_factor(41-j);
        swap_mtm(j) = fixed_leg(j)-floating_leg(j);
    end
        
    for j= 2: row_count
        ran = rand(row_count-1,1);
        IR(j)= IR(j-1)+mean_reversion*(equilibrium-IR(j-1))*0.25+...
               norminv(ran(j-1),0,1)*sqrt(0.25)*volatility;
        discount_factor = zeros((41-j),1);
        for k = 1 : (41-j)
            discount_factor(k) = A(k)*exp(-B(k)*IR(j));
        end
        fixed_leg(j) = swap_rate*sum(discount_factor)*0.25;
        floating_leg(j) = 1-discount_factor(41-j);
        swap_mtm(j) = fixed_leg(j)-floating_leg(j);
    end

    %finalize
    fixed_leg(row_count+1) = 0;
    floating_leg(row_count+1) = 0;
    swap_mtm(row_count+1) = fixed_leg(row_count+1)-...
                            floating_leg(row_count+1);
           
%% display output of monte carlo simulation
IR_monte_carlo_initialize;
time = linspace(0,10,41);

%simulation 250 times
simulation_output = zeros(row_count+1, 250);
for s = 1 : 250
    IR_monte_carlo_cal;
    simulation_output(:,s) = swap_mtm;
end


% calculate EE and PFE
EE = zeros(row_count+1,1);
PFE = zeros(row_count+1,1);

for m = 1 : row_count+1
    c = simulation_output(m,:);
    d = sort(c);
    EE(m) = sum(c(c>0))/250;
    PFE(m) = d(0.9*250);
end

figure
plot(time,EE,time,PFE);
xlabel('Time(years)')
ylabel('Exposure')
figure
plot(time,simulation_output);
xlabel('time(years)')
ylabel('future value')