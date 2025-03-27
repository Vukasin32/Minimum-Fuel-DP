clear
close all
clc

%% KONSTANTE - MODEL AUTOMOBILA
M = 1365; Lambda = 1.089; Rw = 0.316 ; g = 9.81; 
id = 4.31; ig = [3.32, 2, 1.36, 1.01, 0.82]; % ig je konstanta koja se vezuje sa stepene prenosa tj. za 1,2,3,4 i 5 brzinu na menjaču 
nm = 0.95; cr = 0.015; S_cd = 0.65; ro_a = 1.2; 
phi = 0; phi = phi*pi/180; 

% Koeficijenti za potrošnju goriva
a1 = 1.1046e-5; a2 = -7.7511e-8; a3 = 1.6958e-10; 
a4 = 1.7363e-8; a5 = 6.4277e-11; a6 = 1.6088e-10;   
b1 = 1.5545; b2 = -4.8907e-3; b3 = 4.0442e-6;

% Koeficijenti za promenu stepena prenosa
t_shift = 0.3; sigma = 1.5; 

dist_total = 1000; % Ukupna distanca
v0 = 5*1000/3600; % Početna brzina (m/s)
vf = 67*1000/3600; % Krajnja brzina (m/s)
dv = 1*1000/3600; % Inkrement brzine
x = v0:dv:vf;

%% INICIJALIZACIJA MATRICA
n_states = length(x) * length(ig);
dist_covered = inf(n_states,n_states);
m_fuel = inf(n_states, n_states); % Potrošnja goriva za prelazak iz stanja u stanje
distance = inf(n_states, 1); % Pređena distanca
distance(1) = 0; % Početna distanca je 0

%% GENERISANJE STANJA
states = [];
for i = 1:length(x)
    for j = 1:length(ig)
        states = [states; x(i), j]; % Stanje se definiše kao 2D vektor čije su koordinate trenutna brzina automobila i trenutni stepen prenosa 
    end
end

%% POPUNJAVANJE MATRICE m_fuel
for i = 1:size(m_fuel,1)
    for j = i:size(m_fuel,2)
        if states(i,2) == states(j,2) && states(i,1) <= states(j,1) % Trenutna brzina ostaje ista ili raste dok stepen prenosa ostaje isti
            v_curr = states(i,1); v_new = states(j,1);
            ds = v_curr*sigma*t_shift; i_curr = states(i,2);
            v_mean = (v_curr+v_new)/2;

            dv_dt = (v_new - v_curr)*v_mean/ds;    
            F = Lambda * M * dv_dt + (S_cd * ro_a * v_mean^2) / 2 + ...
                cr * M * g * cos(phi) + M * g * sin(phi);
                
            omega = v_mean*id*ig(i_curr)/Rw;
            T = F*v_mean/(nm*omega);
                
            dm_dt = a1*omega + a2*omega^2 + a3*omega^3 + ...
                 a4*omega*T + a5*omega^2*T + a6*omega*T^2;
  
            dm = dm_dt*ds/v_mean;
            Tmax = b1*omega + b2*omega^2 + b3*omega^3;
                
            if T>Tmax % Ukoliko je momenat potreban za realizaciju prelaza veći od Tmax, smatra se da je prelaz nemoguć, tj. dm i ds se postavljaju na vrednost inf 
               dm = inf;
               ds = inf;
            end

        elseif (states(i,2) == states(j,2) - 1) && states(i,1) <= states(j,1) % Trenutna brzina ostaje ista ili raste dok stepen prenosa raste za 1
            v_curr = states(i,1); v_new = states(j,1);
            ds = v_curr*sigma*t_shift; i_curr = states(i,2); i_new = states(j,2);
            
            a_shift = ((S_cd * ro_a * v_curr^2) / 2 + ...
                        cr * M * g * cos(phi) + M * g * sin(phi))/(-Lambda*M);
            v_curr_prim = v_curr + a_shift*t_shift;
            
            v_mean = (v_curr + v_curr_prim)/2;
            ds_shift = v_mean*t_shift;
   
            v_mean = (v_new + v_curr_prim)/2;         
            dv_dt = (v_new - v_curr_prim)*v_mean/(ds-ds_shift);     
            F = Lambda * M * dv_dt + (S_cd * ro_a * v_mean^2) / 2 + ...
                cr * M * g * cos(phi) + M * g * sin(phi);
                
            omega = v_mean*id*ig(i_new)/Rw;
            T = F*v_mean/(nm*omega);
                
            dm_dt = a1*omega + a2*omega^2 + a3*omega^3 + ...
                 a4*omega*T + a5*omega^2*T + a6*omega*T^2;
  
            dm = dm_dt*(ds-ds_shift)/v_mean;
            Tmax = b1*omega + b2*omega^2 + b3*omega^3;
                
            if T>Tmax % Ukoliko je momenat potreban za realizaciju prelaza veći od Tmax, smatra se da je prelaz nemoguć, tj. dm i ds se postavljaju na vrednost inf  
               dm = inf;
               ds = inf;
            end

        else 
            dm = inf;
            ds = inf;
        end
        m_fuel(i,j) = dm; % Matrica m_fuel čuva potrebno gorivo za prelazak iz stanja i u stanje j 
        dist_covered(i,j) = ds; % Matrica dist_covered čuva distancu koju će auto preći pri prelasku iz stanja i u stanje j
    end
end

%% DINAMICKO PROGRAMIRANJE
n_states = size(m_fuel,1);
cost = inf*ones(n_states,1); 
prev_state = zeros(n_states, 1); % Čuvanje prethodnog stanja, potrebno da bi se rekonstruisala optimalna putanja

cost(1) = 0;
for i = 1:n_states
    for j = i:n_states
        if m_fuel(i,j) ~= inf
            new_cost = cost(i) + m_fuel(i,j);
            if new_cost < cost(j)
                cost(j) = new_cost;
                prev_state(j) = i;
            end
        end
    end
end

start_node = 0; end_node = 315; 
path = [end_node]; prev_node = prev_state(end_node); 

while prev_node ~= start_node
    path = [prev_node; path];
    prev_node = prev_state(prev_node);
end

min_cost = cost(end_node)*1000; % Minimalna potrošnja goriva da se dostigne stanje gde je trenutna brzina 67 km/h i stepen prenosa 5, izražena u gramima

disp('Minimalni trošak do finalnog stanja u gramima:');
disp(min_cost);

disp('Optimalan put:');
disp(path);

%% FINALNA POTROŠNJA I GRAFIČKI PRIKAZ REZULTATA
s_opt = 0; s_axis = [s_opt];
for i = 1:length(path)-1
    s_opt = s_opt + dist_covered(path(i),path(i+1)); % Ukupan pređeni put po optimalnoj putanji, od početnog do finalnog stanja
    s_axis = [s_axis; s_opt];
end

cost_remaining = (dist_total - s_opt)*m_fuel(315,315)/dist_covered(315,315)*1000; % Potrošnja goriva od finalnog stanja pa do trenutka kada se pređe ukupno 1000 m, izražena u gramima
cost_accumulated = (min_cost + cost_remaining);

disp('Minimalni potrošnja goriva za celokupan put u gramima:');
disp(cost_accumulated);

s_axis = [s_axis; 1000];
v_axis = states(path,1)*3600/1000; v_axis = [v_axis; v_axis(end)];
gear_axis = states(path,2); gear_axis = [gear_axis; gear_axis(end)];

figure()
plot(s_axis,v_axis)
hold all
scatter(s_axis,v_axis)
title('Speed switching')
grid("on")
xlabel('s[m]')
ylabel('v[km/h]')
ylim([0,70])

figure()
plot(s_axis,gear_axis)
hold all
scatter(s_axis,gear_axis)
title('Gear switching')
xlabel('s[m]')
ylabel('gear')
grid("on")
ylim([0,5.5])