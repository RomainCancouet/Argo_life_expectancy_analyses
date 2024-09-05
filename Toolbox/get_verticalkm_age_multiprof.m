function [vertical_km, vertical_km_mean, floats_age, last_cycle] = get_verticalkm_age_multiprof(Floats, dac_dir)
% EXAMPLE: [vertical_km, vertical_km_mean, floats_age, last_cycle] = get_verticalkm_age_multiprof(Floats, dac_dir)
% calculates vertical km, number of cycles and float age using multi profile files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Floats: struct with al least two fields: WMO.data and dac.data
% dac_dir: path to gdac
%
% OUTPUT
% vertical_km: sum of vertical km (up and down) travelled by one float
%      during all its life time
% vertical_km_mean: mean vertical km (up and down) travelled by a float
%      in one cycle
% floats_age: float age in years from cycle 0
% last_cycle: number of cycles performed by float
%
% NOTES:
% (1) Vertical km are calculated using PRES variable in multiprofile file
%
% AUTHOR: Andrea Garcia Juan, Euro-Argo ERIC
%         (andrea.garcia.juan@euro-argo.eu)
%
% Modified on 2020/03/20

% Modified on 2024/09/01 Romain Cancouët:
% - add QC checks on vertical distance travelled
% - removed twice the vertical distance
% - removed descending profiles
% - add warnings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Calculating vertical km and floats age...')

n_floats = length(Floats.WMO);

% Init
vertical_km = NaN(n_floats,1);
vertical_km_mean = NaN(n_floats,1);
floats_age = NaN(n_floats,1);
last_cycle = NaN(n_floats,1);


%% floats loop
for ifloat = 1: n_floats
    
    disp(' ')
    fprintf('%s\n', Floats.WMO{ifloat})
    % floats directory path string
    dac = Floats.DAC{ifloat};
    % Open multiprofile file
    prof_file = [dac_dir '/dac/' char(dac) '/' Floats.WMO{ifloat} '/' Floats.WMO{ifloat} '_prof.nc'];
    
    try
        
        % Get N_PROF dimension
        vinfo = ncinfo(prof_file,'PRES');
        varSize = vinfo.Size;
        N_PROF = varSize(2);
        N_LEVELS = varSize(1);
        
        % Variables
        % For vertical_km
        pres = ncread(prof_file, 'PRES');
        pres_qc = ncread(prof_file, 'PRES_QC');
        pres_adj = ncread(prof_file, 'PRES_ADJUSTED');
        pres_adj_qc = ncread(prof_file, 'PRES_ADJUSTED_QC');
        p = NaN(N_LEVELS,N_PROF);
        p_qc = NaN(N_LEVELS,N_PROF);        
        
        % For age
        juld = ncread(prof_file, 'JULD');
        juld_qc = ncread(prof_file, 'JULD_QC');
        
        % Extra info to know what the correspondance between N_PROF and other variables
        direction = ncread(prof_file,'DIRECTION');
        cv = ncread(prof_file,'CYCLE_NUMBER');
        % vsc = ncread(prof_file,'VERTICAL_SAMPLING_SCHEME');
        
        %% GOOD PRESSURE VECTOR
        % Build a param vector "p" that is PARAM_ADJUSTED if exists or PARAM otherwise;
        % Do the same for QC
        % <PARAM>_ADJUSTED is mandatory. When no adjustment is performed, the FillValue is inserted.
        % if D files _ADJUSTED is filled, but also sometimes in R files
        for ip = 1:N_PROF
            % check if PARAM_ADJUSTED is full of NaN
            if all(isnan(pres_adj(:,ip)))
                p(:,ip) = pres(:,ip);
                p_qc(:,ip) = pres_qc(:,ip);
            else
                p(:,ip) = pres_adj(:,ip);
                p_qc(:,ip) = pres_adj_qc(:,ip);
            end
            % warning if empty primary profiles
            if all(isnan(p(:,ip)))
                warning(['primary profile of float ' Floats.WMO{ifloat} ' is empty at cycle ' num2str(cv(ip)) num2str(direction(ip))])
            end
        end
        
        % Minimum quality check
        % Avoid PRES_QC = 4 so replace by NaN
        bad_qc = p_qc == '4'; % Create a logical mask
        p(bad_qc) = NaN; % Replace all 4's in with NaN
        
        % Keep only ascent profiles
        i_ascending = find(direction == 'A');
        % Vector of maximum good pressures measured at each cycles (in km)
        vertical_m_vector = max(p(:,i_ascending))/1000;
        
        % Check if deep float
        PLATFORM_TYPE = ncread(prof_file,'PLATFORM_TYPE');
        deep = contains(PLATFORM_TYPE(:,1)','_D');
        if deep % for deep floats
            if any(vertical_m_vector > 6.2) % measured pressures deeper than 6200 dbar
                warning(['Float ' Floats.WMO{ifloat} ' has too deep PRES values'])
            end
        else % for core floats
            if any(vertical_m_vector > 2.2) % measured pressures deeper than 2200 dbar
                warning(['Float ' Floats.WMO{ifloat} ' has too deep PRES values'])
            end
        end
        
        % final step: compute vertical distance performed during a float life
        vertical_km(ifloat) = sum(vertical_m_vector,'omitnan'); % km
        vertical_km_mean(ifloat) = sum(vertical_m_vector,'omitnan')/length(vertical_m_vector); % average distance travelled in km
        
        %% GOOD DATE VECTOR
        % Minimum quality check (JULD_ADJUSTED is in the traj file
        % Avoid JULD_QC = 4 so replace by NaN
        bad_qc = juld_qc == '4'; % Create a logical mask
        juld(bad_qc) = NaN; % Replace all 4's in with NaN
        
        % floats_age(ifloat) = (max(juld,'omitnan') - min(juld,'omitnan'))/365; % Comment RC why max? prefer end and 1
        floats_age(ifloat) = (juld(find(~isnan(juld),1,'last'))-juld(find(~isnan(juld),1)))/365;
        last_cycle(ifloat) = cv(end);
        
        
    catch e
        % if file does not exist
        fprintf(2,'        %s\n', e.message)
        vertical_km(ifloat) = NaN;
        vertical_km_mean(ifloat)= NaN;
        floats_age(ifloat) = NaN;
        last_cycle(ifloat) = NaN;
    end   
    
end
