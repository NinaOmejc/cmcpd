% Initialize empty arrays to hold the data
task_list = {};
region_list = {};
group_list = {};
participant_list = [];
frequency_list = [];
value_list = [];

tasks = {'S', 'C', 'D'};
regions = {'avgpo', 'avgfcen', 'Cz'};
groups = {'H', 'P'};

% Iterate over tasks, regions, and groups to extract the data
for t = 1:length(tasks)
    task = tasks{t};
    for r = 1:length(regions)
        region = regions{r};
        for g = 1:length(groups)
            group = groups{g};
            data = wts.(task).(region).(group);
            [num_participants, num_frequencies] = size(data);
            for p = 1:num_participants
                for f = 1:num_frequencies
                    task_list{end+1} = task; %#ok<*AGROW>
                    region_list{end+1} = region;
                    group_list{end+1} = group;
                    participant_list(end+1) = p;
                    frequency_list(end+1) = f;
                    value_list(end+1) = data(p, f);
                end
            end
        end
    end
end

% Create a table and save it as a CSV file
T = table(task_list', region_list', group_list', participant_list', frequency_list', value_list', ...
    'VariableNames', {'Task', 'Region', 'Group', 'Participant', 'Frequency', 'Value'});
writetable(T, 'wts_all_zscore.csv');
