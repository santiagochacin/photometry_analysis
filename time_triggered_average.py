
import sys 
import pandas as pd
import numpy as np
import sys
sys.path.append('/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/photometry_analysis/photometry_3.0_module/')
import photometry_3 

def align_and_return_events(behavior_df, Fs, channel_to_trigger, channels_to_align=['channel_1_demod_F_F0', 'channel_2_demod_F_F0'], 
    threshold='Default', time_threshold='Default'):

    if channel_to_trigger == 'Acceleration' or channel_to_trigger ==  'Velocity':
        mean_ = np.nanmean(behavior_df[channel_to_trigger].values)
        std_ = np.nanstd(behavior_df[channel_to_trigger].values)
        print('mean is')
        print(mean_)
        print('std dev is')
        print(std_)
        if threshold=='Default':
            threshold = float(mean_) + float(4)*float(std_)
            print('Amplitude threshold is')
            print(threshold)
        if time_threshold=='Default':
            time_threshold = float(0)
            print('Time threshold is')
            print(time_threshold)
        params_df = pd.DataFrame({'amplitude_mean': [mean_], 'amplitude_std': [std_], 'amplitude_threshold': [threshold], 'time_threshold': [time_threshold]})
    elif channel_to_trigger == 'Grooming':
        if threshold == 'Default':
            threshold = 0.5
        if time_threshold == 'Default':
            time_threshold = 1
        params_df = pd.DataFrame({'amplitude_threshold': [threshold], 'time_threshold': [time_threshold]})
    regions = align_on_trigger(behavior_df, Fs, channel_to_trigger, channels_to_align,
        threshold, 'up', time_threshold, 5, (0, 0))
    
    unsorted_regions = return_events(regions)

    return(unsorted_regions, params_df)


def count_events_in_array_threshold(selection, Fs, time_to_wait, threshold=0, up=True):
    """ Counts the number of events (e.g action potentials (AP)) in the current trace.
    Arguments:
    trace -- 1d numpy array
    dt -- sampling interval
    threshold -- (optional) detection threshold (default = 0).
    up -- (optional) True (default) will look for upward events, False downwards.

    Returns:
    An integer with the number of events.
    Examples:
    count_events(500,1000) returns the number of events found between t=500 ms and t=1500 ms
    above 0 in the current trace and shows a stf marker.
    count_events(500,1000,0,False,-10,i) returns the number of events found below -10 in the 
    trace 
    """
    # add in a minimum bout length?
    time_threshold_samples = int(time_to_wait*Fs)
    # algorithm to detect events
    EventCounter,i = 0,0 # set counter and index to zero
    # list of sample points
    sample_points = []
    # choose comparator according to direction:
    if up:
        comp = lambda a, b: a > b
    else:
        comp = lambda a, b: a < b
    # run the loop
    while i<len(selection):

        if comp(float(selection[i]),float(threshold)):

            to_add = i
            bout_start = i

            while i<len(selection) and comp(selection[i],threshold):
               i+=1 # skip values if index in bounds AND until the value is below/above threshold again
       
            bout_end = i 
            if time_threshold_samples < (bout_end - bout_start):
                EventCounter +=1
                sample_points.append(bout_start)   

        else:
            i+=1

    return (EventCounter, sample_points)

def measure_bout_lengths(selection, start_index=0, end_index='end', threshold=0, up=True):
    if end_index =='end':
        end_index = int(len(selection-1))
    # algorithm to detect events
    i = start_index # set index to zero
    # list of sample points
    sample_points = []
    bout_lengths = []
    # choose comparator according to direction:
    if up:
        comp = lambda a, b: a > b
    else:
        comp = lambda a, b: a < b
    # run the loop
    while i<end_index:
        if comp(float(selection[i]),float(threshold)):
            to_add = i
            bout_start = i
            while i<end_index and comp(selection[i],threshold):
               i+=1 # skip values if index in bounds AND until the value is below/above threshold again
            bout_end = i 
            sample_points.append(bout_start)   
            bout_lengths.append(bout_end - bout_start)

        else:
            i+=1

    return (sample_points, bout_lengths)

def return_indicies(df_name, Fs, channel_to_trigger, amplitude_threshold, time_threshold, window_time, up):

    """
    Output: numpy array of sample points from input channel where trace passes above threshold
    """
    indicies_from_trigger_channel_by_sweep = {}
    for sweep in sorted(list(set(df_name.index.get_level_values(0)))):
        trace = df_name.loc[sweep][channel_to_trigger].values
        event_number, sample_indicies = count_events_in_array_threshold(trace,
                                                                 Fs, 
                                                                 time_threshold, 
                                                                 amplitude_threshold, up)
        # remove sample indicies outside of "window"
        sample_indicies_to_add = [i for i in sample_indicies if (window_time[0]*Fs < i < (int(len(trace)) - window_time[1]*Fs))]
        
        if event_number > 0:
            indicies_from_trigger_channel_by_sweep[sweep] = sample_indicies_to_add

    return(indicies_from_trigger_channel_by_sweep)

def align_on_trigger(df_name, Fs, channel_to_trigger, channels_to_align, amplitude_threshold, direction, time_threshold, 
                     time_region, window_time=(0,0)):
    """
    Output: dataframe of trace regions surrounding threshold for all channels
    """
    
    # call return indicies to get sample indicies
    if direction == 'up':
        up = True
    elif direction == 'down':
        up = False 
    indicies_from_trigger_channel = return_indicies(df_name, Fs, channel_to_trigger, 
                                                    amplitude_threshold, time_threshold, window_time, direction)
    
    by_channel_by_sweep = {}
    region_in_samples = int(time_region*Fs)
    
    channels_to_align.append(channel_to_trigger)
    for channel in channels_to_align:
        for sweep in sorted(list(set(df_name.index.get_level_values(0)))):
            if sweep in list(indicies_from_trigger_channel.keys()):
                
                by_index = {}
                for i in indicies_from_trigger_channel[sweep]:
            
                    
                    if int(i-region_in_samples) < 0:
                        trace_region = np.pad(df_name.loc[sweep][channel].values[0:int(i+region_in_samples)], 
                                              (int(region_in_samples)*2-int(i+region_in_samples), 0), 'constant', constant_values=np.nan)
                        print(str(len(trace_region)) + 'earlier')
                    
                    elif int(i+region_in_samples) > len(df_name.loc[sweep][channel].values):
                        trace_region = np.pad(df_name.loc[sweep][channel].values[int(i-region_in_samples):], 
                                              (0, int(i+region_in_samples)-len(df_name.loc[sweep][channel].values))
                                              , 'constant', constant_values=np.nan)
                        print(str(len(trace_region)) + 'later')
                    
                    else:
                        trace_region = df_name.loc[sweep][channel].values[int(i-region_in_samples):int(i+region_in_samples)]
                        print(str(len(trace_region)) + 'middle')
                        
                    by_index[i] = trace_region
                    
                index = list(range(int(-region_in_samples), int(region_in_samples)))
        
                by_index = pd.DataFrame(by_index)
                                                                
                by_channel_by_sweep[(channel, sweep)] = by_index
            else:
                pass

    return(by_channel_by_sweep)
        
def sort_events(array, sort_type, sort_amp_threshold, sort_time_threshold, *argv):
    """
    sort_type is string     e.g  'before' : before_amp_threshold, before_time_threshold
                                 'after'  : after_amp_threshold, after_time_threshold
                                 'before&after' : (before_thresold, after_threshol)}
    args depend on sort type e.g if 'before' or 'after' argv[0] is amp threshold
                                                        argv[1] is time threshold
                                 if 'before&after' argv[0] & [1] are before threshold
                                                   argv[2] & [3] are after threshold
    """
    threshold_index = int(len(array)/2)
    output = True
    if sort_type == 'before':
        for item in array[(threshold_index-sort_time_threshold):threshold_index]:
            if item > sort_amp_threshold: 
                output = False
                break

    elif sort_type == 'after':
        for item in array[threshold_index+1:(threshold_index+sort_time_threshold)]:
            if item < sort_amp_threshold: 
                output = False
                break

    elif sort_type == 'before&after':
        for item in array[(threshold_index-sort_time_threshold):threshold_index]:
            if item > sort_amp_threshold: 
                output = False
                break

        for item in array[threshold_index+1:(threshold_index+argv[0][1])]:
            if item < argv[0][0]: 
                output = False
                break

    return(output)

        
def average_sweep_responses(dict_of_dfs):
    
    channels = list(set((list(regions.keys()))[i][0] for i in range(len(regions.keys()))))
    sweeps = list(set((list(regions.keys()))[i][1] for i in range(len(regions.keys()))))
    
    channel_means = {}
    
    for channel in channels:
        list_to_concat = []
        for sweep in sweeps:
            for key in list(dict_of_dfs.keys()):
                if channel in key:
                    list_to_concat.append(regions[key].mean(axis=1))
        
        channel_means[channel] = pd.concat(list_to_concat, axis=1).mean(axis=1)
        
    mean_df = pd.DataFrame(channel_means)
    
    return(mean_df)

def return_events(dict_of_dfs):
    
    channels = list(set((list(dict_of_dfs.keys()))[i][0] for i in range(len(dict_of_dfs.keys()))))
    sweeps = list(set((list(dict_of_dfs.keys()))[i][1] for i in range(len(dict_of_dfs.keys()))))
    
    concat_by_channels = []
    
    for channel in channels:
        concat_by_sweeps = []
        for sweep in sweeps:
            to_load = (channel, sweep)
            concat_by_sweeps.append(dict_of_dfs[to_load])
        concat_by_channels.append(pd.concat(concat_by_sweeps, keys=sweeps, axis=1))
    concat_by_channels_df = pd.concat(concat_by_channels, keys=channels, axis=1)

    return(concat_by_channels_df)
        
                
    

def return_sorted_events(unsorted_regions, sort_column, sort_type, sort_amp_threshold, sort_time_threshold, *argv):

    # store a list of tuples corresponding to the sweep and sample index of the event passing the sort test
    sorted_events = {}
    for event in range(len(unsorted_regions[sort_column].columns)):
        if sort_events(unsorted_regions[sort_column].iloc[:, event].values, sort_type, sort_amp_threshold, sort_time_threshold, argv) == True:
            sorted_events[unsorted_regions[sort_column].iloc[:, event].name] = ''

    # create data frame with sorted events
    sorted_by_channels = {}
    for channel in unsorted_regions.columns.get_level_values(level=0):
        by_sweeps = {}
        for sweep in list(set([key[0] for key in sorted_events.keys()])):
            events_in_sweep = {}
            for event in list(set([key[1] for key in sorted_events.keys() if key[0] == sweep])):
                events_in_sweep[event] = unsorted_regions[channel][sweep][event]
            by_sweeps[sweep] = pd.DataFrame(events_in_sweep)
        sorted_by_channels[channel] = pd.concat(list(by_sweeps.values()), keys=list(by_sweeps.keys()), axis=1)
        
    all_concat = pd.concat(list(sorted_by_channels.values()), keys=list(sorted_by_channels.keys()), axis=1)

    return(all_concat)


