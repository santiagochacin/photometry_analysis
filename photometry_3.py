import neo 
from neo import io
import pandas as pd
import numpy as np
import scipy.ndimage
from scipy import signal 
from scipy import stats
import matplotlib.pyplot as plt
import tables

## classes to store photometry data

class raw_photometry_file(object):
	## create methods to act on df file (whether a modulated or not modulated signal) for ease of plotting etc...
	def __init__(self, filename, encoder_bin, rotary_output=True):
		#loads an abf file to a df as initialization step
		if filename.endswith('.abf'):
			self.class_df = abf_file_to_df(filename)
		elif filename.endswith('.h5'):
			self.class_df = pd.read_hdf(filename)
		#store sampling interval info 
		self.dt = self.class_df.loc['sweep001']['time'][1]-self.class_df.loc['sweep001']['time'][0]
		self.Fs = 1/self.dt
		self.rotary_output = rotary_output
		#creates a df with encoder data transformed to continuous output if rotary_output=True (default)
		if rotary_output == True and filename.endswith('with_velocity_data.h5'):
			self.cont_df = self.class_df
		elif rotary_output == True and 'with_velocity_data' not in filename:
			self.cont_df = transform_encoder_to_continuous(self.class_df, ['channel_1', 'channel_2'], self.dt, encoder_bin)
			
		
		
	def demodulate_file(self, df_to_demod, demod_params):
		"""
		Inputs:
			demod_params = dictionary with inputs for demod file function
							{ 
							encoder_channel: string e.g. 'channel_1'
							fluorescence_channel: string e.g. 'channel_5'
							modulation_channel_1: string e.g. 'channel_3'
							modulation_channel_2: string e.g. 'channel_4'
							pre_filter_params_ch1: dictionary e.g {'filter_type': 'lowpass' or 'bandpass' ,
    															   'order': integer
    															   'highcut': float --this will apply for both a lowpass and bandpass filter ,
    															   'lowcut': float --will only be used if filter type is 'bandpass' , 
    																}
							pre_filter_params_ch2: dictionary same as for ch_1 
							post_filter_params: 
							window_time: tuple of seconds to window from start and end e.g. (6,1)
							}
			
			file_df = self.class_df
			dt_file = self.dt
			Fs_file = self.Fs
		
		Outputs:
			
			a new instance of demodculated_photometry_file_class 
		
		"""
		demodulated_file_class_output = demodulated_photometry_file(df_to_demod, demod_params, self.dt, self.Fs)
		return(demodulated_file_class_output)
		
	def save_to_excel(self, filename):
		self.class_df.to_excel(filename+'raw_data.xlsx')
		self.cont_df.to_excel(filename+'with_velocity_data.xlsx')
		return()

	def save_to_hdf(self, filename):
		#can also add to an existing file
		#see https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.to_hdf.html
		self.class_df.to_hdf(filename+'raw_data'+'.h5', key='df', mode='w')
		if self.rotary_output == True:
			self.cont_df.to_hdf(filename+'with_velocity_data'+'.h5', key='df', mode='w')
		return(True)
		
class demodulated_photometry_file(object):
	def __init__(self, file_df, demod_params, dt_input, Fs_input, from_file=False):
		
		#sampling interval info passed
		if from_file == False:
			self.dt = dt_input
			self.Fs = Fs_input
		
		#demodulated df
		if from_file == False:
			self.class_df = demod_file(file_df, demod_params['encoder_channel'], demod_params['fluorescence_channel_1'], demod_params['modulation_channel_1'], demod_params['modulation_channel_2'],
				self.dt, self.Fs, demod_params['pre_filter_params_ch1'], demod_params['pre_filter_params_ch2'], demod_params['post_filter_params'], demod_params['seconds_to_window'], demod_params['in_quad'])
		elif from_file == True: 
			self.class_df = pd.read_hdf(file_df)
		
	def save_to_excel(self, filename):
		self.class_df.to_excel(filename+'demodulated.xlsx')
		return()

	def save_to_hdf(self, filename):
		#can also add to an existing file
		#see https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.to_hdf.html
		self.class_df.to_hdf(filename+'demodulated'+'.h5', key='df', mode='w')
		return(True)

class photometry_df_from_file(object):
	def __init__(self, file_name):
		self.class_df = pd.read_hdf(file_name)
		self.dt = self.class_df.loc['sweep001']['time'].values[1]-self.class_df.loc['sweep001']['time'].values[0]
		self.Fs = 1/self.dt
	

class transform_df(object):
	def __init__(self, class_name, channels_to_calc, window_time, transform_params):
		""" 
		transform_params is a dictionary that stores manipulation and params for that manipulation
		e.g. for a gaussian filter
		transform_params = {'gaussian_filter_signal' : [order]}
		"""
		self.class_df = transform_func(class_name, channels_to_calc, window_time, transform_params)
		self.dt = class_name.dt
		self.Fs = class_name.Fs

	#create an output to save to excel
	def save_to_excel(self, filename):
		self.class_df.to_excel(filename+'demodulated.xlsx')
		return()	

	def save_to_hdf(self, filename):
		#can also add to an existing file
		#see https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.to_hdf.html
		self.class_df.to_hdf(filename+'.h5', key='df', mode='w')
		return(True)	


class calc_variance_of_df(object):
	def __init__(self, class_name, df_input, channels_to_calc_variance, window_time):
		self.df_with_variance = get_variance_of_data_frame(class_name, df_input, channels_to_calc_variance, window_time)
		self.dt = class_name.dt
		self.Fs = class_name.Fs

	def save_to_excel(self, filename):
		self.df_with_variance.to_excel(filename+'variance.xlsx')
		return()

## load abf file with raw data and output F/F0

def demod_abf_file(path_to_abf_file, demod_params, encoder_bin=.1, encoder_data=False, 
	channels_to_calc=['channel_1_demod', 'channel_2_demod'], window_time=(5,5), sliding_baseline_length=2):

	recording_class = raw_photometry_file(path_to_abf_file, encoder_bin, encoder_data)
	fname_for_saving = path_to_abf_file.rstrip('.abf')
	recording_class.save_to_hdf(fname_for_saving)
	recording_demod = recording_class.demodulate_file(recording_class.class_df, demod_params)
	recording_demod.save_to_hdf(fname_for_saving)
	print('calculating F/F0')
	recording_demod_F_F0 = transform_df(recording_demod, channels_to_calc, window_time, {'F_F0': sliding_baseline_length})
	recording_demod_F_F0.save_to_hdf(fname_for_saving+'_F_F0') 

	return(recording_demod_F_F0)



## filter functions
# convert  bf_Fc_input from Hz to 0 to 1 (normalized to nyquist frequency, which is 1/2 sampling rate)
#		Fc_input_ch1 = bf_Fc_input_channel_1/(self.Fs/2)
#		Fc_input_ch2 = bf_Fc_input_channel_2/(self.Fs/2)


def butterworth_bandpass_filter(order, Fc_low, Fc_high, signal_unfilt, sampling_rate):
	# convert  input cutoff frequencies from Hz to 0 to 1 (normalized to nyquist frequency, which is 1/2 sampling rate)
	Fc_low_norm = Fc_low/(sampling_rate/2)
	Fc_high_norm = Fc_high/(sampling_rate/2) 
	b, a = signal.butter(order, [Fc_low_norm, Fc_high_norm], btype = 'band')
	signal_filtered = signal.filtfilt(b, a, signal_unfilt)
	return(signal_filtered)


def butterworth_lowpass_filter(order, Fc, signal_unfilt, sampling_rate):
	"""Butterworth filter signal using signal from scipy.signal library"""
	# convert  input cutoff frequencies from Hz to 0 to 1 (normalized to nyquist frequency, which is 1/2 sampling rate)
	Fc = Fc/(sampling_rate/2)
	b, a = signal.butter(order, Fc, 'low')
	zi = signal.lfilter_zi(b, a)
	z, _ = signal.lfilter(b, a, signal_unfilt, zi=zi*signal_unfilt[0])
	z2, _ = signal.lfilter(b, a, z, zi=zi*signal_unfilt[0])
	
	signal_filtered = signal.filtfilt(b, a, signal_unfilt)
	
	return(signal_filtered)
	
def gaussian_filter_signal(signal_to_filter, order):
	"""simple gaussian smoothing"""

	sigma=order
	signal_gaussian_filtered = scipy.ndimage.filters.gaussian_filter(signal_to_filter, sigma)
	
	return(signal_gaussian_filtered)

def gaussian_filter_df(df, order):

	dict_ = {}
	for column in list(df.columns):
		dict_[column] = scipy.ndimage.filters.gaussian_filter(df[column].values, order)
	df_ = pd.DataFrame(dict_, index=df.index)
	return(df_)

	
## misc useful functions 
def count_events_in_array(trace, Fs, time_to_wait, threshold=0, up=True):
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
	samples_to_skip = int(time_to_wait*Fs)
	selection = trace
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
			EventCounter +=1
			sample_points.append(i)

			i += 1 #skip set number of samples from threshold crossing
			
			while i<len(selection) and comp(selection[i],threshold):
				i+=1 # skip values if index in bounds AND until the value is below/above threshold again
		else:
			i+=1

	return (EventCounter, sample_points)

def get_derivative(x_values, y_values):
	derivative = []

	for item in range(len(x_values)-1):
		derivative.append((y_values[item+1]-y_values[item])/(x_values[item+1]-x_values[item]))

	derivative_out = np.pad(derivative, (0,1), 'constant', constant_values=np.nan)

	return(derivative_out)


def sort_events(array, sort_type, sort_amp_threshold, sort_time_threshold, *argv):
	"""
	sort_type is key value  e.g  'before' : before_amp_threshold, before_time_threshold
								 'after'  : after_amp_threshold, after_time_threshold
								 'before&after' : (before_thresold, after_threshol)}
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

		for item in array[threshold_index+1:(threshold_index+argv[1])]:
			if item < argv[0]: 
				output = False
				break

	return(output)
	
## functions for manipulating data sweeps

def window_channel(data_input_array_unwindowed, dt, start_window_time, end_window_time):
	"""windows a sweep of data
	   output: numpy array"""
	start_window_samples = int(start_window_time/dt)
	end_window_samples = int(end_window_time/dt)
	data_output_array_windowed = data_input_array_unwindowed[start_window_samples:(len(data_input_array_unwindowed-end_window_samples))]
	
	return(data_output_array_windowed)
	
def get_z_score(input_data_array):
	"""calculates variance 
	   output: numpy array"""
	df_mean = np.mean(input_data_array)
	df_std = np.std(input_data_array)
	output_data_array = [((point-df_mean)/df_std) for point in input_data_array]
	
	return(output_data_array)

## encoder specific functions

def get_velocity(encoder_input, dt, velocity_bin_in_time):
	"""takes square pulses from an encoder channel and creates continous ouput by binning
	   output: numpy array"""
	bin_length_in_samples = int(velocity_bin_in_time*(1/dt))
	running_count = 0 
	pulses_binned = []
	for selection_start in range(0, np.shape(encoder_input)[0], bin_length_in_samples):
		rotary_outout = count_events_in_array(encoder_input[selection_start:selection_start+bin_length_in_samples], dt, 0 ,2, True)[0]
		pulses_binned.extend((rotary_outout*x) for x in [1] * bin_length_in_samples)
		pulses_binned_out = np.array(pulses_binned)

	return(pulses_binned_out)

## fluorescence specific functions

def demodulate_signal(signal_to_demod, mod_waveform, Fs, Fc_post_demod, post_demod_order):
	"""demodulate signal: multiply by lock-in reference and then lowpass filter"""
	#"manual" lock in multiplication
	lock_in_signal = signal_to_demod*mod_waveform
	#lowpass filter 
	filtered_lock_in = butterworth_lowpass_filter(post_demod_order, Fc_post_demod, lock_in_signal, Fs)
	
	return(filtered_lock_in)

def synthesize_and_demodulate_in_quadrature(file_df, mod_Hz, Fs, window, sweep, fluorescence_signal, fluorescence_demod, modulation_channel, post_filter_params):
	"""
		window: tuple indicating samples (window_start, window_end)
	"""
	##steps to demodulate "in quadrature" with synthesized signal 
	#need to create a separate modulation channel phase shifted by 90 degrees for recombining "in quadrature"
	#samples_sec = _recording.Fs
	#sampling_interval = _recording.dt
	#frequency of modulation is 211Hz, 211 pulses/sec, so 1/211 seconds/cycle
	#shift time base by 1/2 cycle or 1/422 seconds
	#in samples 10000 samples/second, so (10000/422) samples/half cycle 
	samples_to_shift = int(int(Fs)/(int(mod_Hz*2)))
	#shift and windowed modulation signal 
	modulation_signal_shifted_windowed = file_df.loc[sweep][modulation_channel].values[int(window[0]+samples_to_shift):int(window[1]+samples_to_shift)]
	demodulated_quadrature = modulation_signal_shifted_windowed*fluorescence_signal
	sum_in_quadrature = (fluorescence_demod**2+demodulated_quadrature**2)**(1/2)
	post_demod_filtered_sum_quad = butterworth_lowpass_filter(order=post_filter_params['order'], 
																Fc=post_filter_params['highcut'], 
																signal_unfilt=sum_in_quadrature, 
																sampling_rate=Fs)
	return (post_demod_filtered_sum_quad)


def F_f0(fluorescence_power, Fs, sliding_baselne_length):
	"""calculate F/F0 as a sliding baseline
		Inputs: sliding_baseline_length (in seconds)
	"""
	baseline_length_samples = int(Fs*sliding_baselne_length)
	F_F0_array = []
	for point in range(int(baseline_length_samples), len(fluorescence_power)):
		F0 = np.mean(fluorescence_power[int(point-baseline_length_samples):int(point)])
		F_F0 = fluorescence_power[point]/F0
		F_F0_array.append(F_F0)

	return(np.array(F_F0_array))


## functions for loading/saving files
def abf_file_to_df(path_to_abf_file):
	"""Inputs: file path to abf file
	   Outputs: multidimensional data fame containing data organized by sweeps and channels
	"""
	
	r = io.AxonIO(filename=path_to_abf_file)
	bl = r.read_block(lazy=False)
	num_channels = len(bl.segments[0].analogsignals)
	channels = []
	df_list = []
	signals = []
	sweep_list = []

	for seg_num, seg in enumerate(bl.segments):
		channels = ['channel_{0}'.format(str(i)) for i in range(num_channels)]
		signals = []
		
		for i in range(num_channels):
			signals.append(np.array(bl.segments[seg_num].analogsignals[i])[:,0])
			
		data_dict = dict(zip(channels, signals))
		time = seg.analogsignals[0].times - seg.analogsignals[0].times[0]
		data_dict['time'] = time
		df = pd.DataFrame(data_dict)
		df_list.append(df)
		sweep_list.append('sweep' + str(seg_num + 1).zfill(3))
		
	df = pd.concat(df_list, keys=sweep_list, names=['sweep'])
	
	return(df)

## plotting functions	

def plot_multiple_sweeps_channels(df_to_plot, sweeps_list, channels_list, window_samples=(0,0), 
						wspace=0.1, hspace=0.1, X_axis_column=None, X_label=None, Y_labels=None):
	"""
	Inputs: list of sweeps (index, level 0),
			list of channels (columns)
			time(samples) to window from each sweep
			Optional:
			wspace, hspace, column from df_to use for X axis, X axis label (applies to whole plot)
			list of Y labels to use (list 1 for each subplot)
			
	"""
	plt.subplots_adjust(wspace, hspace)
	# sweeps/index level 0 plotted by column
	column_num = int(len(sweeps_list))
	# channels are plotted by row 
	row_num = int(len(channels_list))
	subplot_index = 1
	for channel_to_plot in range(len(channels_list)):
		for sweep in range(len(sweeps_list)):
			plt.subplot(row_num, column_num, subplot_index)
			if X_axis_column:
				plt.plot(df_to_plot.loc[str(sweeps_list[sweep])][str(X_axis_column)].as_matrix(),
						df_to_plot.loc[str(sweeps_list[sweep])][str_channels_list[channel_to_plot]].as_matrix())
			else:
				plt.plot(df_to_plot.loc[str(sweeps_list[sweep])][str_channels_list[channel_to_plot]].as_matrix())
			if X_label:
				plt.xlabel(X_label)
			if Y_labels:
				plt.ylabel(str(Y_labels[subplot_index]))
			else:
				plt.ylabel(str(channels_list[channel_to_plot]))
			subplot_index += 1
	plt.show()
	return()
	
	
## functions for transforming dataframes created from abf files

def demod_file(file_df, encoder_channel_file, fluorescence_channel_1_file, modulation_channel_1_file, 
				modulation_channel_2_file, dt_file, Fs_file, pre_filter_params_ch1, pre_filter_params_ch2, post_filter_params, window_time, in_quad=True):
	labels_list = sorted(list(set(file_df.index.get_level_values(0))))
	sweep_dfs = []
	for sweep in labels_list:
		print('demodulating sweep', sweep)
		sweep_dfs.append(demod_sweep(file_df, sweep, encoder_channel_file, fluorescence_channel_1_file, modulation_channel_1_file, 
									modulation_channel_2_file, dt_file, Fs_file, pre_filter_params_ch1, pre_filter_params_ch2, post_filter_params ,window_time, in_quad))
		print('sweep', sweep, 'done')
	print('concactenating sweeps')	
	demod_file_df = pd.concat(sweep_dfs, keys=labels_list, names=['sweep'])

	return(demod_file_df)

def demod_sweep(file_df, sweep, encoder_channel, fluorescence_channel_1, modulation_channel_1, modulation_channel_2, dt, Fs, pre_filter_params_ch1, pre_filter_params_ch2, post_filter_params,
	 window_time_sweep, in_quad=True):    
    """
    Inputs:
    	pre_filter_params = dictionary e.g {'filter_type': 'lowpass' or 'bandpass' ,
    									'order': integer
    									'highcut': float --this will apply for both a lowpass and bandpass filter ,
    									'lowcut': float --will only be used if filter type is 'bandpass' , 
    									}
    	
    	post_filter_params = dictionary e.g {'filter_type': 'lowpass' or 'bandpass' ,
    									'order': integer
    									'highcut': float --this will apply for both a lowpass and bandpass filter ,
    									'lowcut': float --will only be used if filter type is 'bandpass' , 
    									}							
    """
    #windowing
    windowed_start = int(window_time_sweep[0]*Fs)
    windowed_end = len(file_df.loc['sweep001']['time']) - int(window_time_sweep[1]*Fs)

    # for each sweep
    sweep_dict = {}
    #add time column
    sweep_dict['time'] = file_df.loc[sweep]['time'].as_matrix()[windowed_start:windowed_end]
    
    # window encoder
    encoder_windowed = file_df.loc[sweep][encoder_channel].as_matrix()[windowed_start:windowed_end]
    sweep_dict['encoder_windowed'] = encoder_windowed 

    # window fluorescence channel 1
    fluorescence_channel_1_windowed = file_df.loc[sweep][fluorescence_channel_1].as_matrix()[windowed_start:windowed_end]
    sweep_dict['fluorescence_channel_1_windowed'] = fluorescence_channel_1_windowed 

    # window modulation for channel 1
    modulation_channel_1_windowed = file_df.loc[sweep][modulation_channel_1].as_matrix()[windowed_start:windowed_end]

	# need to create a separate modulation channel phase shifted by 90 degrees for recombining "in quadrature"
	    
    # window modulation for channel 2
    modulation_channel_2_windowed = file_df.loc[sweep][modulation_channel_2].as_matrix()[windowed_start:windowed_end]
    
    # need to create a separate modulation channel phase shifted by 90 degrees for recombining "in quadrature"
    
    ## demodulate
    #channel 1
    #prefilter
    print('prefilter channel 1')
    if pre_filter_params_ch1['filter_type'] == 'lowpass' :
    	channel_1_filt = butterworth_lowpass_filter(pre_filter_params_ch1['order'], pre_filter_params_ch1['highcut'], fluorescence_channel_1_windowed, Fs)
    
    elif pre_filter_params_ch1['filter_type'] == 'bandpass':
    	channel_1_filt = butterworth_bandpass_filter(pre_filter_params_ch1['order'], pre_filter_params_ch1['lowcut'], pre_filter_params_ch1['highcut'], fluorescence_channel_1_windowed, Fs)
    
    #demod
    print('demodulating channel 1')
    channel_1_demod = demodulate_signal(channel_1_filt, modulation_channel_1_windowed, Fs, post_filter_params['highcut'], post_filter_params['order'])
    if in_quad==False:
    	sweep_dict['channel_1_demod'] = channel_1_demod 

    elif in_quad==True:
    	print('demodulating/summing in quadrature channel 1')
    	# entering modulation frequency manually - can also recover from FFT
    	#									synthesize_and_demodulate_in_quadrature(mod_Hz, Fs, window, sweep, fluorescence_signal, fluorescence_demod, modulation_channel, post_filter_params)
    	post_demod_filtered_488_sum_quad = synthesize_and_demodulate_in_quadrature(file_df, 211, Fs, (windowed_start, windowed_end), sweep, channel_1_filt, channel_1_demod, modulation_channel_1, post_filter_params)
    	sweep_dict['channel_1_demod'] = post_demod_filtered_488_sum_quad 

	#channel 2 / second modulation 
    #prefilter
    print('prefilter channel 2')
    if pre_filter_params_ch2['filter_type'] == 'lowpass' :
    	channel_2_filt = butterworth_lowpass_filter(pre_filter_params_ch2['order'], pre_filter_params_ch2['highcut'], fluorescence_channel_1_windowed, Fs)
    
    elif pre_filter_params_ch2['filter_type'] == 'bandpass':
    	channel_2_filt = butterworth_bandpass_filter(pre_filter_params_ch2['order'], pre_filter_params_ch2['lowcut'], pre_filter_params_ch2['highcut'], fluorescence_channel_1_windowed, Fs)
    
    #demod
    print('demodulating channel 2')
    channel_2_demod = demodulate_signal(channel_2_filt, modulation_channel_2_windowed, Fs, post_filter_params['highcut'], post_filter_params['order'])
    if in_quad==False:
    	sweep_dict['channel_2_demod'] = channel_1_demod 

    elif in_quad==True:
    	print('demodulating/summing in quadrature channel 2')
    	# entering modulation frequency manually - can also recover from FFT
    	#synthesize_and_demodulate_in_quadrature(file_df, mod_Hz, Fs, window, sweep, fluorescence_signal, fluorescence_demod, modulation_channel, post_filter_params)
    	post_demod_filtered_405_sum_quad = synthesize_and_demodulate_in_quadrature(file_df, 516, Fs, (windowed_start, windowed_end), sweep, channel_2_filt, channel_2_demod, modulation_channel_2, post_filter_params)
    	sweep_dict['channel_2_demod'] = post_demod_filtered_405_sum_quad 

    sweep_df = pd.DataFrame(sweep_dict)
    
    return(sweep_df)

def window_df(unwindowed_data_frame, dt, start_window_time, end_window_time):
	"""cuts indicated amount of time (in sec) from start and end of data frame"""
	
	start_window_samples = int(start_window_time/dt)
	end_window_samples = int(end_window_time/dt)
	
	windowed_dfs = []
	sweeps = []
	for sweep in sorted(list(set(unwindowed_data_frame.index.get_level_values(0)))):
		windowed_dfs.append(unwindowed_data_frame.loc[sweep].iloc[start_window_samples:(len(unwindowed_data_frame.loc[sweep])-end_window_samples),:])
		sweeps.append(sweep)
	
	windowed_df_out = pd.concat(windowed_dfs, keys=sweeps, names=['sweep'])
	
	return(windowed_df_out)

def transform_encoder_to_continuous(df_input, list_of_encoder_channels, dt, encoder_bin):
	"""takes list of encoder channels (e.g. 'channel_1', 'channel_2')
	returns a dataframe with continuous encoder data substituted in"""
	
	continuous_encoder_dfs = []
	sweeps = []
	for sweep in sorted(list(set(df_input.index.get_level_values(0)))):
		#transform to continuous trace 
		print('reading encoder for', sweep)
		encoder_A_continuous = get_velocity(df_input.loc[sweep][list_of_encoder_channels[0]].as_matrix(), 
		dt, encoder_bin)
		encoder_B_continuous = get_velocity(df_input.loc[sweep][list_of_encoder_channels[1]].as_matrix(), 
		dt, encoder_bin)
		
		angular_velocity_A = angular_velocity = [((encoder_A_continuous[i]*(0.00031415926))/.1) for i in range(len(encoder_A_continuous))]
		angular_velocity_B = angular_velocity = [((encoder_B_continuous[i]*(0.00031415926))/.1) for i in range(len(encoder_B_continuous))]

		continuous_encoder_dfs.append(pd.DataFrame({'encoder_channel_A':encoder_A_continuous, 'encoder_channel_B':encoder_B_continuous,
													'angular_velocity_A':angular_velocity_A,'angular_velocity_B':angular_velocity_B},
												index=np.array(df_input.loc['sweep001','channel_0':'channel_0'].index.tolist())))
												
		sweeps.append(sweep)
	
	print('concactenating')
	continuous_encoder_df = pd.concat(continuous_encoder_dfs, keys=sweeps, names=['sweep'])
	df_out = pd.concat([df_input.loc[:, 'time'],df_input.loc[:, 'channel_0':'channel_0'], 
						continuous_encoder_df.loc[:, 'encoder_channel_A':'angular_velocity_B'],
						df_input.loc[:, 'channel_3':'channel_5']], axis=1)
						
	
	return(df_out)
	

def transform_func(class_name, channels_to_calc, window_time, transform_params):
	"""return a df with columns added for variance info for each sweep loaded from .abf file
	Inputs:
	photometry class object
	window (tuple of time values"""
	
	df_input = class_name.class_df
	
	channels_list = df_input.columns.values
	windowed_start = int(window_time[0]*class_name.Fs)
	windowed_end = len(df_input.loc['sweep001'][str(channels_list[0])]) - int(window_time[1]*class_name.Fs)
	labels_list = sorted(list(set(df_input.index.get_level_values(0))))
	
	sweeps = []
	sweep_data = []
	for sweep in labels_list:
		print(sweep)
		sweep_dict = {}
		#window & add data
		for column_nv in list(channels_list):
			sweep_dict[column_nv] = df_input.loc[sweep][column_nv].as_matrix()[windowed_start:windowed_end]
		
		
		array_length = len(df_input.loc[sweep]['time'].as_matrix()[windowed_start:windowed_end])	
		for column_v in channels_to_calc:
		#transform indicated channels
			transformation = list(transform_params.keys())[0]
			#loop with ways to transform 
			if transformation == 'gaussian filter':
				order = list(transform_params.values())[0]
			
				toadd = gaussian_filter_signal(df_input.loc[sweep][column_v].as_matrix()[windowed_start:windowed_end], 
					order)
				#need to pad initial baseline period to make same length as other columns
				padded = np.pad(toadd, (array_length-len(toadd), 0), 'edge')	
				sweep_dict[str(column_v)+'_'+str(transformation)] = padded 
			
			if transformation == 'z-score':
				order = list(transform_params.values())[0]
			
				toadd = get_z_score(df_input.loc[sweep][column_v].as_matrix()[windowed_start:windowed_end])
				#need to pad initial baseline period to make same length as other columns
				padded = np.pad(toadd, (array_length-len(toadd), 0), 'edge')	
				sweep_dict[str(column_v)+'_'+str(transformation)] = padded 
			
			if transformation == 'F_F0':
				Fs = class_name.Fs
				sliding_baseline_length = list(transform_params.values())[0]
			
				toadd = F_f0(df_input.loc[sweep][column_v].as_matrix()[windowed_start:windowed_end], Fs, sliding_baseline_length)
				#need to pad initial baseline period to make same length as other columns
				padded = np.pad(toadd, (array_length-len(toadd), 0), 'edge')	
				sweep_dict[str(column_v)+'_'+str(transformation)] = padded 
			
		sweep_data.append(pd.DataFrame(sweep_dict))
		sweeps.append(sweep)
		
	file_output = pd.concat(sweep_data, keys=sweeps, names = ['sweep'])
		
	return(file_output)

## functions for aligning data to events in different channels

def return_indicies(class_name, channel_to_trigger, amplitdue_threshold, time_threshold):
	"""
	Output: numpy array of sample points from input channel where trace passes above threshold
	"""
	indicies_from_trigger_channel_by_sweep = {}
	



