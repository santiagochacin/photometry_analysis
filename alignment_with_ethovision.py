import pandas as pd
import sys
import numpy as np
sys.path.append('/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/photometry_analysis/photometry_3.0_module/')
import photometry_3



##alignment steps

# create dataframe of TTL sync pulses with a time column from start of the video

# create dataframe of grooming/behavior information from ethovision file, also with a time column from start of the video

# get the indicies of sync pulses from the abf recording

# form continous array with the intermediate frame numbers between pulses and time stamps to add into the photometry df

# pad the array with zeros at the end of the sweep

# take the photometry data array, and, for each sample point:
#   1) find the time of the closest frame analyzed by ethovision in the ethovision array
#	2) get the time of the frame
#   3) add the grooming/behavior value for this frame to the photometry data
# by iterating over all the samples in the photometry data file, this will "upsample" the behavior analysis and align it 
# to the information already in the photometry data for the abf file

# fold this back into the original photometry df 


def read_ethovision_dataframe(path_to_file, rows_to_skip, list_of_behavior_categories):
	ethovision_df = pd.read_excel(path_to_file, skiprows=rows_to_skip) 

	behavior_dict = {}
	
	behavior_dict['time'] = ethovision_df['Trial time'].values[1:]

	for behavior in list_of_behavior_categories:

		behavior_array = []

		for i in range(1, len(ethovision_df[behavior].values)):
			if ethovision_df[behavior].values[i] == '-':
				behavior_array.append(np.NaN)

			elif (ethovision_df[behavior].values[i] == 0) or (ethovision_df[behavior].values[i] == 1):
				behavior_array.append(int(ethovision_df[behavior].values[i]))

			else:
				behavior_array.append(float(ethovision_df[behavior].values[i]))

		behavior_dict[behavior] = behavior_array

	behavior_df = pd.DataFrame(behavior_dict)

	return(behavior_df)

def read_sync_dataframe(path_to_file):
	
	syn_pulses_df = pd.read_csv(path_to_file, index_col=False, usecols=[0,1], names=['computer time', 'frame_num'], skiprows=1)
	if syn_pulses_df.loc[0]['computer time'] == 500:
		syn_pulses_df = pd.read_csv(path_to_file, names = ['frame num', 'computer time', 'frame_num'])
	start_time = syn_pulses_df['computer time'][0]
	interval = syn_pulses_df['computer time'][1] - syn_pulses_df['computer time'][0]
	video_time = []
	video_time.append(interval)
	for i in range(1, len(syn_pulses_df['computer time'].values)):
		video_time.append(syn_pulses_df['computer time'][i]-start_time+interval)

	syn_pulses_df['video_time'] = pd.Series(video_time, index = syn_pulses_df.index )

	return(syn_pulses_df)

def interpolate_array(indicies_to_match, reference_array, ref_index_start, start_at_1=True):
	#use start at 1 in reference to trim initial time samples--alows for precise alignment over multiple sweeps
	#ref indicies e.g sweep 001 would be (0, 0), sweep 002 would be (1*len(indicies_to_match), 1*len(indicies_to_match) + len(indicies_to_match))
	output_array = np.array([])

	for index in range(0, len(indicies_to_match)-1):
		to_fill = indicies_to_match[index+1] - indicies_to_match[index]
		#print('size to fill ' + str(to_fill))
		#print('reference array start index' + str(index+ref_index_start))
		#print('time at start' + str(reference_array[index+ref_index_start]))
		#print('reference array end index' + str(index+ref_index_start+1))
		#print('time at end' + str(reference_array[index+ref_index_start+1]))
		output_array = np.concatenate((output_array, np.linspace(reference_array[index+ref_index_start], 
																reference_array[index+ref_index_start+1], to_fill)))

	return(output_array)

def push_nearest_index(down_sampled_array, up_sampled_array, secondary_array=[]):

	new_upsampled = []

	point = 0
	while point < len(up_sampled_array):

		diff = float('inf')
		to_push = np.nan

		if point % 10000 == 0:
			print(point)

		if up_sampled_array[point] == np.nan:
			new_upsampled.append(np.nan)

		else:
		#a recursive call here may make this yet faster
			point_downsampled = 0

			while point_downsampled < len(down_sampled_array):
				if (down_sampled_array[point_downsampled] - up_sampled_array[point] >= 0) and (down_sampled_array[point_downsampled] - up_sampled_array[point] < diff):
					diff = down_sampled_array[point_downsampled] - up_sampled_array[point]
					to_push = down_sampled_array[point_downsampled]

					if point_downsampled+1 < len(down_sampled_array):
						if down_sampled_array[point_downsampled+1] - up_sampled_array[point] > diff:
							break
					
				point_downsampled += 1

		#if len(secondary_array) != 0:
		#	new_upsampled.append(secondary_array[point_downsampled])
		#else:
		#	new_upsampled.append(to_push)

		point += 1

	return(np.array(new_upsampled))


def push_nearest(down_sampled_array, up_sampled_array, secondary_array=[]):

	new_upsampled = []

	point = 0

	down_sampled_point = 0

	while point < len(up_sampled_array):

		#if point % 10000 == 0:
			#print(point)

		to_push, diff, down_sampled_point = inner_nearest(up_sampled_array, point, down_sampled_array, down_sampled_point, secondary_array)

		new_upsampled.append(to_push)

		point += 1

	return(new_upsampled)



def inner_nearest(up_sampled_array, up_sampled_point, down_sampled_array, down_sampled_start=0, secondary_array=[]):

	diff = float('inf')
	to_push = 'NaN'

	point_downsampled = down_sampled_start

	while point_downsampled < len(down_sampled_array):

		if (down_sampled_array[point_downsampled] - up_sampled_array[up_sampled_point]) > 0 and (down_sampled_array[point_downsampled] - up_sampled_array[up_sampled_point]) < diff :

			diff = down_sampled_array[point_downsampled] - up_sampled_array[up_sampled_point]
			
			if len(secondary_array) != 0:
				to_push = secondary_array[point_downsampled]
			else:
				to_push = down_sampled_array[point_downsampled]

			if point_downsampled < len(down_sampled_array) - 1:
				if (down_sampled_array[point_downsampled+1] - up_sampled_array[up_sampled_point]) > (down_sampled_array[point_downsampled] - up_sampled_array[up_sampled_point]):
					break

		point_downsampled += 1

	return(to_push, diff, point_downsampled)



def loop_over_photometry_sweeps_and_align(photometry_df, ethovision_df, behavior_list, TTL_pulse_df):
	sweeps_list = sorted(list(set(photometry_df.index.get_level_values(0))))
	dfs_with_behavior_data = []

	
	for sweep in range(len(sweeps_list)):
		#get indicies of sync pulses for each sweep
		sweep_sync_indicies = photometry_3.count_events_in_array(photometry_df.loc[sweeps_list[sweep]]['encoder_windowed'].values,
								photometry_df.loc[sweeps_list[sweep]]['time'][1]-photometry_df.loc[sweeps_list[sweep]]['time'][0],
								1, threshold=3, up=True)[1]
		if len(sweep_sync_indicies)<2:
			break
		else:

			frames_array = interpolate_array(sweep_sync_indicies, TTL_pulse_df['frame_num'].values, sweep*len(sweep_sync_indicies))
			time_array = interpolate_array(sweep_sync_indicies, TTL_pulse_df['video_time'].values, sweep*len(sweep_sync_indicies))

			print('upsampling behavior data for sweep' + str(sweep))

			sweep_df_truncated = photometry_df.loc[sweeps_list[sweep]][sweep_sync_indicies[0]:sweep_sync_indicies[-1:][0]]

			sweep_df_truncated['video_times_from_sync'] = time_array
			sweep_df_truncated['frame'] = frames_array

			for behavior in behavior_list:
				upsampled_with_behavior = push_nearest(ethovision_df['time'].values, 
													sweep_df_truncated['video_times_from_sync'].values, 
													ethovision_df[behavior].values)

				sweep_df_truncated[behavior] = upsampled_with_behavior

			dfs_with_behavior_data.append(sweep_df_truncated)

	behavior_df = pd.concat(dfs_with_behavior_data,keys=sweeps_list, names=['sweeps', 'index'])

	return(behavior_df)
















































































