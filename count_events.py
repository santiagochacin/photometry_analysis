def count_events_in_array(trace, dt, time_to_wait, threshold=0, up=True):
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
	#samples_to_skip = int(time_to_wait/dt)
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
		if comp(selection[i],threshold):
			EventCounter +=1
			sample_points.append(i)	
			while i<len(selection) and comp(selection[i],threshold):
				i+=1 # skip values if index in bounds AND until the value is below/above threshold again
		else:
			i+=1
	
	time_points = [sample_point*dt for sample_point in sample_points];
	return (EventCounter, sample_points, time_points)


