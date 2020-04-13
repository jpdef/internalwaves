# Configuration Paramaters

The configuration parameters written in JSON format describe the simluation of an internal wave field in a cube defined by the range and depth paramaters. The simulation runs as long as the time paramters specify. For efficient storage the user can specify spacial points to sample the cube, instead of all points. The simulation is comprised of a set of waves. As input N frequencies and M modes are provided and NxM waves are produced. Each of those waves has a set of amplitudes and directions (headings). The output is specified by the path and the ftype. 

   * range_end  : the end range in meters
   * range_res  : the increments of range in meters
   * depth_end  : the bottom depth in meters
   * depth_res  : the increaments of depth in meters
   * time_stop  : the time to stop the simulation in seconds
   * time_step  : the increment of time in seconds
   * x_samples  : longitudunal positions to sample the waves
   * y_samples  : latitudunal positions to sample the waves
   * z_samples  : depth positions to sample the waves
   * freqs      : the wave frequencies
   * modes      : the wave mode numbers
   * amps_real  : real part of the wave amplitude
   * amps_imag  : imaginary part of the wave amplitude
   * headings   : the direction in degrees of the wave
   * path       : directory where data output goes
   * ftype      : type of file
       * 0 - feather file (binary datafile)
       * 2 - mp4 movie file
