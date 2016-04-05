#! /usr/bin/python


###################
# IMPORT PACKAGES #
###################

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from scipy import sparse as sp
import scipy.sparse.linalg as lin
import ConfigParser
import sys

class alluvstrat(object):
  """
  All of the concrete implementation codes go here
  """

  ##############################################
  ## DICTIONARIES: VARIABLES AND THEIR METDATA #
  ##############################################
  
  def lists_of_values(self):
    # Long variable names: what the user sees
    # "self" so is visible outside this function
    self.long_var_names = ['vertical cell size',
                           'horizontal (cross-stream) cell size',
                           'model domain width',
                           'channel width',
                           'channel depth',
                           'in-channel aggradation rate',
                           'channel lateral migration rate',
                           'levee exponential decay constant',
                           'stratigraphic array']
    # Variable type - assuming that we are talking about the whole array and 
    # not its elements
    self.types = ['float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'numpy.ndarray']
    # Variable units: have not yet looked at conventions
    self.units = ['meters', 'meters', 'meters', 'meters', 'meters',
                  'meters / second', 'meters / second', 'meters',
                  'indexed raster: 0=overbank deposit, 1=channel deposit, 2=nothing']
    # Array rank
    self.ranks = [0,0,0,0,0,0,0,0,2]
    # Short variable names: what the program sees
    # "self" so is visible outside
    self.short_var_names = ['dz','dy','B','b','h','etadot_ch','zetadot','etadot_ob_xstar','space']

  def create_static_dictionaries(self):
    self.var_type = dict(zip(self.long_var_names, self.types))
    self.var_units = dict(zip(self.long_var_names, self.units))
    self.var_rank = dict(zip(self.long_var_names, self.ranks))
    self.var_name = dict(zip(self.long_var_names, self.short_var_names))

  def update_var_value_dictionary(self):
    # List of values of variables
    # Useful for getter (and maybe setter, too) and input file parser
    # Maybe use an 'eval' on 'short_var_names' in the future...
    self.var_values_list = [self.dz, self.dy, self.B, self.b, self.h, 
                            self.etadot_ch, self.zetadot, self.etadot_ob_xstar, 
                            self.space]
    self.var_values = dict(zip(self.long_var_names, self.var_values_list))

  def additional_metadata(self):
    self.name = 'AlluvStrat' # Component name

  # EXPERIMENTAL! STORE ALL VARIABLES IN A DICTIONARY
  # MAKE SELF.DICTNAME BE CALLED 'v' THAT POINTS TO SELF:
  # If / when I go this route, should import directly into dict, avoiding need 
  # to construct self.var_values_list (above)
  # Lists can hold different types in different entries
  def var_dict(self):
    # Not going to use this actually, per Eric's advice
    self.var_values_shortnames = dict(zip(self.short_var_names, self.var_values_list))
    v = self.var_values_shortnames # Points at it, so now just have to type, e.g., "v['dy']". Same length as "self.dy"! ('v' for variable)
  

  #######################
  ## INPUT FILE PARSER ##
  #######################
  
  def vars_from_input_file(self, infile):
    config = ConfigParser.ConfigParser()
    config.read(infile)
    # Was thinking about doing something with a dict
    # But for first pass just keep as is and import by hand
    self.dz = config.getfloat('grid', 'dz')
    self.dy = config.getfloat('grid', 'dy')
    self.B = config.getfloat('grid', 'B')
    self.dt = config.getfloat('time', 'dt')
    self.nt = config.getint('time', 'nt')
    self.b = config.getfloat('channel', 'b')
    self.h = config.getfloat('channel', 'h')
    self.eta = config.getfloat('channel', 'eta')    
    self.etadot_ch = config.getfloat('channel', 'etadot_ch')    
    self.zetadot = config.getfloat('channel', 'zetadot')
    self.etadot_ob_xstar = config.getfloat('overbank', 'etadot_ob_xstar')

  def var_checker(self):
    """
    Make sure all the vars are integer divisible into the cells
    """
    print "Variable checker not yet written."
    print "Floor division, potentially with rounding, will occur without warning."

  ######################
  ##  INITIALIZATION  ##
  ######################

 
    # NEED TO CHANGE THESE S.T. WE ALWAYS HAVE INTEGER NCELLS ETC.
    # ALSO, MAJOR CHANGES HAVE OCCURED IN THE PROGRAM
    # SO SHOULD LOOK THROUGH THIS WHOLE THING FOR HOW THE CHANGES 
    # WILL AFFECT IT
  
  def calculate_initial_values(self):
    # Grid
    self.y = np.arange(0,self.B,self.dy)
    self.Bcell = self.B // self.dy # The domain width in cells
    # Channel
    self.b2 = self.b//2 # meters to left and right of centerline
    self.b2cell = self.b2//self.dy # cells to left and right of centerline
    self.eta_cell = self.eta*self.dz # Starting bed 
    # Channel motion
    self.etadot_ob_max = self.etadot_ch # max overbank aggradation rate, defaults to in-channel aggradation rate
    self.agg_counter_ch = 0 # counter for channel going up one cell
    self.etadot_ob = np.zeros(self.B)
    self.agg_ob = np.zeros(self.B) # counts up how long a region has been overbank, and
                         # resets to 0 / adds 1 to grid once a region reaches dz
    self.eta_lastavulse = self.eta
    # Calculation space
    self.space = 2 * np.ones((3*self.h/self.dz, self.Bcell)) # 1-meter cells: 100 m vertical w/ 10 km horizontal
    self.space[:self.h/self.dz,:] = 0 # Fill a channel-depth with overbank deposit
    # Initial channel position
    self.channel_centerline = np.floor((np.random.randint(self.b2+1,self.B-self.b2))) # channel centerline can't be with b2 (1/2-width) of edge; top indexing exclusive, bottom inclusive (so +1)
    self.channel_centerline_subgrid = self.channel_centerline # Channel moves on grid based on its resolution, but this tracks its specific position
    self.channel_left = self.channel_centerline - self.b2//self.dy # Inclusive of channel
    self.channel_right = self.channel_centerline + self.b2//self.dy
    self.space[self.eta//self.dz:(self.eta+self.h)//self.dz,self.channel_left:self.channel_right+1] = 1 # Initial channel; "+1" b/c of exclusive indexing
    self.channel_loc = np.zeros(self.B)
    self.channel_loc[self.channel_left:self.channel_right+1] = 1 # +1 to be inclusive
    self.nochannel_loc = 1-self.channel_loc
    # Time steps and friends
    self.t = 0 # Total time counter
    self.grid_update_check_time = self.h / (self.etadot_ch)
    self.npast = 0 # For adjusting the domain
  
  ########################
  ##  CHANNEL DYNAMICS  ##
  ########################

  def lateral_migration(self,edge_condition='wall'):
    # Channel lateral migration: random walk at pace set by
    # migration rate
    self.channel_centerline_subgrid += (2*np.random.randint(0,2) - 1) * self.zetadot * self.dt # 1 or -1 times distance
    # Only update channel position if it has moved
    if self.channel_centerline != np.round(self.channel_centerline_subgrid):
      # Update channel to new subgrid position
      self.channel_centerline = np.round(self.channel_centerline_subgrid)
      # Then handle edges if necessary
      if edge_condition == 'wall':
      # Solid wall at edges
      # Don't let channel move outside of domain
        if self.channel_centerline > (self.Bcell - self.b2cell):
          self.channel_centerline = (self.Bcell - self.b2cell)
        if self.channel_centerline < self.b2cell:
          self.channel_centerline = self.b2cell # no "+1" needed b/c of 0 indexing
        self.channel_left = self.channel_centerline - self.b2cell # Inclusive of channel
        self.channel_right = self.channel_centerline + self.b2cell # Inclusive of channel
      elif edge_condition == 'wraparound':
        print "Placeholder to send one edge to the other and create an infite"
        print "space but with a given width over which things can happen"
        print "Need to make edits in other parts of the code (calculation_space,"
        print "at least) to win"
        print "NOT YET IMPLEMENTED"
        sys.exit()
      else:
        print "Select a valid edge condition: wall or wraparound"
        sys.exit()
      self.space[self.eta//self.dz:(self.eta+self.h)//self.dz,self.channel_left:self.channel_right+1] = 1 # Automatically floors everything for indexing, so could just have eta/dz

  def aggradation_avulsion(self):

    # Channel vertical aggradation
    self.eta += self.etadot_ch*self.dt # Always floor divide, so this is OK

    # Channel avulsion
    if (self.eta - self.eta_lastavulse) >= self.h:
      # Start a grid of possible avulsion locations
      self.avulsion_locs = np.zeros(self.B)
      # Not too close to edges
      self.avulsion_locs[self.b2+1:self.B-self.b2] = 1
      # Get channel location an extra time
      self.channel_loc_extended = np.zeros(self.B)
      self.channel_loc_extended[max(0,self.channel_left-self.b2cell-1):self.channel_right+self.b2cell+2] = 1 # +1 to be inclusive, indexing past end is not problem (will just go to end)
      # Remove present channel position from possible avulsion locations 
      self.avulsion_locs[self.channel_loc_extended.nonzero()] = 0 # Remove current channel position
      # At this point, all points at which avulsion_locs == 1 are potential post-avulsion channel locations
      self.future_channel_locs = ((self.avulsion_locs==1).nonzero())[0]
      
      print "Channel centerline at:", self.channel_centerline
      
      # Find min elevation
      self.min_elev_grid = ((self.space[:,self.future_channel_locs])==2).nonzero()[0].min() # nonzero()[0] for rows: want lowest row
      # Find surface elvation in all columns
      self.elev = -1 * np.ones(self.B) # assuming we'll never have negative elevations
      for i in range(self.future_channel_locs.shape[0]):
        self.elev[self.future_channel_locs[i]] = (self.space[:,self.future_channel_locs[i]]==2).nonzero()[0].min() # Find lowest elevation occupied by nothing   
        self.min_elev_cells = (self.elev == self.min_elev_grid).nonzero()
      # Pick new channel centerline
      # Randomly choose if multiple sites have the minimum elevation; lots of [0] to get array out of tuple
      self.channel_centerline = self.min_elev_cells[0][np.random.randint(0,self.min_elev_cells[0].shape[0])] # min_elev_cells[0] b/c I think it should be a vector
      self.channel_centerline_subgrid = self.channel_centerline # Update the subgrid centerline based on the avulsed-to centerline
      
      # Once I get a good centerline, construct the rest of the channel
      self.channel_left = self.channel_centerline - self.b2//self.dy # Inclusive of channel
      self.channel_right = self.channel_centerline + self.b2//self.dy
      self.channel_loc = np.zeros(self.B)
      self.channel_loc[self.channel_left:self.channel_right+1] = 1 # +1 to be inclusive
      # Make the channel lie at the elevation of the lowest cell with deposits
      self.eta = ((self.space[:,self.channel_loc.nonzero()])==2).nonzero()[0].min()*self.dz - self.h # nonzero()[0] for rows: want lowest row
      self.space[self.eta//self.dz:(self.eta+self.h)//self.dz,self.channel_left:self.channel_right+1] = 1 # Carve channel
      self.space[((self.eta+self.h)//self.dz):,self.channel_left:self.channel_right+1] = 2 # Incise: remove all material above avulsed channel's new location
      # Reset eta_lastavulse
      self.eta_lastavulse = self.eta

    # Overbank aggradation
    self.channel_loc = np.zeros(self.B)
    self.channel_loc[self.channel_left:self.channel_right+1] = 1 # +1 to be inclusive
    self.nochannel_loc = 1 - self.channel_loc
    
    self.etadot_ob[:self.channel_left] = self.etadot_ob_max * np.exp(-(self.channel_left-self.y[:self.channel_left])/self.etadot_ob_xstar)
    self.etadot_ob[self.channel_right:] = self.etadot_ob_max * np.exp(-(self.y[self.channel_right:]-self.channel_right)/self.etadot_ob_xstar)

    self.agg_ob += self.nochannel_loc * self.etadot_ob*self.dt # Removes all in channel; important b/c I can still have aggradation on right end if channel pinned on R side
    self.agg_ob[self.channel_loc.nonzero()] = 0 # remove all progress where channel passes
    self.agg_cols = (self.agg_ob >= self.dz).nonzero()[0] # Find cols that will aggrade on main grid in this time-step
    for i in range(self.agg_cols.shape[0]):
      self.space[(self.space[:,self.agg_cols[i]]==2).nonzero()[0].min(),self.agg_cols[i]] = 0 # new overbank cell
      self.agg_ob[self.agg_cols[i]] = 0

  def update_t(self):
    self.t += self.dt

  def raise_the_roof(self):
    """
    Increases domain size if needed.
    The problem is solved upside-down because [0,0] is upper left of array
    So vstack to bottom
    """
    # If we have just past another time check (npast and the floor divide keep 
    # this a binary outcome)
    if self.t//self.grid_update_check_time - self.npast == 1:
      self.npast += 1 # Update counter
      # Then check if anything is in striking distance (i.e. within h) of the 
      # top of the domain
      # "+1" for safety! Should floor divide and be OK, but don't want to 
      # take the chance - maybe not actually! need the 1.
      # If any cells are not empty
      if (self.space[(self.space.shape[0] - (self.h//self.dz) - 1):,:] != 2).any():
        # Then add another channel depth to the top of the array
        self.space = np.vstack((self.space,2*np.ones(((self.h//self.dz,self.space.shape[1])))))
        #print "Raising the roof"


class structure(alluvstrat):
  """
  The abstract structure of the program that talks to the concrete 
  implementation goes here.
  """


  def __init__(self):
    super(structure, self).lists_of_values()

  def initialize(self):
    print "Initializing class alluvstrat"
    super(structure, self).vars_from_input_file('input')
    super(structure, self).lists_of_values()
    super(structure, self).create_static_dictionaries()
    super(structure, self).calculate_initial_values()
    super(structure, self).additional_metadata()
    super(structure, self).var_checker()
    super(structure, self).calculate_initial_values()
    """
    self.lists_of_values()
    self.create_static_dictionaries()
    self.calculate_initial_values()
    self.additional_metadata()
    self.vars_from_input_file()
    self.var_checker()
    self.calculate_initial_values()
    """
    print "class alluvstrat initialized"
  
  def update(self):
    self.update_t()
    self.lateral_migration()
    self.aggradation_avulsion()
    self.raise_the_roof()

  def finalize(self, rmv_start=True):
    print "Finalizing."
    if rmv_start:
      # Remove first channel-depth of rows: this is affected by initial condition
      self.space = self.space[self.h/self.dz:,:]
    

