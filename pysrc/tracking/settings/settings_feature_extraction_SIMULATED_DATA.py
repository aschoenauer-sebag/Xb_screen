result_folder = '../resultData/simulated_traj/simres'
outputFile = 'dist_sim_{}.pkl'
ctrl_exp_filename = 'ctrl_exp_{}.pkl'
data_folder = '../resultData/simulated_traj/simres/plates'
##filename for trajectory features
#filename='hist_tabFeatures_{:>03}.pkl'

##filename format for the wells
#format_='{:>03}'

mitocheck_file = None
quality_control_file = '/cbio/donnees/aschoenauer/workspace2/Xb_screen/data/qc_simulated.txt'

#Option to say that we're only dealing with numerical features for the moment
histDataAsWell=False

#Option to say if we want to redo experiments anyway
redo=True
#to tell it's not the xb screen
xb_screen=False