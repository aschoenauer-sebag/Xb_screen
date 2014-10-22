import simulator

from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-w", "--what", dest="what",
                      help="what to do: either nb_traj or a feature (check feature list)")
    parser.add_option("-s", "--stat", dest="stat",
                      help="stat: either avg or full")
    
    (options, args) = parser.parse_args()
    
    rs = simulator.RealData()
    rs(options.what, options.stat)
    
    