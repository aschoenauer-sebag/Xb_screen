import os


def check_h5_files_exist(rawDataFolder = '/share/data20T/mitocheck/compressed_data', 
                         h5Folder = '/share/data40T/aschoenauer/drug_screen/results_August_2016/mito_joint_classifier'):
    
    plates = filter(lambda x: 'Valid' not in x, os.listdir(rawDataFolder))
    missing_h5=[]
    past_len=0
    
    for plate in plates:
        well_h5files = ['00{}_01'.format(el.split('--')[0]) for el in os.listdir(os.path.join(rawDataFolder, plate))]
        try:
            existing_h5files=set(os.listdir(os.path.join(h5Folder, plate, 'hdf5')))
        except:
            missing_h5.extend([(plate, el) for el in well_h5files])
        else:
            missing_h5.extend([(plate, el) for el in well_h5files if '{}.ch5'.format(el) not in existing_h5files])
        if len(well_h5files)>past_len:
            print plate, ' missing ', len(well_h5files)-past_len
        past_len = len(well_h5files)
        
    return missing_h5 