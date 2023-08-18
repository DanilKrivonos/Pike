import numpy as np

def filter_cluster(clt_dat):
    
    #_____________FILTERING_BY_LENGTH_____________________________________
    median = np.median(clt_dat['Length'].values)
    st_div = np.std(clt_dat['Length'].values)
    left = median - 10
    right = median + 10
    filtered_data = clt_dat[clt_dat['Length'] > left]
    filtered_data = filtered_data[filtered_data['Length'] < right]
    #_____________________________________________________________________
    #_______________FILTERING_BY_GC_______________________________________
    median = np.median(clt_dat['GC content'].values)
    st_div = np.std(clt_dat['GC content'].values)
    left = median - 2*st_div
    right = median + 2*st_div
    filtered_data = filtered_data[filtered_data['GC content'] > left]
    filtered_data = filtered_data[filtered_data['GC content'] < right]
    #_____________________________________________________________________
    
    return  filtered_data