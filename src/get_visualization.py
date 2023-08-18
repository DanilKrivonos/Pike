import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

def get_color(tax_list):
    
    color = ''
    
    while color not in tax_list.values() and color == '':
        
        color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    
    return color

def get_visualisation(output, 
                      barcode, 
                      filtered_add):

    fig, axs = plt.subplots(figsize=(10, 10))
    colordict = {}
    
    for clt in filtered_add.keys():
        
        mean_dists = np.median(squareform(pdist(filtered_add[clt][['1 UMAP COMPONENT', '2 UMAP COMPONENT']].values), 'braucyrtis'), axis=1)
        ref_idx = filtered_add[clt][mean_dists == np.min(mean_dists)].index[0]

        #____Looking for cluster centroids for figure______
        x =  filtered_add[clt].loc[ref_idx]['1 UMAP COMPONENT']
        y = filtered_add[clt].loc[ref_idx]['2 UMAP COMPONENT']
        colordict[clt] = get_color(colordict)
        
        plt.text(x, y, 
                 f'Cluster {clt}', 
                 fontsize=3, 
                 fontweight='book')
        sns.scatterplot(data=filtered_add[clt], 
                        x='1 UMAP COMPONENT', 
                        y='2 UMAP COMPONENT', 
                        s=1,
                        alpha=.3,
                        color=colordict[clt],
                        ax=axs)
        sns.scatterplot(data=filtered_add[clt], 
                        x='1 UMAP COMPONENT', 
                        y='2 UMAP COMPONENT', 
                        s=3,
                        alpha=.7,
                        color=colordict[clt],
                        ax=axs)
    #sns.despine(offset=10, trim=True)
    plt.xlabel('1 UMAP COMPONENT', fontweight='book')
    plt.ylabel('2 UMAP COMPONENT', fontweight='book')
    plt.grid()
    plt.legend([],[], frameon=False)
    plt.savefig(f'{output}/{barcode}/reults/{barcode}.pdf')
    plt.savefig(f'{output}/{barcode}/reults/{barcode}.pdf', dpi=800)