'''
Load the Science News data and save it as a .mat file for use with MATLAB
'''



import numpy as np
import scipy.io as sio



def load_SN():
    '''
    Load the Science News data and save it as a .mat file for use with MATLAB
    '''
    
    with np.load('./ScienceNews.npz') as SN:
        cmatrix = SN['matrix']
        articles = SN['doc_titles']
        words = np.array([x[0][0] for x in SN['words']])
        classes = np.array(SN['doc_class'][0,:])
        class_names = [x[0] for x in SN['score_titles'][0]]

    sio.savemat('../ScienceNews.mat', {'cmatrix':cmatrix, 'articles':articles, 'words':words, 'classes':classes, 'class_names':class_names})

    return



if __name__ == '__main__':
    load_SN()