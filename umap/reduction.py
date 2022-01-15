import umap
from sklearn.datasets import load_digits
from scipy.sparse.csgraph import connected_components
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.manifold import TSNE
import time

def main():
    digits = load_digits()
    digits.target = [float(digits.target[i]) for i in range(len(digits.target))]

    # UMAP
    start_time = time.time()
    embedding = umap.UMAP().fit_transform(digits.data)
    interval = time.time() - start_time
    plt.scatter(embedding[:,0],embedding[:,1],c=digits.target,cmap=cm.tab10)
    plt.colorbar()
    plt.savefig('umap.png')

    # t-SNE
    plt.clf()
    start_time2 = time.time()
    tsne_model = TSNE(n_components=2)
    tsne = tsne_model.fit_transform(digits.data)
    interval2 = time.time() - start_time2
    plt.scatter(tsne[:,0],tsne[:,1],c=digits.target,cmap=cm.tab10)
    plt.colorbar()
    plt.savefig('tsne.png')

    print('umap : {}s'.format(interval))
    print('tsne : {}s'.format(interval2))

if __name__ == "__main__":
    main()
