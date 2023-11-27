# This version reads the new json file format.

import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from sklearn import metrics
import json

from sklearn.preprocessing import StandardScaler


def asse(a, b):
    val = b.groupby(a)
    mean = b.groupby(a).mean()

    j = 0
    sum = 0
    l = 0
    for a, b in val:
        temp = b.tolist()
        ln = len(temp)
        l += ln
        for i in range(ln):
            sum += ((temp[i] - mean[j]) ** 2)
        j += 1
    print(sum / l)
    return (sum / l)

def create_array_for_childs(data, dataset):
    have_children = True
    children = None
    if "children" not in data.keys():
        have_children = False
    else:
        children = data["children"]
        if children is None or children == []:
            have_children = False

    if have_children == False:
        name = data["name"]
        # Get all values
        data_list = str(name).strip().split("_")
        # print(data_list)
        if (len(data_list) > 0):
            profile = []
            for value in data_list:
                if "=" in value:
                    key_val = value.strip().split("=")
                    val = float(key_val[1])

                    # if val == 0.0:
                        # continue
                    profile.append(val)
                else:
                    profile.append(float(value))
            # Append the profile to the final dataset
            if len(profile) == 23:
                dataset.append(profile)
    else:
        if children != None:
            for child in children:
                # Recursion
                create_array_for_childs(child, dataset)


def filter_data(dims, dataset):
    new_data = []

    for d in dataset:
        new_d = []
        for dim in dims:
            new_d.append(d[dim])
        new_data.append(new_d)

    return new_data

# dataset = np.delete(dataset, 6, axis=1)

# X_embedded = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=3).fit_transform(arr)

def Reduce_Dimensionality(dataset, n = 0, perp = 0):
    if n != 0:
        tsne = TSNE(n_components=n, learning_rate='auto', init='random', method='exact',
                    random_state=1, perplexity=perp, n_jobs=-1)
    else:
        tsne = TSNE(n_jobs=-1, random_state=1)
    return tsne.fit_transform(dataset)


def Optimize_DBSCAN(reduced_data):
    s_max = -np.inf
    d_e = -1
    d_i = -1
    min_count = len(dataset)
    n_clusters = 0
    # reduced_data = np.array(dataset)
    # reduced_data = np.delete(reduced_data, 6, 1)
    # print(reduced_data[0][0:10])
    df = pd.DataFrame(reduced_data)
    # ss = MinMaxScaler(feature_range=(-1, 1))
    ss = StandardScaler()
    # print(reduced_data[0][0:10])
    X = ss.fit_transform(reduced_data)
    # print(X[0][0:10])
    max_k = int(np.sqrt(len(X)))
    min_k = int(max_k/2)
    print("max_k", max_k)
    for eps in np.linspace(0.01, 1, 100):
        for k in range(1, max_k, 1):
            # print(i)
            # for xi in np.linspace(1, 100, 100):
            try:
                # print("Here", eps, i)
                clustering = DBSCAN(eps=eps, min_samples=k, n_jobs=-1).fit(X)
                # clustering = MeanShift(min_bin_freq=eps, n_jobs=-1, max_iter=i).fit(reduced_data)
                # clustering = OPTICS(eps=eps, min_samples=i, xi=xi, n_jobs=-1).fit(reduced_data)
                # clustering = AffinityPropagation(damping=eps, max_iter=i).fit(reduced_data)

                cluster_labels = clustering.labels_
                df['label'] = cluster_labels

                c = 0
                for j in clustering.labels_:
                    if j == -1:
                        c += 1

                cluster_labels = clustering.labels_
                unique, counts = np.unique(cluster_labels, return_counts=True)
                ele = dict(zip(unique, counts))

                l = len(unique)
                try:
                    count_one = ele[-1]
                except:
                    count_one = len(dataset)

                if l > 5 and l < 15:
                    try:
                        s = metrics.silhouette_score(df, df['label'])
                    except:
                        continue

                    if count_one < min_count and s > s_max and s > 0:
                        d_e = eps
                        d_i = k
                        s_max = s
                        min_count = count_one
                        n_clusters = len(unique)

                        print(d_e, d_i, s_max, min_count, n_clusters)
            except:
                continue

    return (d_e, d_i, s_max, min_count, n_clusters)

def perform_clustering(dataset, orig_dataset, fname):
    # ss = MinMaxScaler(feature_range=(-1, 1))
    ss = StandardScaler()
    X = ss.fit_transform(dataset)

    reduced_data = Reduce_Dimensionality(X)
    eps, n_samples, score, n_outliers, n_clusters = Optimize_DBSCAN(reduced_data)

    # ss = MinMaxScaler(feature_range=(-1, 1))
    ss = StandardScaler()
    X = ss.fit_transform(reduced_data)

    clustering = DBSCAN(eps=eps, min_samples=n_samples, n_jobs=-1).fit(X)
    # print(clustering)
    df = pd.DataFrame(reduced_data)
    cluster_labels = clustering.labels_
    df['label'] = cluster_labels

    metrics.silhouette_score(df, df['label'])

    labels = clustering.labels_
    # print(labels)
    # n_clusters_ = 6
    core_samples_mask = np.zeros_like(labels, dtype=bool)
    X = reduced_data

    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]

        class_member_mask = labels == k

        xy = X[class_member_mask & core_samples_mask]
        plt.plot(
            xy[:, 0],
            xy[:, 1],
            "o",
            markerfacecolor=tuple(col),
            markeredgecolor="k",
            markersize=18,
        )

        xy = X[class_member_mask & ~core_samples_mask]
        plt.plot(

            xy[:, 0],
            xy[:, 1],
            "o",
            markerfacecolor=tuple(col),
            markeredgecolor="k",
            markersize=18,
        )

    orig_dataset = np.array(orig_dataset)
    fname = fname + "_clusters"
    labels = np.array(labels)
    labels = labels[..., np.newaxis]
    print(orig_dataset.shape, labels.shape)
    results = np.concatenate((X,  np.array(labels)), axis=1)
    results = np.concatenate((dataset, results), axis=1)
    results = np.concatenate((orig_dataset,  results), axis=1)
    data_file_path = "../channel1_data/datasets2/" + fname + ".txt"
    np.savetxt(fname=data_file_path, X=results)
    # plt.title("Number of clusters: %d" % )
    plt.show()


if __name__ == "__main__":
    fname = 'channel_I256_J192_K256_Full'
    ext = '.json'
    filepath = "../channel1_data/datasets2/" + fname + ext
    file = open(filepath)
    data = json.load(file)
    dataset = []
    create_array_for_childs(data, dataset)

    # Best results with these dims
    # dims_to_keep = [0, 5, 6, 7]
    red_dataset = np.asarray(dataset)[:, 4:]
    # dims_to_keep = [0, 5, 13, 14, 15, 16, 17]     Couette
    # dims_to_keep = [0, 13, 14, 15, 16, 17]        Cylinder2
    # dims_to_keep = [0, 6, 14, 16]                 Bernard
    # dims_to_keep = [0, 5, 13, 15]                 Bernard 2
    # dims_to_keep = [0, 12, 14, 15, 16]         §     Fort
    # dims_to_keep = [2, 13, 14, 15, 16]
    dims_to_keep = [5, 6, 12, 13, 14, 15, 16, 17]
    indices = np.asarray(dataset)[:, 0:3]
    new_dataset = filter_data(dims=dims_to_keep, dataset=red_dataset)
    new_dataset = np.array(new_dataset)
    print(new_dataset[0])
    perform_clustering(new_dataset, indices, fname)