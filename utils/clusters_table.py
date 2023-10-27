import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

nodesdf = pd.read_csv("data/results/networks/mmseqs_27_09_22/data/merged_net_nodes_table.csv", usecols=['name', 'classe', 'family', 'genus', 'order', 'phylum', 'species', 'superkingdom'])
edgesdf = pd.read_csv("data/results/networks/mmseqs_27_09_22/data/merged_net_edge_table.csv", usecols=['cluster','name'])
edgesdf[['key','target']] = edgesdf.name.str.split(r'\s\(.+\)\s', expand=True)
edgesdf = edgesdf.drop(columns=['name'])

# add key row to target column
rows = pd.DataFrame(
    [[edgesdf[edgesdf.key == k]['cluster'].iloc[0], k, k] for k in edgesdf.key.unique()],
    columns=['cluster','key','target']
    )
edgesdf = pd.concat([edgesdf, rows], axis=0, ignore_index=True)


def diversity(row):
    return ",".join(list(set(row.to_list())))
def n_diversity(row):
    return len(list(set(row.to_list())))

def network(df):
    """Returns a the cluster/component number as a new column, given a dataframe in with 'key'
    and 'target' columns."""
    G = nx.from_pandas_edgelist(df, source='key', target='target')
    components = list(nx.connected_components(G))
    cluster_num = {}
    main_nodes = {}
    for i, component in enumerate(components):
        for node in component:
            cluster_num[node] = i+1
        subgraph = G.subgraph(component)
        degrees = dict(subgraph.degree())
        main_node = max(degrees, key=degrees.get)
        main_nodes[i+1] = main_node
    df['cluster_network'] = df['key'].map(cluster_num)
    df['main_node'] = df['cluster_network'].map(main_nodes)
    return df

def generate_colors(labels: list) -> dict:
    """Given a list of labels, generates len(list) distinct colors."""
    # Define a colormap
    cmap = plt.cm.get_cmap('Set3')

    # Map each string to a unique integer using the hash function
    hashes = [hash(label) for label in labels]

    # Normalize the hashes to [0, 1]
    norm = plt.Normalize(min(hashes), max(hashes))

    # Generate a list of colors based on the normalized hashes using the colormap
    colors = [cmap(norm(hashval)) for hashval in hashes]

    # Print the list of colors
    return dict(zip(labels, colors))
    
def draw_net(df, cluster):
    sk_colors = {
        'eukaryota':'green',
        'bacteria':'red'
    }
    netdf_filtered = df[df.cluster_network == cluster]
    G = nx.from_pandas_edgelist(netdf_filtered, source='key', target='target')
    G.remove_edges_from(nx.selfloop_edges(G))
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=30, alpha=0.8, node_color=[sk_colors.get(nodesdf[nodesdf.name == n]['superkingdom'].values[0], 'blue') for n in G.nodes()])
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    nx.draw_networkx_labels(G, pos, labels={n: n for n in G.nodes()}, font_size=5)
    plt.axis('off')
    plt.savefig("data/results/networks/mmseqs_27_09_22/test_cluster.png", dpi=450)

def main():
    netdf = network(edgesdf)

    # How many nodes in cluster
    nclusters_table = edgesdf.groupby('cluster_network').agg(counts=pd.NamedAgg(column="target", aggfunc='count'))
    print(nclusters_table)

    # Diversity
    mergedf = nodesdf.merge(edgesdf, left_on='name', right_on='target')
    divclusters_table = (
        mergedf.groupby('cluster_network').agg(
            diversity=pd.NamedAgg(column="superkingdom", aggfunc=diversity),
            n_sk=pd.NamedAgg(column="superkingdom", aggfunc=n_diversity))
        .sort_values(by='n_sk',ascending=False)
    )
    print(divclusters_table)
    
    draw_net(netdf, cluster=62)

    # full table
    clusters_table = pd.concat([divclusters_table, nclusters_table, netdf.drop_duplicates(subset='cluster_network').set_index('cluster_network')['main_node']], axis=1).sort_values('counts', ascending=False)
    clusters_table.to_csv("data/results/networks/mmseqs_27_09_22/clusters_information_table.csv")

    
if __name__ == "__main__":
    main()