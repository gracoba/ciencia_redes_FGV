# Write a function that identifies all nodes in a triangle relationship with a given node.
def node_triplets(G, n):
    """
    Returns the triplets in a graph `G` that are involved in a triangle relationship with the node `n`.
    """
    triangle_nodes = []

    # Iterate over all possible triangle relationship combinations
    for n1, n2 in combinations(G.neighbors(n), 2):
        triangle_nodes.append([n1,n2])
    return triangle_nodes
# ==============================================================================
def strengh_barret(G, node = None, weight = 'weight', sentido = 'out', vizinhos = None):
  nos = list(G.nodes)
  
  a_weight = nx.to_numpy_matrix(G,weight = weight)
  a_matrix = nx.to_numpy_matrix(G,weight = vizinhos)
  
  # AB = np.array([[0, 0, 1, 0],
  #               [1, 0, 1, 0],
  #               [1, 0, 0, 0],
  #               [0, 1, 0, 0],])
  # print('k_2 out = 1:', np.sum(AB,axis = 0))
  # print('k_2 in  = 2:', np.sum(AB,axis = 1))

  # nos saindo
  s_i = np.sum(a_weight, axis = 0)  
  k_i = np.sum(a_matrix, axis = 0)  
  
  # nos chegando
  s_j = np.sum(a_weight, axis = 1)  
  k_j = np.sum(a_matrix, axis = 1) 
  
  if node == None:
    if sentido == sentido:
      return k_i, s_i
    else:
      return k_j, s_j
  else: 
    index = nos.index(node)
    
    if sentido == 'out':
      return k_i[0,index], s_i[0,index]
    else:
      return k_j[index,0], s_j[index,0]

# ==============================================================================    
def Cw_barret(G,n, weight = 'weight', sentido = 'out', vizinhos = None):
# Barrat clustering coefficient
# https://www.pnas.org/doi/epdf/10.1073/pnas.0400087101
# https://toreopsahl.com/tnet/weighted-networks/
# http://jponnela.com/web_documents/a9.pdf

  tmp0 = nx.get_edge_attributes(G,weight)
  tmp1 = len(list(tmp0.keys())[0])
   
  # inicia variavel 
  Cw_barret = 0

  # executa em uma direcao dos nos do triplet
  for a in node_triplets(G, n):
    if tmp1 < 3:
      w_ij = G[n][a[0]][weight]
      w_ih = G[n][a[1]][weight]
    else:
      w_ij = G[n][a[0]][0][weight]
      w_ih = G[n][a[1]][0][weight]

    a_ij = G.has_edge(n, a[0])
    a_ih = G.has_edge(n, a[1])
    a_hj = G.has_edge(a[0], a[1])
    Cw_barret = Cw_barret + ((w_ij + w_ih)/2)*(a_ij * a_ih * a_hj)
  
  # executa na outr direcao dos nos do triplet
  for a in node_triplets(G, n):
    if tmp1 < 3:
      w_ij = G[n][a[0]][weight]
      w_ih = G[n][a[1]][weight]
    else:
      w_ij = G[n][a[0]][0][weight]
      w_ih = G[n][a[1]][0][weight]
      
    a_ij = G.has_edge(n, a[0])
    a_ih = G.has_edge(n, a[1])
    a_hj = G.has_edge(a[1], a[0])
    Cw_barret = Cw_barret + ((w_ij + w_ih)/2)*(a_ij * a_ih * a_hj)
  
  k_i, s_i = strengh_barret(G, n, weight, sentido, vizinhos = vizinhos)
  if (s_i == 0) | (k_i <= 1) :
    Cw_barret = 0
  else:
    Cw_barret = (1/((k_i - 1) * s_i))  * Cw_barret
  
  return Cw_barret



def degree_histogram_directed(G, in_degree=False, out_degree=False):
    """Return a list of the frequency of each degree value.

    Parameters
    ----------
    G : Networkx graph
       A graph
    in_degree : bool
    out_degree : bool

    Returns
    -------
    hist : list
       A list of frequencies of degrees.
       The degree values are the index in the list.

    Notes
    -----
    Note: the bins are width one, hence len(list) can be large
    (Order(number_of_edges))
    """
    nodes = G.nodes()
    if in_degree:
        in_degree = dict(G.in_degree())
        degseq=[in_degree.get(k,0) for k in nodes]
    elif out_degree:
        out_degree = dict(G.out_degree())
        degseq=[out_degree.get(k,0) for k in nodes]
    else:
        degseq=[v for k, v in G.degree()]
    dmax=max(degseq)+1
    freq= [ 0 for d in range(dmax) ]
    for d in degseq:
        freq[d] += 1
    return freq

  
  
def opsahl(G, u=None, distance=None, normalized=True):
  """Calculates closeness centrality using Tore Opsahl's algorithm.
    See: https://toreopsahl.com/2010/03/20/closeness-centrality-in-networks-with-disconnected-components/
    The parameters are just like those in networkx.closeness_centrality()
    
    Parameters
    ----------
    G (graph) –
        A NetworkX graph.
    u (node, optional) –
        Return only the value for node u.
    distance (edge attribute key, optional (default=None)) –
        Use the specified edge attribute as the edge distance
        in shortest path calculations.
    normalized (bool, optional) –
        If True (default) normalize by the number of nodes
        in the connected part of the graph.
    
    Returns
    -------
    nodes –
        Dictionary of nodes with closeness centrality as the value.
    """
    
  # Which function to use for Dijkstra's algorithm
  if distance is not None:
      path_length = functools.partial(nx.single_source_dijkstra_path_length,
                                      weight=distance)
  else:
      path_length = nx.single_source_shortest_path_length

  # Whether to calculate for all nodes or just one
  if u is None:
      nodes = G.nodes()
  else:
      nodes = [u]
  
  closeness_centrality = {}
  for n in nodes:
      sp = dict(path_length(G, n))
      totsp = functools.reduce(lambda x, y: x + 1/y if y!=0 else x, sp.values(), 0.0)
      
      if totsp > 0.0 and len(G) > 1:
          closeness_centrality[n] = totsp
          # normalize to number of nodes-1 in connected part
          if normalized:
              s = 1 / ( len(G) - 1 )
              closeness_centrality[n] *= s
      else:
          closeness_centrality[n] = 0.0
  
  # Return centralit(y/ies)
  if u is not None:
      return closeness_centrality[u]
  else:
      return closeness_centrality
     
def plot_basics(data, data_inst, fig, units, discrete = True, xmin = 1,linear_bins=True):
    from powerlaw import plot_pdf, Fit, pdf
    annotate_coord = (-.4, .95)
    ax1 = fig.add_subplot(n_graphs,n_data,data_inst)
    x, y = pdf(data, linear_bins=linear_bins)
    ind = y>0
    y = y[ind]
    x = x[:-1]
    x = x[ind]
    ax1.scatter(x, y, color='r', s=1)
    plot_pdf(data[data>0], ax=ax1, color='b', linewidth=2)
    from pylab import setp
    setp( ax1.get_xticklabels(), visible=False)

    if data_inst==1:
        ax1.annotate("A", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)

    
    from mpl_toolkits.axes_grid.inset_locator import inset_axes
    ax1in = inset_axes(ax1, width = "30%", height = "30%", loc=3)
    ax1in.hist(data, density=True, color='b')
    ax1in.set_xticks([])
    ax1in.set_yticks([])

    
    ax2 = fig.add_subplot(n_graphs,n_data,n_data+data_inst, sharex=ax1)
    plot_pdf(data, ax=ax2, color='b', linewidth=2)
    fit = Fit(data, xmin=xmin, discrete=discrete)
    fit.power_law.plot_pdf(ax=ax2, linestyle=':', color='g')
    p = fit.power_law.pdf()

    ax2.set_xlim(ax1.get_xlim())
    
    fit = Fit(data, discrete=discrete)
    fit.power_law.plot_pdf(ax=ax2, linestyle='--', color='g')
    from pylab import setp
    setp( ax2.get_xticklabels(), visible=False)

    if data_inst==1:
       ax2.annotate("B", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)        
       ax2.set_ylabel(u"p(X)")# (10^n)")
        
    ax3 = fig.add_subplot(n_graphs,n_data,n_data*2+data_inst)#, sharex=ax1)#, sharey=ax2)
    fit.power_law.plot_pdf(ax=ax3, linestyle='--', color='g')
    fit.exponential.plot_pdf(ax=ax3, linestyle='--', color='r')
    fit.plot_pdf(ax=ax3, color='b', linewidth=2)
    
    ax3.set_ylim(ax2.get_ylim())
    ax3.set_xlim(ax1.get_xlim())
    
    if data_inst==1:
        ax3.annotate("C", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)

    ax3.set_xlabel(units)
