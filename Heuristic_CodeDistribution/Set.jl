# Le seguenti sono istruzioni che permettono operazioni su insiemi,
# necessarie per la euristica sulle cliques

# set = [1, 5, 15, 0, 0]: ovvero l insieme non deve contenere altri elementi dopo lo 0
# isin(v, set) v = singoletto; set = insieme RESULT true/false
# notemptyintersect(neigh, clust) neigh, clust = insieme, RESULT true/false
# whereis(v, set) v = singoletto, set = insieme, RESULT = theindex i such that set[i] = v
# delfromset(v, set) v = singoletto, set = insieme RESULT set2 = set - v
# isGconnected(P) P = matrice di adiecenze di un grafo, RESULT = true/false

function isin(v, set) # calcola se v sta dentro linsieme vettore set
  in = false
  i = 1
  while ( (in == false) && (i <= length(set)) && (set[i] > 0) )
    if (v == set[i])
      in = true
    else
      i = i + 1
    end
  end
  return in
end

function setminus(a, b) # calcola a\b insiemistico
  a1 = a + 0
  for (i in 1:length(b))
    k = whereis(b[i],a)
    if (k > 0)
      a1[k] = 0
    end
  end
  a2 = a1[find(x -> x > 0, a1)]
  return a2
end

function notemptyintersect(neigh, clust)  # stabilisce se i due insiemi hanno almeno un elemento in comune
  cnt = false                             # funziona anche se uno dei due e vuoto
  i = 1
  while ( (cnt == false) && (i <= length(neigh) && (neigh[i] > 0)) )
    if (isin(neigh[i], clust) == true)
      cnt = true
    else
      i = i + 1
    end
  end
  return cnt
end

function whereis(v, set)
  set1 = vcat(set, 0)  # aggiungo uno 0 per avere un puntatore alle fine di set
  posv = 0
  it = 1
  f = false
  while (f == false && set1[it] > 0)
    if (set1[it] == v)
        posv = it
        f = true
      else
        it = it + 1
      end #if
    end# while
  return posv
end

function delfromset(v, set)  # attento che un vettore con degli 0 finali
  set2 = set + 0
  z = whereis(v, set2)
  #if (z == length(set2))
  set2[z] = 0
  #else
  if (z < length(set2))
    while ((z + 1) <= length(set2) && set[z + 1] > 0)
      set2[z] = set[z + 1]
      set2[z + 1] = 0
      z = z + 1
    end
  end
  return set2
end

function isGconnected(C) # C = nxm, i nodi devono essere numerati da 1 a n, la riga i contiene le connessioni del nodo j
  n = size(C)[1]
  if (n == 1)
    return true
  end
  m = size(C)[2]
  labeled = zeros(n)
  queue = zeros(n)
  queue[1] = 1
  labeled[1] = 1
  qbegin = 1 #primo termine > 0 non analizzato
  qend = 2 #primo termine = 0
  while (qbegin < qend)
    i = queue[qbegin]
    qbegin = qbegin + 1
    j = 1
    while (j <= m && C[i,j] > 0)
      if (labeled[C[i,j]] == 0)
        labeled[C[i,j]] = 1
        queue[qend] = C[i,j]
        qend = qend + 1
      end # if
      j = j + 1
    end  # while
  end # if
  if (sum(labeled) < (n - 0.1))
    return false
  else
    return true
  end
end

function findSubGconnected(C) # C = nxm, i nodi devono essere numerati da 1 a n, la riga i contiene le connessioni del nodo j
  n = size(C)[1]
  m = size(C)[2]
  labeled = zeros(n)
  queue = zeros(n)
  queue[1] = 1
  labeled[1] = 1
  qbegin = 1 #primo termine > 0 non analizzato
  qend = 2 #primo termine = 0
  while (qbegin < qend)
    i = queue[qbegin]
    qbegin = qbegin + 1
    j = 1
    while (j <= m && C[i,j] > 0)
      if (labeled[C[i,j]] == 0)
        labeled[C[i,j]] = 1
        queue[qend] = C[i,j]
        qend = qend + 1
      end # if
      j = j + 1
    end  # while
  end # if
  subG = 1
  for (i in 2:n)
    if (labeled[i] == 1)
      subG = vcat(subG, i)
    end # if
  end #for
  return subG
end

function makeGsym(G)
  G1 = int32(G + 0)
  n = size(G1)[1]
  for (i in 1:n)
    for (j in 1:n)
      if (G1[i,j] == 1)
        G1[j,i] = 1
      end
    end
    G1[i, i] = 0
  end
  return G1
end

function rmLastCol(A)
  n = size(A)[1]
  B = A[1:n;1:n]
  return B
end

function ARI(cl1, cl2)
  nrow = maximum(cl1)
  ncol = maximum(cl2)
  table = fill(0, nrow, ncol)
  tbin = fill(0, nrow, ncol)
  for (i in 1:length(cl1))
      table[cl1[i],cl2[i]] = table[cl1[i],cl2[i]] + 1
  end
  rsum = sum(table, 1)
  csum = sum(table, 2)
  for (i in 1:size(table)[1])
   for (j in 1:size(table)[2])
    tbin[i,j] = binomial(table[i,j], 2)
   end
  end
  rbin = fill(0, length(rsum))
  for (i in 1:length(rsum))
    rbin[i] = binomial(rsum[i],2)
  end
  cbin = fill(0, length(csum))
  for (i in 1:length(csum))
    cbin[i] = binomial(csum[i],2)
  end
  nbin = binomial(length(cl1), 2)
  num1 = sum(tbin)
  num2 = (sum(cbin)*sum(rbin))/nbin
  den1 = (sum(cbin) + sum(rbin))/2
  den2 = num2
  ari = (num1 - num2)/(den1 - den2)
  return ari
end # function

function cliqueaggregator(D, cl)
  n = maximum(cl)
  #Dagr = copy(D)
  Dagr = fill(0.0, size(D))
  for (i in 1:n)
  d = find(x -> x==i, cl)
  if (length(d) > 1)
    Drdc = D[:, d]
    medDr = median(Drdc, 2)
      for (j in 1:length(d))
        Dagr[:,d[j]] = medDr
      end
    else
      Dagr[:,d[1]] = D[:,d[1]]
    end # if
  end # for
  #for (j in 1:length(d))
  #  Dagr[:,d[j]] = medDr
  return Dagr
end

