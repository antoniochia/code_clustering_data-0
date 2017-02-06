# Varie funzioni per generare cose


function signrand()
  z = rand()
  if z > 0.5
    b = 1
  else
    b = -1
  end
  return b
end

function rbinom(prob)
  z = rand()
  if  z < prob
    b = 1
  else
    b = 0
  end
  return b
end

function randsurvey(units, range)
  A = -range*eye(units)
  for i in 1:(units - 1)
    for j in (i+1):units
      A[i,j] = range*signrand()*rand()
      A[j,i] = A[i,j]
    end
  end
  return A
end

function symsurvey(units, relations, prob)
  S = zeros(units, units)
  for k in 1:relations
    A = zeros(units, units)
    for i in 1:(units-1)
    for j in (i + 1):units
        A[i,j] = rbinom(prob)
        A[j,i] = A[i,j]
      end
    end
    S = S + A
  end
  M = fill(relations, units, units)
  G = 2*relations*eye(units)
  S = M - 2*S - G
  return S
end

function cliquematrix(Survey)
  units = size(Survey)[1]
  feats = size(Survey)[2]
  S = zeros(units, units)
  for (i in 1:(units - 1))
    for (j in (i + 1):units)
      S[i,j] = sum(abs(Survey[i,:] - Survey[j,:]))
      S[j,i] = S[i,j]
    end
  end
  M = fill(feats, units, units)
  #G = 2*feats*eye(units)
  C = -M + 2*S #- G
  return C
end

function randbinetwork(n; prob = 0.5)  # calcola un grafo casuale classico
  N = zeros(n, n)
  for (i in 1:(n-1))
    for (j in (i + 1):n)
      N[i,j] = rbinom(prob)
      N[j,i] = N[i,j]
    end
  end
  return N
end

function rand2gn(ngroup, ncluster; pwithin = 0.15, pout = 0.05)
  G1 = randbinetwork(ngroup; prob = pwithin)
  G2 = randbinetwork(ngroup; prob = pwithin)
  G3 = zeros(ngroup, ngroup)
  for (i in 1:(ngroup - 1))
       for(j in (i + 1):ngroup)
         G3[i,j] = rbinom(pout)
       end
     end
  G1 = vcat(G1, G3)
  G2 = vcat(transpose(G3), G2)
  Gtot = hcat(G1, G2)
  return Gtot
end

function IncToAdj(G)  # funziona anche se G ha un solo nodo
  if (length(G) > 1)
  s = sum(G, 1)
  N = zeros( length(s), int32(maximum(s)))
  for (i in 1:length(s))
    pos = 1
    for (j in 1:length(s))
      if (G[i,j] > 0.01)
        N[i, pos] = j
        pos = pos + 1
      end
    end
  end
  else # nel caso il grafo sia formato da un singoletto
    N = zeros(1, 1)
  end
  return int32(N)
end

function AdjToInc(N)
  m = size(N)[1]
  n = size(N)[2]
  G = zeros(Int32, m, m)
  for (i in 1:m)
    for (j in 1:n)
	    if (N[i,j] > 0)
	      G[i, N[i,j]] = 1
        G[N[i,j], i] = 1
	    end
    end
  end
  return(G)
end

function GenRandProb(ngroup, nfeat, prob, pr_in, pr_out)
  cl = [fill(1, ngroup), fill(2, ngroup)]
  p1 = [prob, 1-prob]
  p2 = [1-prob, prob]
  D1 = Categorical(p1)
  D2 = Categorical(p2)
  MData = rand(D1, nfeat)
  for (i in 2:ngroup)
    MData = hcat(MData, rand(D1, nfeat))
  end
  for (i in 1:ngroup)
    MData = hcat(MData, rand(D2, nfeat))
  end
  C = cliquematrix(transpose(MData))
  GC = fill(1, size(C))
  G = rand2gn(ngroup, 2; pwithin = pr_in, pout = pr_out)
  return int32(G), int32(C)
end

