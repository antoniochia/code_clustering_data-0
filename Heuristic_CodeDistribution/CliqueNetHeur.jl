
function CNsetAeqB(A)
  B = A + 0
  return B
end


function CNrandsol(v, kl, G)
  x = rand(1:kl, v)
  x = CNremempclust(x)
  kl2 = maximum(x)
  ncl = kl2
  for (i in 1:kl2)
    ind = find(a -> a == i, x)
    GRdc = G[ind, ind]
    NRdc = IncToAdj(GRdc)
    t = ind[findSubGconnected(NRdc)]
    while (length(t) < length(ind))
      ind = setminus(ind, t)
      x[ind] = ncl + 1
      ncl = ncl + 1
      GRdc = G[ind, ind]
      NRdc = IncToAdj(GRdc)
      t = ind[findSubGconnected(NRdc)]
    end # while
  end # for
  return x
end

function CNneighvnsA1(N, x, nexch)
  G = AdjToInc(N)
  newx = x + 0
  ncl = maximum(x)
  t = randperm(length(newx))
  for (it in 1:nexch)
    newx[t[it]] = ncl + it
  end # for
  l = ncl + nexch
  for (i in 1:ncl)
    ind = find(a -> a == i, newx)
    if (length(ind) > 0)
	    GRdc = G[ind,ind]
    	NRdc = IncToAdj(GRdc)
    	t = ind[findSubGconnected(NRdc)]
    	while (length(t) < length(ind))
      		ind = setminus(ind, t)
      		newx[ind] = l + 1
      		l = l + 1
      		NRdc = IncToAdj(G[ind,ind])
      		t = ind[findSubGconnected(NRdc)]
    	end # while
    end # if
  end # for
  newx = CNremempclust(newx)
  return newx
end

function CNfo(C, x) # x vettore che dice a quale cluster appartiene
                      # ciascun elemento
                      # funziona sia che i singoletti siano 0
                      # sia che abbiano un indicazione del cluster (testato)
                      # funziona anche se i cluster non sono numerati in sequenza
  v = length(x)
  if (size(C)[1] != v)
    return("Dimensione errata")
  end
  f = 0
  for i in 1:(v - 1)
    for (j in (i + 1):v)
      if (x[j] == x[i])
        f = f + C[i, j]
      end
    end
  end
  return f
end

function CNclustercomposition(x) # restituisce una matrice i cui ogni riga
                               # elenca i nodi di quella partizione
                              # x[i] = k significa che il nodo i (non ammette x[i] = 0)
  ncl = maximum(x)
  P = zeros(ncl, length(x))
  d = find(x -> x == 1, x)
  if( length(d) > 0)
    for ( j in 1:length(d))
      P[1,j] = d[j]
    end
  end
  for (i in 2:ncl)
   d = find(x -> x == i, x)
   if( length(d) > 0)
    for ( j in 1:length(d))
      P[i,j] = d[j]
    end
   end
  end
  return int32(P)
end


function CNFeasOut(v, N, ClustOut)  # da testare
  st = delfromset(v, ClustOut)
  ind = st[find(x -> x>0, st)]
  if (length(ind) < 1)
    return true
  else
    G = AdjToInc(N)
    GRdc = G[ind,ind]
    NRdc = IncToAdj(GRdc)
    ch2 = isGconnected(NRdc)
    return ch2
  end #if
end

function CNbestreloc(C, Net, x)  # Net[i,j] = k, significa che k è adiacente a i
  v = length(x)
  ncl = maximum(x)
  feasout = trues(v)
  feasin = falses(v, ncl)
  costout = zeros(v)       # costi uscita cluster
  costin = zeros(v, ncl)   # costi entrata cluster (matrice)

  P = CNclustercomposition(x) # mi da la composizione di ogni cluster
  for (i in 1:v)
      feasout[i] = CNFeasOut(i, Net, vec(P[x[i], :]))
      if (feasout[i] == true)
        it = 1
        while ( (it <= size(P)[2]) && (P[x[i], it] > 0) )
          if (P[x[i], it] != i) # considera C[i,j] solo se i neq j
            costout[i] = costout[i] + C[i, P[x[i], it]]
          end
          it = it + 1
        end # end while
      end # end if
  end # end for
  for (i in 1:v)
    if (feasout[i] == true)
      for (k in 1:ncl)
        if (x[i] != k)
          cnt = notemptyintersect(Net[i,:], P[k, :])  # il nodo i deve essere connesso ad almeno un nodo di P[k, ]
          feasin[i,k] = cnt
          if (cnt == true)
              costk = 0
              it = 1
              while ( (it <= size(P)[2]) && (P[k, it] > 0) )
                costk = costk + C[i, P[k, it]]
                it = it + 1
              end # while
              costin[i,k] = costk
          end#if
        end# if
      end # for k
    end# if
  end# for i
  return costout, costin, feasout, feasin
  end

function CNcompilemove(costout, costin, feasout, feasin)
  d = size(costin)
  ncl = d[2]
  v = d[1]
  delta = zeros(v, ncl + 1) # ultima colonna corrisponde a creare un singoletto
  for (i in 1:v)
      if (feasout[i] == false)
        delta[i, :] = +Inf16
      else
        delta[i, ncl + 1] = -costout[i]
        for (j in 1:ncl)
          if (feasin[i,j] == false)
            delta[i,j] = +Inf16
          else
            delta[i,j] = costin[i,j] - costout[i]
          end # if
        end # for
      end# if
  end # for
 return delta # l
end

function CNchoosemove(delta, x)
    dlt = delta + 0
    v = size(dlt)[1]
    moves = zeros(v, 2)
    gains = zeros(v)
    P = CNclustercomposition(x)
    it = 1
    while (it <= v)
      t = ind2sub(size(dlt), indmin(dlt)) # t = i,k: i entra in k
      if (dlt[t[1], t[2]] >= 0)
        return int32(moves), gains, dlt
      else
        moves[it,1] = t[1]
        moves[it,2] = t[2]
        gains[it] = delta[t[1],t[2]]
        #dlt[t[1],:] = Inf16 # i non puo essere piu spostato
        dlt[:, x[t[1]]] = Inf16 # kout non puo piu entrare nessuno
        if (t[2] < size(dlt)[2]) # se kin non deve essere creato allora...
          dlt[:,t[2]] = Inf16 # kin non puo piu essere modificato
        end # if
        for (i in 1:size(P)[2]) # tutti gli elementi di kout non si possono piu muovere
          if (P[x[t[1]], i] > 0)
            dlt[P[x[t[1]],i],:] = Inf16
          end # if
        end # for
        if (t[2] < size(dlt)[2]) # tutti gli elementi di kin non possono spostarsi
          for (i in 1:size(P)[2])
            if (P[t[2],i] > 0)
             dlt[P[t[2],i],:] = Inf16
            end # if
         end#for
        end# if
      end # if
      it = it + 1
    end # while
    return int32(moves), gains, dlt
  end

function CNupdatesol(moves, gains, x, fo)
    fn = fo + 0
    y = x + 0
    ncl = maximum(x)
    it = 1
    while ( it <= size(moves)[1] && moves[it, 1] > 0 )
      fn = fn + gains[it]
      if (moves[it, 2] <= maximum(x))
        y[moves[it,1]] = moves[it,2]
      else
        ncl = ncl + 1
        y[moves[it,1]] = ncl
      end #if
      it = it + 1
    end # while
    return fn, y
  end # function

function CNremempclust(x)
  ncl = maximum(x)
  uno = fill(1,length(x))
  crd = zeros(ncl)
  for (i in 1:ncl)
    crd[i] = sum(x -> x==i, x)
  end
  fll = find(x -> x >=1, crd) # which(crd >= 1)
  for (i in 1:length(fll))
    if (fll[i] != i)
      d = find(x -> x == fll[i], x)
      x[d] = i
    end
  end
  ncl = maximum(x)
  for (i in 1:length(x))
    if (x[i] == 0)
      ncl = ncl + 1
      x[i] = ncl
    end
  end
  return x
end

function CNlclopt(C, Net, x)
  locopt = false
  it = 1
  fo = CNfo(C, x)
  flag = 0 #message = "Tutto OK"
  while (locopt == false)
    # f = CNfo(C, x) # questo si puo togliere
    u = CNbestreloc(C, Net, x)
    delta = CNcompilemove(u[1], u[2], u[3], u[4])
    m = CNchoosemove(delta, x)
    if (m[2][1] == 0)
      locopt = true
    else
      #cnt = checksinglemove(C, m[1], m[2], x)
      #if (cnt[1] == 0)
      #  return cnt
      #end
      v = CNupdatesol(m[1], m[2], x, fo)
      #newx = CNremempclust(v[2])
      x = CNremempclust(v[2])
      fo = v[1]
      #checkfo = CNfo(C,newx)
      #if ( abs(fo - checkfo) > 0.1)
      #  flag = 1 # message = "Ancora problemi"
      #  locopt = true
      #  return flag, checkfo, newx, it
      #else
      #  x = newx
      #end # if
      it = it + 1
    end # if
  end #while
  return flag, fo, x, it
end


function CNneighvns(N, x, nexch)
  # newx = x + 0
  v = length(x)
  newx = fill(0, v)
  newx = x + 0
  ncl = maximum(x)
  it = 1
  while (it <= nexch)
    i = rand(1:v)
    olcluster = find(a -> a == newx[i], newx)
    feasout = CNFeasOut(i, N, olcluster)
    if (feasout == true)
      ncl = maximum(newx)
      k = rand(1:ncl)
      if (k == newx[i])
        newx[i] = ncl + 1
        it = it + 1
      else
        newcluster = find(a -> a == k, newx)
        feasin = notemptyintersect(N[i,:], newcluster)
        if (feasin == true)
          newx[i] = k
          it = it + 1
        end
      end
    end
    newx = CNremempclust(newx)
  end
  return newx
end

function CNlocHeur(C, Net; version = "RR", maxstart = 10)
  v = size(C)[1]
  G = AdjToInc(Net)
  x = CNrandsol(v, rand(3:7), G)
  sol = CNlclopt(C, Net, x)
  xbest = CNsetAeqB(sol[3])
  itbest = sol[4]
  it = sol[4]
  fbest = sol[2]
  if (sol[1] > 0)
      return sol
  end
  flag = 0
  itstart = 2
  vnsngb = 3
  maxvnsngb = int32(ceil(v/2))
  while(itstart <= maxstart)
    if (version == "RR")
      x = CNrandsol(v,rand(3:7), G)
    end
    if (version == "VNS")
      #xb = CNremempclust(xbest)
      x = CNneighvnsA1(Net, xbest, vnsngb)
    end
    sol = CNlclopt(C, Net, x)
    if (sol[1] > 0)
      flag = sol[1]
      return flag, fbest, xbest, itbest
    end
    itstart = itstart + 1
    it = it + sol[4]
    if (sol[2] < fbest)
      fbest = sol[2]
      xbest = CNsetAeqB(sol[3])
      itbest = it
      vnsngb = 3
    else
      vnsngb = vnsngb + 1
      if (vnsngb > maxvnsngb)
        vnsngb = 3
      end
    end
  end
  return flag, fbest, xbest, itbest
end

function checksinglemove(C, moves, gains, x)
  fo = CNfo(C, x)
  it = 1
  while (gains[it] < 0)
    newx = x + 0
    newx[moves[it,1]] = moves[it,2]
    newfo1 = CNfo(C, newx)
    newfo2 = fo + gains[it]
    if (abs(newfo1 - newfo2) > 0.1)
      message = "Mossa sbagliata"
      return 1, message, it
    else
      it = it + 1
    end
  end
  return 0, "Mosse singole ok"
end

function checkfeasol(G, x)
  a = false
  P = CNclustercomposition(x)
  for (i in 1:size(P)[1])
    ind = vec(P[i, find(x -> x>0, P[i,:])])
    GRdc = G[ind,ind]
    NRdc = IncToAdj(GRdc)
    a = isGconnected(NRdc)
    #if (a == false)
    #  return a, i
    #end
  end
  return a #, 0
  end

