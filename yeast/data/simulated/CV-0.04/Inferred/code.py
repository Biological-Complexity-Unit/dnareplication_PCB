#!/usr/bin/env python3
import sys
import math
import os
import random
import time
import itertools
import bz2

# Range of grow rate Lambda / v
# in choosing above doubling time is between 50 minuts to 150 minutes and speed between 10 bps-1 to 1000 bps-1, ln(2) is involved
kmin=0.000000077
kmax=0.000023105

# Range of initiation rates I / v
Rmin=0.000000111
Rmax=0.0016666667
#in choosing above fork firing time is between 1 minuts to 150 minutes and speed between 10 bps-1 to 1000 bps-1, ln(2) is not involved

# Chromosome names
#Chromosomes = [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
#                "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16" ]
Chromosomes = [ "chr4" ]

# Chromosome lengths
Length = {
  "chr1":  252221,
  "chr2":  844051,
  "chr3":  341087,
  "chr4":  1575129,
  "chr5":  572496,
  "chr6":  277781,
  "chr7":  1092105,
  "chr8":  564939,
  "chr9":  430780,
  "chr10": 762303,
  "chr11": 683383,
  "chr12": 1084257,
  "chr13": 929095,
  "chr14": 793416,
  "chr15": 1108329,
  "chr16": 986200,
  "chrp2": 6300,
  "chrM":  94871
}

# Genome length
L = sum([Length[chr] for chr in Chromosomes])
print("Genome length:", L)

# Chromosome starting points
Start = dict()
start = 1
for chr in Chromosomes:
    Start[chr] = start
    start = start + Length[chr]

# Find chromosome for genomic location
def getchr(l):
  chri = 0
  chr = Chromosomes[chri]
  while (l >= Start[chr] + Length[chr]):
    chri=chri+1
    chr = Chromosomes[chri]
  return chr

# Window size
dL=1000

# Check whether the list contains elements in non-decreasing order
def isNonDecreasing(l):
  for i, el in enumerate(l):
    if (i > 0) and (l[i-1] > el):
      return False
  return True

# Sort origins by distance to a given location l.
# Prepares vectors y, rate that are ready to be passed to abundance()
# Arguments:
#   N: Number of origins
#   z: Vector of length N containing origin locations x_n
#   rate: Vector of length N containing origin initiation rates I_n / v [bp]
#   l: Location to sort by distance to
# Results:
#     y: Sorted vector of length N containing distances between l and origins z
#     z: Vector of length N containing corresponding origin locations x_n
#     rate: Vector of length N containing corresponding initiation rates I_n / v
def ordered(N, z, rate, f, w, wp, k, l):
  y=[ abs(l-zz) for zz in z ]
  if isNonDecreasing(y) and f is not None and w is not None and wp is not None:
    # Fast path, no need to sort
    return y, z, rate, f, w, wp
  else:
    # Have to sort. Create list of indices in ascending order of y
    idx = sorted(range(N), key=lambda i: y[i])
    # Generate sorted y, z and rate lists
    y = [ y[i] for i in idx ]
    z = [ z[i] for i in idx ]
    rate = [ rate[i] for i in idx ]
    f = [ r/k for r in rate ]
    w = [ 1.0 + r/k for r in itertools.accumulate(rate) ]
    wp = [ 1/w[i] - 1/w[i+1] for i in range(len(w)-1) ]
    return y, z, rate, f, w, wp
  # Old O(N^2) code:
  #y=[0 for i in range(N)]
  #for n in range (0, N):
  #  for m in range (n+1, N):
  #    c=abs(l-z[m])
  #    d=abs(l-z[n])
  #    if (c<d):
  #      z[n],z[m]=z[m],z[n]
  #      rate[n],rate[m]=rate[m],rate[n]
  #  y[n]=abs(z[n]-l)
  #return y,z,rate

# Compute abundance for a single location l
# Arguments:
#  N: Number of origins
#  k: Growthrate Lambda / v [bp]
#  y: Sorted vector of length N containing distances between l and origins
#  rate: Vector of length N containing origin initiation rates I_n / v [bp]
# Result: Abundance P(l)
def abundance1(N, y, rate, f, w, wp, k):
  assert N == len(y)
  assert N == len(rate)
  H = math.exp(-k*y[0]) / (1.0 + 1.0/f[0])
  n = 1
  Yn = y[0] * f[0]
  while n < N:
    Yn = Yn + y[n] * f[n]
    H = H + wp[n-1] * math.exp(-k*(y[n] * w[n] - Yn))
    n = n + 1
  return H
  #H=math.exp(-k*y[0])
  #aa=0.0
  #bb=0.0
  #for n in range (0,N-1):
  #  assert y[n] <= y[n+1]
  #  aa=aa+rate[n]
  #  bb=bb+rate[n]*y[n]
  #  a=(k+aa)*y[n]-bb
  #  b=(k+aa)*y[n+1]-bb
  #  H=H-((math.exp(-a)-math.exp(-b))*(k/(aa+k)))
  #aa=aa+rate[N-1]
  #bb=bb+rate[N-1]*y[N-1]
  #a=(k+aa)*y[N-1]-bb
  #H=H-(k*math.exp(-a)/(aa+k))
  #return H

# Update abundance P for selected chromosomes
# Arguments:
#  C: Chromosomes of genomic locations
#  I: Genomic locations
#  P: Abundances at genomic locations
#  chrs: Chromosomes to update
#  n: Dictionary of number of origins on each chromosome
#  z: Vector of length n containing origin locations (chr_n, x_n)
#  rate: Vector of length n containing origin initiation rates I_n / v [bp]
#  k: Growthrate Lambda / v [bp]
# Result: Abundance P(l)
def update_abundance(C,I,P,chrs,n,z,rate,k):
  assert len(C) == len(I)
  assert len(P) == len(I)
  assert len(z) == len(rate)
  # Evaluate P(l) on dl-spaced grid
  Pnorm = 0.0
  chr = None
  chr_update = None
  for i in range(len(I)):
    # If the current chromosome changed, update chromosome origin lists
    if C[i] != chr:
      chr = C[i]
      chr_update = (chr in chrs)
      if chr_update:
        # Extract origins on current chromosome
        chr_idx = [i for i in range(len(z)) if z[i][0]==chr]
        chr_n = len(chr_idx)
        assert chr_n == n[chr]
        chr_z = [z[i][1] for i in chr_idx]
        chr_rate = [rate[i] for i in chr_idx]
        chr_f = None
        chr_w = None
        chr_wp = None
      else:
        chr_idx, chr_n, chr_z, chr_rate, chr_f, chr_w, chr_wp = None, None, None, None, None, None, None
    # Update abundances if necessary
    if chr_update:
      l = I[i]
      chr_y, chr_z, chr_rate, chr_f, chr_w, chr_wp = ordered(chr_n, chr_z, chr_rate, chr_f, chr_w, chr_wp, k, l)
      P[i] = abundance1(chr_n, chr_y, chr_rate, chr_f, chr_w, chr_wp, k)
    # Update norm
    Pnorm = Pnorm + (P[i]*dL)
  return Pnorm

def initial(Nc, Nu):
  # Place Nc origins on each chromosome
  z=[ (chr, random.randint(Start[chr], Start[chr]+Length[chr]-1)) for chr in Chromosomes for i in range(Nc) ]
  # Place additional Nu origins randomly on the genome
  if Nu > 0:
    x=[ random.randint(1,L) for i in range(Nu) ]
    z=z + [ (getchr(l), l) for l in x ]
  # Compute number of origins per chromosome    
  n={ chr: sum([1 for zz in z if zz[0]==chr]) for chr in Chromosomes }
  # Draw random rates
  rate=[ random.uniform(Rmin,Rmax) for zz in z ]
  # Draw random growth rate
  k=random.uniform(kmin,kmax)
  return n,z,rate,k

# Randomly mutate configuration comprising origin locations, rates and growth rate k
# Arguments
#   nn: Current number of origins
#   zz: Current origin locations x_n
#   ratee: Current origin initiation rates I_v / v
#   kk: Current growth rate Lambda / v
# Results:
#   nn: Updated number of origins
#   zz: Updated origin locations x_n
#   ratee: Updated origin initiation rates I_v / v
#   kk: Updated growth rate Lambda / v
#   chrs: Chromosomes on which abundace predictions changed
def dynamics(nn,zz,ratee,kk) :
  Z1=0.24
  Z2=0.48
  Z3=0.72
  Z4=0.96
  r=random.random()
  chrs = []
  if(r<Z1):
    #a randomly chosen origin is relocated to randomly chosen location on the same chromosome (without affecting its rate)
    zz = zz.copy()
    m=random.randint(0,len(zz)-1)
    chr = zz[m][0]
    l=random.randint(Start[chr],Start[chr]+Length[chr]-1)
    zz[m] = (chr, l)
    chrs = [ chr ]
  elif(r>=Z1 and r<Z2):
    #rate increasing or decreasing
    ratee = ratee.copy()
    m=random.randint(0,len(zz)-1)
    ratee[m] = random.uniform(Rmin,Rmax)
    chrs = [ zz[m][0] ]
  elif(r>=Z2 and r<Z3):
    #removing an origin
    if(len(zz) > 1) :
      zz, ratee, nn = zz.copy(), ratee.copy(), nn.copy()
      m = random.randint(0,len(zz)-1)
      chr = zz[m][0]
      if (nn[chr] > 1):
        zz[m] = zz[len(zz)-1]
        ratee[m] = ratee[len(ratee)-1]
        zz = zz[:len(zz)-1]
        ratee = ratee[:len(ratee)-1]
        nn[chr] = nn[chr]-1
        chrs = [ chr ]
  elif(r>=Z3 and r<Z4):
    #adding an origin
    zz, ratee, nn = zz.copy(), ratee.copy(), nn.copy()
    l=random.randint(1,L)
    chr = getchr(l)
    zz.append((chr, l))
    ratee.append(random.uniform(Rmin,Rmax))
    nn[chr] = nn[chr]+1
    chrs = [ chr ]
  elif(r>=Z4):
    #changing the growth rate
    kk = random.uniform(kmin,kmax)
    chrs = Chromosomes
  return nn,zz,ratee,kk,chrs

def energy(C, I, P, Pnorm, Q, sigma, n, full_likelihood=True):
  # Evaluate energy for initial configuration.
  # We use the AIC as our energy function, and assume independent Gaussian errors.
  # The likelihood is thus the product of
  #   exp( - (P[l]/rho - Q[l])^2 / sigma[l])^2 / sqrt(2 * pi * sigma[l]^2)
  # over all locations l, and the AIC is therefore 
  #   2*(2N + 1) + sum [l] log(2 * pi * sigma[l]^2) + ((P[l] - Q[l]) / sigma[l])^2.
  
  # AIC penalty for parameter k
  akaike = 2.0*1.0
  chr = None
  for i in range(len(I)):
    # If the current chromosome changed, include AIC penalty for nr. of origins
    if C[i] != chr:
      chr = C[i]
      akaike = akaike + 2.0 * (2.0 * n[chr])
    # Sum negative likelihoods
    a = (P[i]/Pnorm-Q[i])/sigma[i]
    akaike = akaike+(a*a)
    if full_likelihood:
      akaike = akaike + math.log(2.0*math.pi*sigma[i]*sigma[i])
  return akaike

def read_observed(inputfnamepattern, empirical_error_est=False):
  C=[]
  I=[]
  Q=[]
  sigma=[]
  for chr in Chromosomes:
    with bz2.open(inputfnamepattern % (chr), "r") as inputfile:
      for j in inputfile:
        aa=j.split()
        l=Start[chr]+int(float(aa[0]))-1
        C.append(chr)
        I.append(l)
        Q.append(float(aa[1]) * float(Length[chr]) / float(L))
        sigma.append(float(aa[2]) * float(Length[chr]) / float(L))
  print("Observed abundances AUC: ", "%.5f" % (sum(Q)*dL))
  avgsd = math.sqrt(sum([ pow(Q[i]-Q[i+1], 2) / 2 for i in range(1, len(Q)-1) ]) / (len(Q)-1))
  avgsigma = sum(sigma) / len(sigma)
  print("Average of provided error estimate (sigma): ", "%.3E" % avgsigma)
  print("Empirical error estimate: ", "%.3E" % avgsd)
  if empirical_error_est:
    print("Using empirical errors")
    sigma = [avgsd] * len(I)
  else:
    print("Using provided error estimate")
  return C,I,Q,sigma

def output_results(C, I, P, Pnorm, Q, sigma, n, z, rate, k, prefix=""):
  abdfname = "abundance.dat"
  originsfname = "origins.dat"
  paramsfname = "parameters.dat"

  with open("%s%s" % (prefix, abdfname), "w") as abdfile:
    print("x\tseq\tlocation\tP\tA\tQ\tsigma", file=abdfile)
    for i in range(len(I)):
      chr = C[i]
      l = I[i]
      print("%d\t%s\t%d\t%f\t%E\t%E\t%E" % (l, chr, l-Start[chr]+1, P[i], P[i]/Pnorm, Q[i], sigma[i]), file=abdfile)

  with open("%s%s" % (prefix, originsfname), "w") as originsfile:
    print("x\tseq\tlocation\tI/v", file=originsfile)
    for i in range(len(z)):
      chr = z[i][0]
      loc = z[i][1]
      print("%d\t%s\t%d\t%E" % (loc, chr, loc-Start[chr]+1, rate[i]), file=originsfile)
    
  with open("%s%s" % (prefix, paramsfname), "w") as paramsfile:
      print("norigins\tk", file=paramsfile)
      print("%d\t%E" % (len(z), k), file=paramsfile)

def main():
  # Arguments
  if len(sys.argv) < 7:
    print("Arguments: <seed> <T> <Tfinal> <jmax> <javr> <errest> <prefix>", file=sys.stderr)
    print("  seed: Random seed", file=sys.stderr)
    print("  T: Initial temperature", file=sys.stderr)
    print("  Tfinal: Final temperature", file=sys.stderr)
    print("  jmax: Number of interations", file=sys.stderr)
    print("  javr: Reporting interval", file=sys.stderr)
    print("  errest: Error estimate to use, 'provided' or 'empirical'", file=sys.stderr)
    print("  prefix: Output file prefix", sys.stderr)
  seed = int(sys.argv[1])
  T = float(sys.argv[2])
  Tfinal = float(sys.argv[3])
  jmax = int(sys.argv[4])
  javr = int(sys.argv[5])
  errest = str(sys.argv[6])
  if errest not in ["provided", "empirical"]:
    raise "Invalid error estimate"
  if len(sys.argv) >= 8:
    prefix = str(sys.argv[7])
  else:
    prefix = ""

  # Parameters
  print("Seed: %d" % seed)
  print("T: %E" % T)
  print("Tfinal: %E" % Tfinal)
  print("jmax: %E" % float(jmax))
  print("javr: %d" % javr)
  print("errest: %s" % errest)
  mult = math.exp(math.log(Tfinal/T)/jmax)
  print("T factor per step: %.8f" % mult)
  inputfnamepattern = "../Truth/%s.dat.bz2"
  print("Output prefix: %s" % prefix)
  logfname = prefix + "data.dat"

  # Make reproducible
  random.seed(seed)

  # Read observed abundances
  C,I,Q,sigma = read_observed(inputfnamepattern, empirical_error_est=errest == "empirical")

  # Generate initial configuration
  n,z,rate,k = initial(2, 0)
  P=[0] * len(I)
  Pnorm = update_abundance(C,I,P,Chromosomes,n,z,rate,k)
  akaike = energy(C, I, P, Pnorm, Q, sigma, n, full_likelihood=False)

  # Run simulated annealing
  logfile = open(logfname,"w")
  line = "j\tn\tpaccept\ttemperature\tavgenergy\tavgdelta\tk"
  print(line, file=logfile)
  logfile.flush()
  print(line)
  sys.stdout.flush()

  ja=0
  frac=0.0
  asum=0.0
  for j in range(1, jmax):
    # Save current configuration
    nold, zold, rateold, kold, Pold, Pnormold, akaikeold = n, z, rate, k, P, Pnorm, akaike
    # Randomly perturb configuration
    n, z, rate, k, chrs = dynamics(n, z, rate, k)
    # Evaluate energy
    P = P.copy()
    Pnorm = update_abundance(C, I, P, chrs, n, z, rate, k)
    akaike = energy(C, I, P, Pnorm, Q, sigma, n, full_likelihood=False)
    # Compute acceptance probability p
    if(akaike < akaikeold):
      p=1.0
    else:
      p=math.exp(-(akaike-akaikeold)/T) 
    # Update temperature-scaled improvement average
    frac=frac+((akaike-akaikeold)/T)
    # Accept with probability p
    r=random.random()
    if(r>p) :
      # Reject perturbation
      n, z, rate, k, P, Pnorm, akaike = nold, zold, rateold, kold, Pold, Pnormold, akaikeold
    else:
      # Accept perturbation
      ja = ja+1
    # Update energy average
    asum += akaike
    # Log average acceptance probability every javr steps
    if (j % javr == 0) :
      line = "%d\t%d\t%.4f\t%.5E\t%.3E\t%.3E\t%.3E" % (j, len(z), ja/javr, T, asum/javr, frac/javr, k)
      print(line, file=logfile)
      logfile.flush()
      print(line)
      sys.stdout.flush()
      output_results(C, I, P, Pnorm, Q, sigma, n, z, rate, k, prefix="%sj%d-" % (prefix, j))
      ja = 0
      frac = 0.0
      asum = 0.0
    # Scale temperature
    T = T * mult

  # Output final fit
  output_results(C, I, P, Pnorm, Q, sigma, n, z, rate, k, prefix=prefix)

if __name__ == '__main__':
    main()
