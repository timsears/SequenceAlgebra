module SequenceAlgebra where

-- Code from the the author of the following
-- K. Clenaghan, Code for "In Praise of Sequence (co)Algebra
-- and its implementation in Haskell"  (2018)
-- The core part is based on M.D. McIlroy, Power series, power serious, 
-- J Functional Prog. 9(3) pp325-327, 1999

-- This code is directly derived from a file shared by K. Clenaghan, adding cabal and nix support.

import Prelude hiding (cycle)
import Data.Ratio
default (Integer, Rational)

-------------- Sequences as Num ------------------------

instance (Eq a, Num a) => Num [a]  where
     negate            = map negate 
     f+[]              = f
     []+g              = g
     (f0:f')+(g0:g')   = f0+g0 : f'+g'
     []*_              = []
     (0:f')*g          = 0: f'*g
     _*[]              = []
     (f0:f')*g@(g0:g') = f0*g0:(f0*|g' +f'*g)
     fromInteger c     = [fromInteger c]
     abs _             = error "abs not defined on sequences"
     signum _          = error "signum not defined on sequences"

-------------- x, Catalans, scalar product -------------

x::Num a => [a]
x=[0,1]
c = 1 + x*c^2

infix 7 *|
--(*|)  :: Num a => a -> [a] -> [a]
a *| f = map (a*) f

-------------- Sequences as Fractional -----------------

instance (Eq a, Fractional a) => Fractional [a] where
    recip f           = 1/f
    _/[]              = error "divide by zero."
    []/_              = []
    (0:f')/(0:g')     = f'/g'
    (_:f')/(0:g')     = error "divide by zero"
    (f0:f')/g@(g0:g') = let q0=f0/g0 in q0:((f' - q0*|g')/g)
    fromRational c    = [fromRational c]

-------------- Square root -----------------------------

sqroot (0:0:f'') = 0:sqroot f''
sqroot f@(1:f')  = 1:(f'/(1+sqroot f))

-------------- Composition and converse ----------------

-- o:: (Eq a, Num a) => [a]->[a]->[a]
[] `o` _              = []
(f0:f') `o` g@(0:g')  = f0: g'*(f' `o` g)
(f0:f') `o` g@(g0:g') = [f0] + (g0*|(f' `o` g))+ (0:g'*(f' `o` g))

converse(0:f') = g where g = 0: 1/(f' `o` g)

-------------- Exponential-to-Ordinary -----------------

e2o f    = zipWith (*) f facs 
o2e f    = zipWith (/) f facs
from     :: Num a=>a->[a]
from     = iterate (+1)
nats, pos, zeros, facs :: Num a=>[a]
nats     = from 0
pos      = from 1
zeros    = 0:zeros
facs     = scanl (*) 1 pos

-------------- Derivation and Integration ---------------

deriv f = zipWith (*) pos (tail f)
integ f = 0:zipWith (/) f pos

-------------- Exp and Star -----------------------------

expx  :: (Eq a,Fractional a) => [a]
expx  = 1 + integ expx
starx :: (Eq a, Num a) => [a]
starx = 1 : starx

-------------- Rational-to-Whole -----------------------

makeWhole r  = case properFraction r of
                (n,0)         -> n
                otherwise     -> error "not whole"
makeAllWhole = map makeWhole 
takeW n      = take n . makeAllWhole

-------------- Core Sequences --------------------------

lgnx,  tanx, cosx, sinx, secx, coshx, tanhx::  (Eq a, Fractional a) => [a]

lgnx     = integ (1/(1+x))
sinx     = integ cosx
cosx     = 1-(integ sinx)
tanx     = integ (1+(tanx^2))
secx     = 1+integ (secx * tanx) -- 1/cosx
coshx    = 1+integ sinhx
sinhx    = integ coshx
tanhx    = sinhx/coshx
gdx      = integ (1/coshx) 

arctanx, sinhx :: (Eq a, Fractional a) => [a]

arctanx  = integ (1/(1+x^2)) 
arcsinx  = integ (1/sqroot(1-x^2))

-------------- Bivariate u, z, pascal, powerSums -------
    
u,z, pascal :: (Eq a, Num a) => [[a]]
u      = [[0],[1,0]]
z      = [[0],[0,1]]
pascal = starx `o` (u+z)

ebinom     =  expx `o` (z+u*z)
powerSums  = ((expx `o` (u*z))-1)/((expx `o` z)-1)
list, pluralList :: (Eq a, Num a) => [a]
list       = starx
pluralList = list - x - 1
schroeder  =  z + u*(pluralList `o` schroeder)

-------------- Bivariate utilities ---------------------

select s t   = zipWith take s t
selectW s t  = map makeAllWhole (zipWith take s t)
unDiag       :: Num a=> [[a]]->[[a]]
unDiag       = transpose . padTri
unDiage2o    :: Fractional a=> [[a]]->[[a]]
unDiage2o    = unDiag . (map e2o)
padTri t     = zipWith padRight t [1..]
padRight r k = r++(take (k-(length r)) zeros)    
    
transpose []   = []
transpose m    = c:transpose m'
  where (c,m') = foldr detachHead ([],[]) m
        detachHead   ([r0]) b  = (r0:fst b,snd b)
        detachHead   (r0:r') b = (r0:fst b,r':snd b)

-------------- Partial Differentiation and Test case ---------------

dz s = map deriv (tail s)
du s = map (reverse . deriv . reverse) (tail (padTri s))

set, nonEmptySet, emptySet :: (Eq a, Fractional a) => [a]  
set         = expx
emptySet    = 1
nonEmptySet = set - emptySet
parts       = set `o` (u*(nonEmptySet `o` z))

allEq c r    = foldr (\a b-> (a==c) && b) True r
allZeros s t = allEq True (map (allEq 0) (select s t))
-- allZeros [1..6] ((dz parts) - (u*parts +u*du parts))

-------------- Examples --------------------------------------------

cycle, perm :: (Eq a, Fractional a) => [a]

perm  = starx
lg g  = lgnx `o` (g-1)
cycle = lg starx

cayleyTree            = x*(set `o` cayleyTree)
connectedAcyclicGraph = cayleyTree - cayleyTree^2 / 2

-------------- Exponentiation and example -------------------------

infix 7 |^
--(|^)      :: (Eq a, Fractional a) => [a] -> a -> [a]
f |^ r    = expx `o` (r *| lg f)
legendre  = (1-2*u*z+z^2) |^ [-1%2]
hermite   = expx `o` (2*u*z-z^2)

-------------- Shuffle Product, Inverse, and Factorial examples ---

fac = 1+x*fac + x^2*(deriv fac)

cf_fac  = 1/(cfdenom 1)
  where cfdenom n = 1 - x*(2*n-1)-x^2*n^2*(1/(cfdenom (n+1))) 

infix 7 |><|    -- shuffle product
f@(f0:f') |><| g@(g0:g') = (f0 * g0): ((f' |><| g)+(f |><| g'))
_ |><| []                = []
[] |><| _                = []
-- f |><| g        = e2o ((o2e f) * (o2e g))
shInv f@(f0:f') = (1/f0): (-f' |><| ((shInv f) |><| (shInv f)))
-- fac = makeAllWhole (shInv (1-x))

-------------- Difference, Delta and Sigma   ----------------------

h2i s          = (1/(1+x)) |><| s
i2h s          = (1/(1-x)) |><| s
rh2i s@(s0:s') = s0: rh2i (delta s)
delta s        = (tail s) - s
sigma s        = x*starx*s
recur          :: (Eq a, Num a) => a -> [a]
recur a        = a:recur a
--fib          = 1/(1-x-x^2)
fib            = (take 2 (rb*[1,1]))/rb where rb = reverse (x^2-x-1)
-------------- Hadamard and Infiltration products -----------------

infix 7 |*|
f |*| g = zipWith (*) f g

infix 7 |^|
f@(f0:f') |^| g@(g0:g') = (f0*g0): ((f'|^|g)+(f|^|g'))+(f'|^|g')
_ |^| []                = []
[] |^| _                = []

infProd f g = h2i (i2h f |*| i2h g)

-------------- Falling factorial polynomials ----------------------

cycles        = set `o` (u* (cycle `o` z))
fall n m      = product [n-i|i<-take m nats]
alt           :: Num a => a->[a]
alt r         = r:alt (-r)
altMat m      = zipWith op (alt 1) m
                where op sign r = zipWith (*) (alt sign) r
monom2FacPoly = select [1..] (unDiage2o parts)
facMonom2Poly = select [1..] (altMat (unDiage2o cycles))

toFacPoly p   = sum (zipWith (*|) p monom2FacPoly)
fromFacPoly p = sum (zipWith (*|) p facMonom2Poly)
squaresFacPoly= o2e (take 4 (h2i [0,1,5,14,30,55]))

-------------- Logan polynomials ----------------------------------

logan      = (((sinx `o` z)+u*(cosx `o` z))/
             ((cosx `o` z)-u*(sinx `o` z)))

loganPolys = iterate (\p -> (1+x^2) * deriv p) x

-------------- Entringer triangle and ue2o ------------------------

prefixSums      = scanl (+) 0
suffixSums      = reverse.prefixSums.reverse
alternate f g a = a:alternate g f (f a)
entringer       = alternate suffixSums prefixSums [1]

gkpEnt =   (((sinx `o` u) + (cosx `o` u))/(cosx `o` (u+z)))
ue2o   :: Fractional a => [a]->[a]
ue2o   = reverse.e2o.reverse

-------------- Kozen, Silva, Moessner, Long, Paasche --------------

x2z rho     = sum (zipWith (\c zn->[[c]]*zn) rho zPowers)
              where zPowers = (iterate (*z) 1) 

moessnerT r  = iterate (\s -> (x2z (s!!r))*pascal) pascal
nats2Power r = [sn!!r!!r | sn <- moessnerT r]
moessnerH r  = iterate (\rho -> ((x2z rho)*pascal)!!r) 1

ksmlp h0 d = map last (scanl step h0 d)
 where step hn dn = ((x2z hn)*pascal)!!((length hn - 1)+dn)
                               
moessner r = ksmlp 1 (r:zeros)
long a b r = ksmlp [b,a-b] (r:zeros)
paascheFac = ksmlp 1 [1,1..]
superFac   = ksmlp 1 [1,2..]

-------------- xcot, and similar -----------------------------------

xcotx, xcothx :: (Eq a, Fractional a) => [a]
xcotx   = (x*cosx)/sinx
xcothx  = (x*coshx)/sinhx
xcth r  = [r]*(x*(coshx `o` ([r]*x)))/(sinhx `o` ([r]*x))


-------------- Gaussian Rationals ----------------------------------

infix  6  :&
data Gaussian a   = a :& a
                    deriving (Eq, Read, Show)

i    :: (Eq a, Num a) => Gaussian a
i    = 0:&1
ix   :: (Eq a, Num a) => [Gaussian a]
ix = [0,i]

instance  (Num a) => Num (Gaussian a)  where
    negate (x:&y)       = negate x :& negate y
    (x:&y) + (x':&y')   = (x+x') :& (y+y')
    (x:&y) * (x':&y')   = (x*x'-y*y') :& (x*y'+y*x')
    fromInteger n       = fromInteger n :& 0
    abs (x:&y)          = abs x :& abs y
    signum _            = error "signum not defined for Gaussians"

instance  (Fractional a) => Fractional (Gaussian a)  where
    (x:&y) / (x':&y')   =  (x*x'+y*y') / d :& (y*x'-x*y') / d
                           where d  = x'*x' + y'*y'
    fromRational a      =  fromRational a :& 0

conj        :: (Num a)  => Gaussian a -> Gaussian a
conj (x:&y) = x :& (-y)

makeReal (r:&0)  = r
makeReal _       = error "not real"
makeAllReal g    = map makeReal g

-------------- Shuffle ``ring'' ----------------------------------

newtype Shuffle a = S [a] deriving (Eq, Read, Show)
unS (S s) = s

instance (Eq a, Num a) => Num (Shuffle a) where
    negate (S s)  = S (negate s)
    (S s) + (S t) = S (s+t)
    (S s) * (S t) = S (s |><| t)
    fromInteger n = S [fromInteger n]
    abs _         = error "abs undefined on Shuffle seq"
    signum _      = error "signum undefined on shuffle seq"

secNums = 1:unS (S secNums*(S (e2o tanx)))

-------------- Maclaurin and Taylor expansions -------------------

maclaurin f = o2e (map head (iterate deriv f))
taylor f    = map o2e (zp (map (`o` u) (iterate deriv f)))
              where zp (g0:g') = g0 + z* (zp g')

bsinx, tsinx :: [[Rational]]
bsinx = sinx `o` (u+z)
tsinx = taylor sinx
 
------------------ univarite counting sequences ------------------

-- commented-out defns have already been introduced above

nonLoopCycle, pluralSet, singletonSet, derangement :: (Eq a, Fractional a) => [a]
 
nonEmptyList :: (Eq a, Num a) => [a]
-- emptySet       = 1
-- nonEmptySet    = set - emptySet
-- cycle          = lg starx
singletonList     = x
nonEmptyList      = list - 1
singletonSet      = x

ordPair           = x^2
fibonacci         = list `o` (singletonList + ordPair)
oneOrTwoCycle     = x + x^2/2
involution        = set `o` oneOrTwoCycle
nonLoopCycle      = cycle - x
derangement       = set `o` nonLoopCycle
permutation       = derangement * set
pluralSet         = nonEmptySet - singletonSet
setPartition      = set `o` nonEmptySet
oddNumberOfParts  = coshx `o` nonEmptySet
evenSizedParts    = set `o` sinhx
catalanTree       = x*(list `o` catalanTree)
motzkinTree       = x*(1+motzkinTree +motzkinTree^2)
connectedMapping  = cycle `o` cayleyTree
mapping           = set `o` connectedMapping
fixedPointFree    = set `o` (nonLoopCycle `o` cayleyTree)
idempotent        = set `o` (x*set)
partialMapping    = mapping*(set `o` cayleyTree)
surjection        = list `o` nonEmptySet
acyclicGraph      = set `o` connectedAcyclicGraph
zigzag            = 2*(tanx + secx)

-- pluralList        = list - x - 1
-- cayleyTree        = x*(set `o` cayleyTree)
-- connectedAcyclicGraph = cayleyTree - cayleyTree^2 / 2

------------------ bivarite sequences ----------------------------

intCompositions = list `o` (u*(nonEmptyList `o` z))
catalanLeaves   = u*z + z*(nonEmptyList `o` catalanLeaves)
cayleyLeaves    = u*z + z*(nonEmptySet `o` cayleyLeaves)
permFixedPts    = (derangement `o` z) * (set `o` (u*z))
ascents         = list `o` (z+(pluralSet `o` (u*z-z))/(u-1))
valleys         = r/(r - (tanhx `o` (z*r))) where r = sqroot (1-u)
bernoulli       = (z*(expx `o` (u*z)))/((expx `o` z)-1)
chebyshev       = (1-u*z)/(z^2-2*u*z+1)
laguerre        = (1/(1-z))*(expx `o` ((-u*z)/(1-z)))
meixner         = ((1+z^2) |^ [-1%2])*(expx `o` (u*(arctanx `o` z)))

-- ebinom      =  expx `o` (z+u*z)
-- powerSums   = ((expx `o` (u*z))-1)/((expx `o` z)-1)
-- parts       = set `o` (u*(nonEmptySet `o` z))
-- cycles      = set `o` (u* (cycle `o` z))
-- schroeder   =  z + u*(pluralList `o` schroeder)
-- legendre    = (1-2*u*z+z^2) |^ [-1%2]
-- hermite     = expx `o` (2*u*z-z^2)

