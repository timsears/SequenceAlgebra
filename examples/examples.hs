module Main where
import Prelude hiding (cycle)
import SequenceAlgebra

main :: IO ()
main = sequence_ $ examples

printN :: Show a => Int -> [a] -> IO ()
printN n v = print $ take n v

examples = [ printN 8 ( 3 * x^2 - x + 4)
           , printN 8 expx 
           , printN 8 lgnx
           , printN 8 (lg list)
           , printN 10
             ( 1 / 4 *  ( 1 + x - sqroot (1 - 6 * x + x^2)))
           -- from C
           , print $ select [1..6] (unDiag schroeder)
           , print "From Table 2"
           , printN 5 emptySet
           , printN 5 singletonSet
           , printN 5 nonEmptyList
           ]

