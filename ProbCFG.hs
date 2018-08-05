module ProbCFG where

---------------------------------------------------------------

import qualified Data.Map as Map

-- A useful helper for debugging with Map. Feel free to ignore 
-- the implementation of this.
printMap :: (Show k, Show v) => Map.Map k v -> IO ()
printMap = putStr . unlines . map show . Map.toList

---------------------------------------------------------------

data Cat = S | NP | VP | N | D | V | PP | P | Adv deriving (Show,Eq,Ord)

data StrucDesc = Leaf Cat String | Binary Cat StrucDesc StrucDesc
                 deriving Show

type ProbCFG = ([(Cat,Double)],
                [((Cat,String),Double)],        -- terminal rules
                [((Cat,(Cat,Cat)),Double)],     -- nonterminal rules
                [Cat])

-- A minor variant of the grammar on page 384 of 
-- Manning and Schutze's "Foundations of Statistical NLP"
pcfg1 :: ProbCFG
pcfg1 = (   [(S,1.0)] ,
            [((P,"with"), 1.0), 
             ((V,"saw"), 1.0), 
             ((NP,"dogs"), 0.1), ((NP,"telescopes"), 0.18), ((NP,"saw"), 0.04), ((NP,"cats"), 0.18), ((NP,"hamsters"), 0.1)] ,
            [((S,(NP,VP)), 1.0),
             ((PP,(P,NP)), 1.0),
             ((VP,(V,NP)), 0.7), ((VP,(VP,PP)), 0.3),
             ((NP,(NP,PP)), 0.4)] ,
            [S,NP,VP,PP,P,V]
        )

-- Like above but reversed probabilities on the rules for expanding VP
pcfg2 :: ProbCFG
pcfg2 = (   [(S,1.0)] ,
            [((P,"with"), 1.0), 
             ((V,"saw"), 1.0), 
             ((NP,"dogs"), 0.1), ((NP,"telescopes"), 0.18), ((NP,"saw"), 0.04), ((NP,"cats"), 0.18), ((NP,"hamsters"), 0.1)] ,
            [((S,(NP,VP)), 1.0),
             ((PP,(P,NP)), 1.0),
             ((VP,(V,NP)), 0.3), ((VP,(VP,PP)), 0.7),
             ((NP,(NP,PP)), 0.4)] ,
            [S,NP,VP,PP,P,V]
        )

--------------------------------------------------
-- Utility functions for getting information from grammars.

probLookup :: (Eq a) => [(a,Double)] -> a -> Double
probLookup []           key = 0.0
probLookup ((x,y):rest) key = if key == x then y else probLookup rest key

allCats :: ProbCFG -> [Cat]
allCats (starting,ending,transitions,cats) = cats

startProb :: ProbCFG -> Cat -> Double
startProb (starting,ending,transitions,cats) = probLookup starting

endProb :: ProbCFG -> Cat -> String -> Double
endProb (starting,ending,transitions,cats) c s = probLookup ending (c,s)

trProb :: ProbCFG -> Cat -> (Cat,Cat) -> Double
trProb (starting,ending,transitions,cats) c (c1,c2) = probLookup transitions (c,(c1,c2))

-------------------------------------------------------------
-- Simple recursive definition of inside probabilities
leftChildCat :: ProbCFG -> Cat -> [Cat]
leftChildCat (starting,ending,transitions,cats) cat = 
    map (\((c,(l,r)),_) -> l) (filter (\((c,(l,r)),_) -> c == cat) transitions)

rightChildCat :: ProbCFG -> Cat -> [Cat]
rightChildCat (starting,ending,transitions,cats) cat = 
    map (\((c,(l,r)),_) -> r) (filter (\((c,(l,r)),_) -> c == cat) transitions)

naiveInside :: ProbCFG -> [String] -> Cat -> Double
naiveInside g [] c = undefined
naiveInside g [w] c = endProb g c w
naiveInside g ls c =
    sum (map (\(c1, c2) -> sum (map (\i -> trProb g c (c1,c2) * naiveInside g (take i ls) c1 * naiveInside g (drop i ls) c2) [1 .. (length ls - 1)])) [(x,y) | x <- leftChildCat g c, y <- rightChildCat g c]) 

-------------------------------------------------------------

-- A name for the specific Map type that we will use to represent 
-- a table of inside probabilities.
type InsideTable = Map.Map ([String],Cat) Double

fastInside :: ProbCFG -> [String] -> Cat -> Double
fastInside g sent c =
    Map.findWithDefault 0 (sent,c) (buildInsideTable g sent)

buildInsideTable :: ProbCFG -> [String] -> InsideTable
buildInsideTable g sent =
    fillCellsInside g Map.empty (cellsToFill g sent)

cellsToFill :: ProbCFG -> [String] -> [([String],Cat)]
cellsToFill g sent = [(chunk,cat) | chunk <- chunks sent, cat <- allCats g]

chunks :: [String] -> [[String]]
chunks ws = [take n (drop x ws) | n <- [1 .. (length ws)], x <- [0 .. (length ws - n)] ]

fillCellsInside :: ProbCFG -> InsideTable -> [([String],Cat)] -> InsideTable
fillCellsInside g tbl [] = tbl
fillCellsInside g tbl ((chunk,c):rest) = 
    -- First calculate the probability for (sentenceChunk,st)
    let result =
            case chunk of
            [] -> undefined
            [w] -> endProb g c w
            ls -> let probInside = \ys -> \cat -> Map.findWithDefault 0 (ys,cat) tbl in 
                sum (map (\(c1, c2) -> sum (map (\i -> trProb g c (c1,c2) * probInside (take i ls) c1 * probInside (drop i ls) c2) [1 .. (length ls - 1)])) [(x,y) | x <- leftChildCat g c, y <- rightChildCat g c])
    in
    -- Now add the calculated probability to the table, if nonzero
    let updatedTbl =
            if (result > 0) then
                Map.insert (chunk,c) result tbl
            else
                tbl
    in
    -- Continue on with the rest of the cells to be filled
    fillCellsInside g updatedTbl rest

-------------------------------------------------------------

-- A name for the specific Map type that we will use to represent 
-- a table of viterbi probabilities and backpointers.
type ViterbiTable = (Map.Map ([String],Cat) Double, Map.Map ([String],Cat) (Cat,Cat,Int))

getTriples :: [a] -> [b] -> [c] -> [(a,b,c)]
getTriples as bs cs = [(a,b,c) | a <- as, b <- bs, c <- cs]

fillCellsViterbi :: ProbCFG -> ViterbiTable -> [([String],Cat)] -> ViterbiTable
fillCellsViterbi g tbl [] = tbl
fillCellsViterbi g (tblProb, tblPtr) ((chunk,c):rest) =
    case chunk of
        [] -> undefined
        [w] -> let result = endProb g c w in
            if result > 0 then
                fillCellsViterbi g (Map.insert (chunk, c) result tblProb, tblPtr) rest
            else 
                fillCellsViterbi g (tblProb, tblPtr) rest
        ls -> 
            let viterbiProb = \ys -> \cat -> Map.findWithDefault 0 (ys,cat) tblProb in
            let triples = getTriples (leftChildCat g c) (rightChildCat g c) [1..(length ls-1)] in
            if length triples > 0 then 
                let (result, backptr) = maximum(map(\(c1, c2, i) -> (trProb g c (c1,c2) * viterbiProb (take i ls) c1 * viterbiProb (drop i ls) c2, (c1,c2,i))) triples) in
                if result > 0 then
                    fillCellsViterbi g (Map.insert (chunk, c) result tblProb, (Map.insert (chunk,c) backptr tblPtr)) rest
                else
                    fillCellsViterbi g (tblProb, tblPtr) rest
            else 
                fillCellsViterbi g (tblProb, tblPtr) rest

buildViterbiTable :: ProbCFG -> [String] -> ViterbiTable
buildViterbiTable g ls = 
    fillCellsViterbi g (Map.empty, Map.empty) (cellsToFill g ls)

extractStrucDesc :: ViterbiTable -> ([String],Cat) -> StrucDesc
extractStrucDesc (tblProb, tblPtr) (ls, c)  =
    case ls of 
        [] -> undefined
        [w] -> Leaf c w
        ls -> let (leftCat, rightCat, p) = Map.findWithDefault undefined (ls,c) tblPtr in
            let (prob, leftStr, rightStr) = maximum(map(\(left,right) -> (Map.findWithDefault 0 (left,leftCat) tblProb * Map.findWithDefault 0 (right,rightCat) tblProb, left, right)) (map (\i -> (take i ls, drop i ls)) [1..(length ls - 1)])) in
            Binary c (extractStrucDesc (tblProb, tblPtr) (leftStr, leftCat)) (extractStrucDesc (tblProb, tblPtr) (rightStr, rightCat))

