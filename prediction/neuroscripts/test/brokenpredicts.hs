{-# LANGUAGE NoMonomorphismRestriction #-} -- Don't let The Man keep you down
import FastA
import System.IO
import HSH
import Control.Monad
import Control.Applicative ((<$>),(<*>))
import Data.List
import System.Environment (getArgs)
import System.Directory (setCurrentDirectory)
---------------------------------
--constants
predDir = "/srv/data2/pipeline/prediction/neuroscripts/"
newExt oldName ext = (takeWhile (/= '.') oldName) ++ ext
fA projName = newExt projName "_final_assembly.fasta"
--------------------------------------------------
--Main routine
main = do
 projName <- head <$> getArgs
 fileNames <- formatForPredictors (fA projName)
 sequence_ $ [runSignalP,runTargetP,runTMHMM] <*> [fileNames]
-------------------------------------
runSignalP fileNames = do
 mapM_ runIO $ map sigP fileNames
  where sigP file = predDir ++ "signalp-3.0/signalp -t euk -f short " ++ file ++ " > " ++ (newExt file "_signalp.out")
-----------------------------------------
runTargetP fileNames = do
 mapM_ runIO $ map targP fileNames 
  where targP file = predDir ++ "targetp-1.1/targetp -N " ++ file ++ " > " ++ (newExt file "_targetp.out")
----------------------------------------------------- 
runTMHMM fileNames = do
 mapM_ runIO $ map tmhmm fileNames
tmhmm file = "cat " ++ file ++ " | decodeanhmm -f " ++ predDir ++ "tmhmm-2.0c/lib/TMHMM2.0.options -modelfile " ++ predDir ++ "tmhmm-2.0c/lib/TMHMM2.0.model | " ++ predDir ++ "tmhmm-2.0c/bin/tmhmmformat.pl > " ++ (newExt file "_tmhmm.out")
-----------------------------------------------------

--Parsing
data TargetPResult = Tp String Int deriving (Show,Eq)
instance Prediction TargetPResult where
 contgName (Tp name _) = name

data SignalPResult = Sp String Float Float deriving (Show,Eq)
instance Prediction SignalPResult where
 contgName (Sp name _ _) = name

data TMHMMResult = TMHMM String Int deriving (Show,Eq)
instance Prediction SignalPResult where
 contgName (Sp name _ _) = name


class Prediction a where
 contgName :: a -> String

ParsePred lineToPred cutoffFun rawParse file = ((map lineToPred) . cutoffFun . (map words) . rawParse) <$> readFile file

-----------------------------------------------------
--TargetP parsing
rawTP = tail . init . init . (drop 8) . lines


parseTP file = parsePred toPredictionT chooseSuccess rawTP file
 where chooseSuccess = filter (\x -> x!!5 == "S")
       toPredictionT [name,len,mtp,sp,other,loc,rc] = Tp {contgName = name, qual = read rc}
------------------------------------------------------
--SignalP parsing
rawSP = (drop 2) . lines
parseSP file = parsePred toPredS chooseSuccess rawSP file
 where chooseSuccess = filter (\x -> x!!13 == "Y" && (last x) == "Y")
       toPredS xs = Sp {contgName = (head xs), nnQual = read $ xs!!12, hmmQual = read $ xs!!19}
------------------------
--TMHMM parsing
rawTMHMM = (filter (\x -> elem "Number" $ words x)) . lines
parseTMHMM file = parsePred toPredTMHMM chooseSuccess rawTMHMM file
 where chooseSuccess = filter (\x -> (read . last $ x) > 0)
       toPredTMHMM xs = TMHMM {contgName = xs!!1, qual = read . last $ xs}
--------------------------------
--merge predictions

groupByName = (groupBy (\a b -> (contgName a) == (contgName b))) . (sortBy (\a b -> compare (contgName a) (contgName b)))




--Formatting
sigPLimit = 2000

formatForPredictors fileName = do
 seqs <- (fitLines (sigPLimit `div` 6)) <$> readFasta fileName
 let fileNames = map (\indx -> newExt fileName (('-':(show indx))++".prots")) [1..length seqs]
 zipWithM_ writeTranslate fileNames seqs
 return fileNames


breakLines n str
 | length str >= n = Just (splitAt n str)
 | length str > 0 = Just (str,[])
 | otherwise = Nothing

fitLines n = unfoldr $ breakLines n
--------------------------
--FIN
