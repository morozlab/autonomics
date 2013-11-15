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
 sequence_ $ [runSignalP,runTargetP,runTMHMM,mergeAllPreds parses] <*> [fileNames]
parses = [parseSP,parseTP,parseTMHMM]
 
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
  where tmhmm file = "cat " ++ file ++ " | decodeanhmm -f " ++ predDir ++ "/tmhmm-2.0c/lib/TMHMM2.0.options -modelfile " ++ predDir ++ "/tmhmm-2.0c/lib/TMHMM2.0.model | " ++ predDir ++ "/tmhmm-2.0c/bin/tmhmmformat.pl > " ++ (newExt file "_tmhmm.out")
------------------------------------------------------

--Parsing
data Prediction = Tp {contgName :: String, qual :: Int} | Sp {contgName :: String, nnQual :: Float, hmmQual :: Float} | TMHMM {contgName :: String, qual :: Int} deriving (Show,Eq)

parsePred :: ([String] -> Prediction) -> ([[String]] -> [[String]]) -> (String -> [String]) -> FilePath -> IO [Prediction]
parsePred lineToPred cutoffFun rawParse file = ((map lineToPred) . cutoffFun . (map words) . rawParse) <$> readFile file
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
 where chooseSuccess = filter (\x -> (read . last $ x) == 0)
       toPredTMHMM xs = TMHMM {contgName = xs!!1, qual = read . last $ xs}
--------------------------------
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
--merger
merge predGroups = (filter ((==3) . length)) . (map sortPreds) . groupByName $ concat predGroups
groupByName = (groupBy (\a b -> (contgName a) == (contgName b))) . (sortBy (\a b -> compare (contgName a) (contgName b)))
sortPreds preds = sortBy (\a b -> compare (orderPred a) (orderPred b)) $ preds
orderPred pred = case pred of
 Tp _ _ -> 1 
 Sp _ _ _ -> 2
 TMHMM _ _ -> 3 

comparePreds pred1 pred2 = orderPred pred1 == orderPred pred2

unBy chr ws = foldr1 (\w s -> w ++ chr:s) ws
mergeFiles a b c = do
 x <- parseTP a
 y <- parseSP b
 z <- parseTMHMM c
 return $ merge [x,y,z]


formatPreds preds = tabDelim lined
 where lined = map (\pred -> (contgName . head $ pred):(map formatPred pred)) $ preds
       formatPred pred = case pred of
        Tp _ q -> show q
        Sp _ q1 q2 -> (show q1) ++ ('\t':(show q2))
        TMHMM _ q -> show q

tabDelim = (unBy '\n') . (map $ unBy '\t')
       
mergeAllPreds prsrs fileNames = do
 allPredicts <- prsrs <*> fileNames
 return $ groupBy comparePreds allPredicts
 (groupBy comparePreds =<<) <$> (prsrs <*> fileNames)

--FIN
