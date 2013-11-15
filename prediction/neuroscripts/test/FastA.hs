module FastA (readFasta, grokFasta, convertFasta, translateAllFrames, writeTranslate, Contig, reading, header, Translate) where
{-# LANGUAGE NoMonomorphismRestriction #-} -- Don't let The Man keep you down
import Text.Parsec
import Control.Monad
import Control.Monad.Trans
import Data.List
import Control.Applicative ((<$>),(<*>))
import System.Environment
import System.IO
---------------------------------------
-----------------------------------------------------
data Contig = Contig {header :: String, reading :: String} deriving (Show,Eq)
data Translate = Translate {protHeader :: String, protReading :: String} deriving (Show,Eq)
type ParsecIO a = ParsecT String () IO a
--------------------------------------------------------
maxLength = 60

rawParse prsr str = unEither <$> runParserT prsr () "" str
unEither a = case a of
 Right b -> b
 Left _ -> error "failed parse"
getHeader :: ParsecIO String
getHeader = (char '>') >> (many alphaNum)
getSeq :: ParsecIO String
getSeq = concat `liftM` endBy (many letter) (char '\n')
getContig :: ParsecIO Contig
getContig = liftM2 Contig getHeader getSeq
getReading :: ParsecIO [Contig]
getReading = many getContig
----------------------------------------------------
translateByLongestFrame :: Contig -> Translate
translateByLongestFrame contg = Translate (header contg) (concat . (map $ assureLookup trans_tbl) . longest $ map (toCodon . (flip ($) (reading contg))) [id,drop 1, drop 2, reverse, (drop 1) . reverse, (drop 2) . reverse])
 where longest = foldl' (\a x -> if length a > length x then a else x) []

toCodon = unfoldr (\str -> if length str >= 3 then Just (splitAt 3 str) else Nothing)

translateAllFrames :: Contig -> [Translate]
translateAllFrames contg = zipWith Translate (headers . header $ contg) (((map (assureLookup trans_tbl <=< toCodon)) . frames . reading) contg)
 where headers str = (map (flip (++)) ["Xf1","Xf2","Xf3","Xr1","Xr2","Xr3"]) <*> [str]
       frames str = [id,drop 1,drop 2, reverse, (drop 1) . reverse, (drop 2) . reverse] <*> [str]


assureLookup tbl key = case (lookup key tbl) of
 Nothing -> "X"
 Just amino -> amino

trans_tbl :: [(String,String)]
trans_tbl = [("GCT","A"),("GCC","A"),("GCA","A"),("GCG","A"),
             ("CGT","R"),("CGC","R"),("CGA","R"),("CGG","R"),("AGA","R"),("AGG","R"),
             ("AAT","N"),("AAC","N"),
             ("GAT","D"),("GAC","D"),
             ("TGT","C"),("TGC","C"),
             ("CAA","Q"),("CAG","Q"),
             ("GAA","E"),("GAG","E"),
             ("GGT","G"),("GGC","G"),("GGA","G"),("GGG","G"),
             ("CAT","H"),("CAC","H"),
             ("ATT","I"),("ATC","I"),("ATA","I"),
             ("TTA","L"),("TTG","L"),("CTT","L"),("CTC","L"),("CTA","L"),("CTG","L"),
             ("AAA","K"),("AAG","K"),
             ("TTT","F"),("TTC","F"),
             ("CCT","P"),("CCC","P"),("CCA","P"),("CCG","P"),
             ("TCT","S"),("TCC","S"),("TCA","S"),("TCG","S"),("AGT","S"),("AGC","S"),
             ("ACT","T"),("ACC","T"),("ACA","T"),("ACG","T"),
             ("TAT","Y"),("TAC","Y"),
             ("TGG","W"),
             ("GTT","V"),("GTC","V"),("GTA","V"),("GTG","V")]
stopCodons = ["TAG","TGA","TAA"]
startCodon = "ATG"


{-this version enforces start/stop codons
toCodon contg = (drop 1) . (dropWhile (/= startCodon)) $ unfoldr readCodon contg
readCodon str
 | length str >= 3 = if take 3 str `elem` stopCodons then Nothing else Just (splitAt 3 str)
 | otherwise = Nothing
-}




breakLines n str 
 | length str >= n = Just (splitAt n str)
 | length str > 0 = Just (str,"")
 | otherwise = Nothing

fitLines n = unfoldr $ breakLines n

formatFasta p = unlines $ ('>':(protHeader p)) : (fitLines maxLength $ protReading p)

writeTranslate name seqs = writeFile name (seqs >>= (translateAllFrames >=> formatFasta))
grokFasta = rawParse getReading
readFasta = grokFasta <=< readFile

convertFasta name = do
 let fileName = (++ ".prots") . (takeWhile (/= '.')) $ name
 seqs <- grokFasta <=< readFile $ name
 writeTranslate fileName seqs
