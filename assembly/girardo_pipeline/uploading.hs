import HSH
import System.Directory
import System.Environment
main = do
 [signal,projName] <- getArgs
 case signal of
  "-s" -> sendFun projName
  "-c" -> cleanupFun projName

sendFun projName = do
 setCurrentDirectory $ "/srv/data2/pipeline/" ++ projName ++ "/"
 files <- map (\x -> "/srv/data2/pipeline/" ++ projName ++ "/" ++ x) `fmap` (run "ls" :: IO [String])
 writeFile "listFile.txt" (unlines $ projName:files)
 runIO "curl -T listFile.txt ftp://ftpguest:7002ram@74.252.103.104"

cleanupFun projName = do
 setCurrentDirectory $ "/srv/data2/pipeline/" ++ projName ++ "/"
 writeFile "listFile.txt" ""
 runIO "curl -T listFile.txt ftp://ftpguest:7002ram@74.252.103.104"
