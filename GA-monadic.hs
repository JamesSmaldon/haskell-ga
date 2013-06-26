import Control.Monad
import Control.Monad.State
import System.Random
import Data.Char
import Data.List

data GenotypeInfo = GenotypeInfo { minSize :: Int, 
                                   maxSize :: Int } deriving (Show)

class Individual a where
    mutate :: (RandomGen g) => PopulationInfo -> a -> State g a
    crossover :: (RandomGen g) => a -> a -> State g (a, a)
    createRandom :: (RandomGen g) => GenotypeInfo -> State g a

type Fitness = Float

class Phenotype a where
    assignFitness :: (Individual a) => a -> Fitness

data Population a = Population [(a, Fitness)] deriving (Show)

data PopulationInfo = PopulationInfo { size :: Int,
                                       mutProb :: Float,
                                       tournamentSize :: Int }

data StrGene = StrGene String deriving (Show)

randomVal :: (RandomGen g, Random a) => State g a
randomVal = do
                x <- get
                let (newC, x') = random x
                put x'
                return newC

randomRVal :: (RandomGen g, Random a) => (a,a) -> State g a
randomRVal rng = do
                        x <- get
                        let (newC, x') = randomR rng x
                        put x'
                        return newC

dropVal :: Int -> [a] -> (a, [a])
dropVal i l = (l !! i, xs ++ tail ys)
                where (xs, ys) = splitAt i l

dropRand :: (RandomGen g) => [a] -> State g (a, [a])
dropRand l = randomRVal(0, length l - 1) >>= return . (flip dropVal $ l)

chooseRand :: (RandomGen g) => [a] -> State g a
chooseRand l = randomRVal(0, length l - 1) >>= return . ((!!) l)


instance Functor Population where
    fmap f (Population inds) = Population (map g inds)
        where g = \(ind,fit) -> (f ind, fit)

                    
createRandomPopulation :: (RandomGen g, Phenotype a, Individual a) => PopulationInfo -> GenotypeInfo -> State g (Population a)
createRandomPopulation pinfo ginfo = newPop >>= return . Population . map (\i -> (i, assignFitness i))
                                        where newPop = (replicateM (size pinfo) (createRandom ginfo))

cmpFitness :: (Individual a) => (a, Fitness) -> (a, Fitness) -> Ordering
cmpFitness (_, f1) (_, f2) = compare f1 f2 

tournamentSelect :: (RandomGen g, Individual a) => Int -> Population a -> State g a
tournamentSelect count (Population inds) = replicateM count (chooseRand inds) >>= (return . fst . head . (sortBy cmpFitness))

mate :: (RandomGen g, Individual a) => PopulationInfo -> Population a -> State g (a, a)
mate popinfo pop = do
                        p1 <- tournamentSelect tsize pop
                        p2 <- tournamentSelect tsize pop
                        (c1, c2) <- crossover p1 p2
                        liftM2 (,) (mutate popinfo c1) (mutate popinfo c2)
                        where tsize = tournamentSize popinfo

evolve :: (RandomGen g, Individual a, Phenotype a) => PopulationInfo -> Population a -> State g (Population a)
evolve popinfo (Population inds) = replicateM numMates (mate popinfo (Population inds)) >>= return . Population . (foldr ttoa []) 
                                    where numMates = length inds `div` 2
                                          ttoa = (\(f,s) a -> (f, assignFitness f) : (s, assignFitness s) : a)

generateRun :: (RandomGen g, Individual a, Phenotype a) => Int -> PopulationInfo -> (Population a -> State g (Population a))
generateRun 0 popinfo = evolve popinfo 
generateRun n popinfo = (evolve popinfo) >=> generateRun (n-1) popinfo 

optimise :: (Individual a, Phenotype a) => Int -> PopulationInfo -> GenotypeInfo -> Population a
optimise iters popinfo ginfo = evalState (generateRun iters popinfo pop) r
                                    where (pop, r) = runState (createRandomPopulation popinfo ginfo) (mkStdGen 1)

mutateChar :: (RandomGen g) => PopulationInfo -> Char -> State g Char
mutateChar pinfo c = do
                        decider <- randomVal
                        newC <- randomRVal (0,255)
                        if decider > (mutProb pinfo) then return c else return $ chr newC

instance Individual StrGene where
    mutate pinfo (StrGene str) = liftM StrGene (mapM (mutateChar pinfo) str)
    crossover (StrGene str1) (StrGene str2) = do
                                                break <- randomRVal(0, min (length str1 - 1) (length str2 - 1))
                                                let i1 = ((take break str1) ++ (drop break str2))
                                                let i2 = ((take break str2) ++ (drop break str1))
                                                return $ (StrGene i1, StrGene i2)
                                                      
    createRandom ginf = do
                            size <- randomRVal (minSize ginf, maxSize ginf)
                            nums <- replicateM size (randomRVal (0,255))
                            return $ StrGene (map chr nums)


instance Phenotype StrGene where
    assignFitness (StrGene str) = (fromIntegral . sum) (zipWith charCmp "Hello World!" str)
                                                where charCmp x y = abs (ord x - ord y)
                    
{-
Example run

x = PopulationInfo 100 0.03 3
y = GenotypeInfo 12 12
z = optimise 1 x y :: Population StrGene 
-}
