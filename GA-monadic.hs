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

class Phenotype a where
    assignFitness :: (Individual a) => a -> Float

data Population a = Population [a] deriving (Show)

data PopulationInfo = PopulationInfo { size :: Int,
                                       mutProb :: Float }

data GAInfo = GAInfo { generations :: Int }

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

mutateChar :: (RandomGen g) => PopulationInfo -> Char -> State g Char
mutateChar pinfo c = do
                        decider <- randomVal
                        newC <- randomRVal (0,255)
                        if decider > (mutProb pinfo) then return c else return $ chr newC

instance Functor Population where
    fmap f (Population inds) = Population (map f inds)

                    
createRandomPopulation :: (RandomGen g, Individual a) => PopulationInfo -> GenotypeInfo -> State g (Population a)
createRandomPopulation pinfo ginfo = liftM Population (replicateM (size pinfo) (createRandom ginfo))

tournamentSelect :: (RandomGen g, Individual a, Ord a) => Int -> Population a -> State g a
tournamentSelect count (Population inds) = replicateM count (chooseRand inds) >>= (return . head . sort)

mate :: (RandomGen g, Individual a, Ord a) => PopulationInfo -> Population a -> State g (a, a)
mate popinfo pop = do
                        p1 <- tournamentSelect 3 pop
                        p2 <- tournamentSelect 3 pop
                        (c1, c2) <- crossover p1 p2
                        liftM2 (,) (mutate popinfo c1) (mutate popinfo c2)

evolve :: (RandomGen g, Individual a, Ord a) => PopulationInfo -> Population a -> State g (Population a)
evolve popinfo (Population inds) = replicateM numMates (mate popinfo (Population inds)) >>= return . Population . (foldr ttoa []) 
                                    where numMates = length inds `div` 2
                                          ttoa = (\(f,s) a -> f : s : a)

generateRun :: (RandomGen g, Individual a, Ord a) => Int -> PopulationInfo -> (Population a -> State g (Population a))
generateRun 0 popinfo = evolve popinfo 
generateRun n popinfo = (evolve popinfo) >=> generateRun (n-1) popinfo 

optimise :: (Individual a, Ord a) => Int -> PopulationInfo -> GenotypeInfo -> Population a
optimise iters popinfo ginfo = evalState (generateRun iters popinfo pop) r
                                    where (pop, r) = runState (createRandomPopulation popinfo ginfo) (mkStdGen 1)

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
                    
instance Eq StrGene where
    (==) lhs rhs = assignFitness lhs == assignFitness rhs

instance Ord StrGene where
    (<=) lhs rhs = assignFitness lhs <= assignFitness rhs

x = PopulationInfo 100 0.03
y = GenotypeInfo 12 12
z = optimise 1 x y :: Population StrGene

get2 :: Population StrGene -> (StrGene, StrGene)
get2 (Population inds) = (head inds , head . tail $ inds) 
