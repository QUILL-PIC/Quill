module Chameleon (configure, get, Particle(..), readParticles) where

{-# LANGUAGE ForeignFunctionInterface #-}

import Foreign.C.Types
import Foreign.C.String
import Foreign.Marshal.Array
import Foreign.Ptr
import Foreign.Storable
-- import Foreign.ForeignPtr -- https://hackage.haskell.org/package/base-4.9.0.0/docs/Foreign-ForeignPtr.html

foreign import ccall "configure" cConfigure :: CString -> Bool -> IO ()

configure :: [Char] -> Bool -> IO ()
configure logName shouldPrint = withCString logName (\x -> cConfigure x shouldPrint)

foreign import ccall "get" cGet :: CString -> IO CDouble

get :: [Char] -> IO Double
get varName = withCString varName (\s -> cGet s >>= return . realToFrac)
-- if cGet :: CString -> CDOUble (i.e., without IO) and get varName = withCString varName (return .
-- realToFrac . cGet) -- then bug with (sometimes) wrong varName passed arises

data Particle = Particle {
    q :: Double,
    x :: Double,
    y :: Double,
    z :: Double,
    ux :: Double,
    uy :: Double,
    uz :: Double,
    g :: Double,
    chi :: Double
    } deriving(Show)

-- C structure can be treated as "array" somehow, with 8 byte alignment in mind.
-- https://wiki.haskell.org/FFI_complete_examples#Reading_structs_in_Haskell
foreign import ccall "read_particles" cReadParticles :: CString -> Ptr CParticles

data CParticles = CParticles CInt (Ptr Particle)

instance Storable CParticles where
    sizeOf _ = 16
    alignment _ = 8
    peek ptr = peekByteOff ptr 0 >>=
               \a -> peekByteOff ptr 8 >>=
               \b -> return $ CParticles a b

instance Storable Particle where
    sizeOf _ = 8 * 9
    alignment _ = 8
    peek ptr = peekByteOff ptr 0 >>=
               \q -> peekByteOff ptr 8 >>=
               \x -> peekByteOff ptr 16 >>=
               \y -> peekByteOff ptr 24 >>=
               \z -> peekByteOff ptr 32 >>=
               \ux -> peekByteOff ptr 40 >>=
               \uy -> peekByteOff ptr 48 >>=
               \uz -> peekByteOff ptr 56 >>=
               \g -> peekByteOff ptr 64 >>=
               \chi -> return $ Particle {q = q, x = x, y = y, z = z, ux = ux, uy = uy, uz = uz, g = g, chi = chi}

readParticles :: [Char] -> IO [Particle]
readParticles fileName = toParticles $ withCString fileName (return . cReadParticles)
    where toParticles :: IO (Ptr CParticles) -> IO [Particle]
          toParticles ptr = ptr >>= peek >>=
              \(CParticles n ptr') -> peekArray (fromEnum n) ptr'
