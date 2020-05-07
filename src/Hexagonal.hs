-- |
-- Copyright: (c) 2020 tim put
-- SPDX-License-Identifier: GPL-3.0-only
-- Maintainer: tim put <timput@gmail.com>
--
-- Utilites for data sampled on regular hexagonal grids.
module Hexagonal where

import Linear.Matrix
import Linear.Metric
import Linear.V2
import Linear.Vector

-- /‾\
-- \_/
-- Basis for 'wide' hexagons
hBasis :: (Floating a) => M22 a
hBasis = V2 (V2 (sqrt 3 / 2) 0) (V2 (0.5) 1)

-- /‾\
-- \_/
-- Inverse of basis for 'wide' hexagons
hBasis' :: (Floating a) => M22 a
hBasis' = inv22 $ hBasis

--- |‾‾|
--- |__|
--- Square basis
sBasis :: (Floating a) => M22 a
sBasis = V2 (V2 1 0) (V2 0 1)

-- Convert square coordinates to axis aligned hex coordinates
sToH :: Floating a => V2 a -> V2 a
sToH v = hBasis' !* v

-- Convert axis aligned hex coordinates to square coordinates
hToS :: Floating a => V2 a -> V2 a
hToS v = hBasis !* v

-- convert offset row hex coordinates to square coordinates
rhToS :: (Integral a, Floating b) => V2 a -> V2 b
rhToS (V2 x y) = V2 (fromIntegral x + if even y then 0 else 0.5) (fromIntegral y * bh)
  where
    bh = sqrt 3 / 2

-- convert square coordinates to offset row hex coordinates
-- since offset row hex coordinates are only well defined on integers
-- this function is necessarily badly behaved. The intended use is computing
-- the size of hex which will fit in a given rectangle.
sToRh :: (Integral a1, Floating a2) => V2 a1 -> V2 a2
sToRh (V2 x y) = V2 (fromIntegral x - if even y then 0 else 0.5) (fromIntegral y * bh)
  where
    bh = 2 / sqrt 3

-- convert offset row hexagonal coordinates to axis aligned coordinates
rhToH :: Integral a => V2 a -> V2 a
rhToH (V2 x y) = V2 (x - y `div` 2) y

-- convert offset row hexagonal coordinates to axis aligned coordinates
hToRh :: Integral a => V2 a -> V2 a
hToRh (V2 x y) = V2 (x + y `div` 2) y

--- a       b
---  |‾‾‾‾‾|
---  |  v  |
---  |_____|
--- c       d
---
--- Computes the bilinear interpolation of a point v in a square of values a,b,c,d
bilinear :: Num a => V2 a -> a -> a -> a -> a -> a
bilinear (V2 x y) a b c d = (V2 (1 - y) y) `dot` (V2 (V2 a b) (V2 c d) !* (V2 (1 - x) x))

-- sample rectangular space with axis aligned hexagonal coords using bilinear interpolation
bilinearStoH :: (V2 Int -> Float) -> V2 Int -> Float
bilinearStoH arr hex =
  let v = hToS . fmap fromIntegral $ hex
      V2 i j = fmap floor v
   in bilinear
        v
        (arr (V2 i j))
        (arr (V2 (i + 1) j))
        (arr (V2 i (j + 1)))
        (arr (V2 (i + 1) (j + 1)))

-- sample rectangular space with offset row hexagonal coords using bilinear interpolation
bilinearStoRh ::
  (Floating a1, RealFrac a1, Integral a3, Integral a2) =>
  (V2 a2 -> a1) ->
  V2 a3 ->
  a1
bilinearStoRh arr hex =
  let v = rhToS $ hex
      V2 i j = fmap floor v
   in bilinear
        (v - fmap fromIntegral (V2 i j))
        (arr (V2 i j))
        (arr (V2 (i + 1) j))
        (arr (V2 i (j + 1)))
        (arr (V2 (i + 1) (j + 1)))

-- generate list of indices of axis aligned regular hexagon of side
-- length l, offset from the origin by offset, in counter clockwise
-- order with local origin at the lower left corner.
hexagonIndices :: (Integral a1, Num a2) => a1 -> V2 a2 -> [V2 a2]
hexagonIndices l offset =
  for [0 .. l -1] (\i' -> let i = fromIntegral i' in offset + i *^ b1) -- right
    ++ for [0 .. l -1] (\i' -> let i = fromIntegral i' in offset + i *^ b2 + (l') *^ b1) -- up-right
    ++ for [0 .. l -1] (\i' -> let i = fromIntegral i' in offset + i *^ (b2 - b1) + (l') *^ b1 + (l') *^ b2) -- up-left
    ++ for [0 .. l -1] (\i' -> let i = fromIntegral i' in offset - i *^ b1 + (l') *^ b1 + (l') *^ b2 + (l') *^ (b2 - b1)) -- left
    ++ for [0 .. l -1] (\i' -> let i = fromIntegral i' in offset - i *^ b2 + (l') *^ b2 + (l') *^ (b2 - b1)) -- down-left
    ++ for [0 .. l -1] (\i' -> let i = fromIntegral i' in offset - i *^ (b2 - b1) + (l') *^ (b2 - b1)) -- down-right
    ++ [offset] -- close the loop
  where
    l' = fromIntegral l
    -- b1 = V2 1 0
    -- b2 = V2 (sqrt 3 / 2) 0.5
    b1 = V2 1 0
    b2 = V2 0 1
    for = flip map

-- In row offset coordinates
hexagonIndicesRh ::
  (Integral a1, Integral a2) => a2 -> V2 a1 -> [V2 a1]
hexagonIndicesRh l offset = fmap hToRh $ hexagonIndices l offset

-- determine if a point, given in cartesian coordinates lies inside of
-- a regular hexagon of side length l with local origin offset.
pointInHex :: (Ord a, Floating a) => a -> V2 a -> V2 a -> Bool
pointInHex l offset p = aboveBottom && belowTop
  where
    V2 x y = ((p - offset) ^/ l)
    aboveBottom = y >= max (sqrt 3 * (abs (x -0.5) - 1 / 2)) 0
    belowTop = y <= min (negate (sqrt 3 * (abs (x -0.5) + 1 / 2)) + (2 * sqrt 3)) (sqrt 3)
