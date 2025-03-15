mod add;
mod div;
mod gcd;
mod matrix;
mod mul;
mod shift;
#[cfg(test)]
mod test_util;
use std::{fmt::Display, str::FromStr};

use smallvec::SmallVec;

pub type Limb = std::ffi::c_ulong;
type LimbStorage = SmallVec<[Limb; 2]>;

const BITS: usize = Limb::BITS as usize;

const HALF_MASK: Limb = interleave_mask(BITS.ilog2() - 1);

fn normalized_len(limbs: &[Limb]) -> usize {
    limbs.len() - limbs.iter().rev().take_while(|&&x| x == 0).count()
}

fn deg(limbs: &[Limb]) -> u64 {
    let Some(last) = limbs.last() else {
        return 0;
    };
    limbs.len() as u64 * BITS as u64 - 1 - last.leading_zeros() as u64
}

const fn limbs_for_deg(n: u64) -> usize {
    (n as u64 / BITS as u64) as usize + 1
}

#[derive(Clone, Default, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub struct Gf2Poly {
    deg: u64,
    limbs: LimbStorage,
}

impl Gf2Poly {
    fn resize_to_deg_unnormalized(&mut self, deg: u64) {
        self.deg = deg;
        self.limbs.resize(limbs_for_deg(deg), 0);
    }

    /// Degree of the polynomial. Returns degree 0 for the zero polynomial.
    pub fn deg(&self) -> u64 {
        self.deg
    }

    /// Whether the polynomial is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.limbs.is_empty()
    }

    /// Gets the nth coefficient of the polynomial,
    /// with the 0th coefficient being the constant part.
    pub fn get(&self, n: u64) -> bool {
        let Ok(idx) = usize::try_from(n / BITS as u64) else {
            return false;
        };
        match self.limbs.get(idx) {
            Some(limb) => {
                let res = n % BITS as u64;
                (limb >> res) & 1 == 1
            }
            None => false,
        }
    }

    /// Sets the nth coefficient of the polynomial to 1,
    /// with the 0th coefficient being the constant part.
    pub fn set(&mut self, n: u64) {
        if n > self.deg || self.is_zero() {
            self.resize_to_deg_unnormalized(n);
            self.deg = n;
        }
        let idx = n / BITS as u64;
        let res = n % BITS as u64;
        self.limbs[usize::try_from(idx).unwrap()] |= 1 << res;
    }

    /// Converts a vector of limbs into a polynomial, with least significant limb first.
    fn from_limb_storage(mut limbs: LimbStorage) -> Self {
        limbs.truncate(normalized_len(&limbs));
        let deg = deg(&limbs);
        Gf2Poly { limbs, deg }
    }

    /// Converts a slice of limbs into a polynomial, with least significant limb first.
    pub fn from_limbs(limbs: &[Limb]) -> Self {
        let len = normalized_len(&limbs);
        let limb_vec = LimbStorage::from(&limbs[..len]);
        let deg = deg(&limbs);
        Gf2Poly {
            limbs: limb_vec,
            deg,
        }
    }

    /// Converts a slice of bytes into a polynomial, with least significant byte first.
    pub fn from_bytes(bytes: &[u8]) -> Gf2Poly {
        let limb_size = std::mem::size_of::<Limb>();
        let mut limbs = LimbStorage::with_capacity(bytes.len() / limb_size + 1);
        for chunk in bytes.chunks(std::mem::size_of::<Limb>()) {
            let mut limb = 0;
            for (i, &byte) in chunk.iter().enumerate() {
                limb |= (byte as Limb) << (i << 3);
            }
            limbs.push(limb);
        }
        Gf2Poly::from_limb_storage(limbs)
    }

    /// 0 (Zero).
    /// The polynomial with all coefficients 0.
    /// The unique polynomial that is divisible by every other polynomial.
    /// The unique fixpoint of the derivative operator on GF(2)[x].
    pub fn zero() -> Self {
        Default::default()
    }

    /// 1
    pub fn one() -> Self {
        Gf2Poly {
            deg: 0,
            limbs: std::iter::once(1).collect(),
        }
    }

    /// Returns the x of GF(2)[x].
    pub fn x() -> Self {
        Gf2Poly {
            deg: 1,
            limbs: std::iter::once(2).collect(),
        }
    }

    /// Calculates the derivative of self.
    pub fn derivative(&self) -> Self {
        const MASK: Limb = Limb::MAX / 3;
        let mut limbs: LimbStorage = self.limbs.iter().map(|x| x >> 1 & MASK).collect();
        limbs.truncate(normalized_len(&limbs));
        let deg = deg(&limbs);
        Gf2Poly { deg, limbs }
    }

    /// Calculates the number of coefficients with value 1.
    pub fn hamming_weight(&self) -> u64 {
        self.limbs.iter().map(|x| x.count_ones() as u64).sum()
    }

    /// Calculates the polynomial with all coefficients reversed in place.
    /// If f is of degree n, this is x^n * f(1/x).
    pub fn reverse(&self) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        let mut limbs;
        let bit_offset = (self.deg + 1) % BITS as u64;
        if bit_offset == 0 {
            limbs = self
                .limbs
                .iter()
                .rev()
                .copied()
                .map(Limb::reverse_bits)
                .collect();
        } else {
            limbs = LimbStorage::with_capacity(self.limbs.len());
            let mut prev = 0;
            let mut first = true;
            for x in self.limbs.iter().rev() {
                let rev = x.reverse_bits();
                if !first {
                    limbs.push(rev << bit_offset | prev);
                } else {
                    first = false;
                }
                prev = rev >> (BITS as u64 - bit_offset);
            }
            if prev != 0 {
                limbs.push(prev);
            }
        }
        limbs.truncate(normalized_len(&limbs));
        let deg = deg(&limbs);
        Gf2Poly { deg, limbs }
    }

    /// Calculates n for f such that x^n || f, i.e. the number of trailing zeros.
    /// Returns None if f is zero.
    pub fn trailing_zeros(&self) -> Option<u64> {
        if self.is_zero() {
            return None;
        }
        let mut zeros = 0;
        for limb in self.limbs.iter() {
            if *limb == 0 {
                zeros += BITS as u64;
            } else {
                zeros += limb.trailing_zeros() as u64;
                break;
            }
        }
        Some(zeros)
    }

    /// Calculates f mod x^n, i.e. truncates f to the n least significant bits.
    pub fn truncate_mut(&mut self, n: u64) {
        if n == 0 {
            *self = Self::zero();
            return;
        }
        if self.deg < n {
            return;
        }
        let n_limbs = (n - 1) as usize / BITS;
        let n_bits = (n - 1) as usize % BITS;
        self.limbs.truncate(n_limbs + 1);
        self.limbs[n_limbs] &= (1 << n_bits) | (1 << n_bits) - 1;
        self.limbs.truncate(normalized_len(&self.limbs));
        self.deg = deg(&self.limbs);
    }

    /// Calculates f mod x^n, i.e. truncates f to the n least significant bits.
    pub fn truncated(&self, n: u64) -> Self {
        if n == 0 {
            return Self::zero();
        }
        let mut cpy = self.clone();
        cpy.truncate_mut(n);
        cpy
    }

    /// Given f, calculates f * f in linear time.
    pub fn square(&self) -> Self {
        if self.is_zero() {
            return Self::zero();
        }

        let deg = self.deg * 2;
        let mut limbs = LimbStorage::with_capacity(limbs_for_deg(deg));
        let (last_limb, start_limbs) = self.limbs.split_last().unwrap();
        let lo_half = |x: Limb| x & HALF_MASK;
        let hi_half = |x: Limb| (x >> BITS / 2) & HALF_MASK;
        for x in start_limbs {
            limbs.push(spacen_bits(lo_half(*x)));
            limbs.push(spacen_bits(hi_half(*x)));
        }
        limbs.push(spacen_bits(lo_half(*last_limb)));
        let hi = hi_half(*last_limb);
        if hi != 0 {
            limbs.push(spacen_bits(hi));
        }
        Gf2Poly { deg, limbs }
    }

    /// Evaluates the polynomial at the point x.
    pub fn eval(&self, x: bool) -> bool {
        if !x {
            1 & *self.limbs.first().unwrap_or(&0) != 0
        } else {
            1 & self.hamming_weight() != 0
        }
    }

    /// Splits the polynomial into hi * x^n + lo where lo.deg() < n.
    /// Returns (hi, lo).
    pub fn split_at(&self, n: u64) -> (Self, Self) {
        if self.is_zero() {
            return (Self::zero(), Self::zero());
        }
        let hi = self.clone() >> n;
        let lo = self.truncated(n);
        (hi, lo)
    }

    #[cfg(test)]
    fn is_normalized(&self) -> bool {
        self.limbs.is_empty() && self.deg == 0 || self.limbs.len() == limbs_for_deg(self.deg)
    }
}

const fn interleave_mask(n: u32) -> Limb {
    !0 / ((1 << (1 << n)) + 1)
}

fn spacen_bits(mut x: Limb) -> Limb {
    for i in (0..(Limb::BITS.ilog2() - 1)).rev() {
        x = (x | (x << (1 << i))) & interleave_mask(i);
    }
    x
}

impl Display for Gf2Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Some((last, limbs)) = self.limbs.split_last() else {
            write!(f, "0")?;
            return Ok(());
        };
        write!(f, "{last:x}")?;
        for limb in limbs.iter().rev() {
            write!(f, "{limb:016x}")?;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum Gf2PolyConversionError {
    UnexpectedChar(char),
}

impl FromStr for Gf2Poly {
    type Err = Gf2PolyConversionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut limbs = LimbStorage::with_capacity(s.len() / BITS + 1);
        let mut limb = 0;
        let mut limb_idx = 0;
        for c in s.chars().rev() {
            let Some(digit) = c.to_digit(16) else {
                return Err(Gf2PolyConversionError::UnexpectedChar(c));
            };
            limb |= (digit as Limb) << (limb_idx % BITS);
            limb_idx += 4;
            if limb_idx % BITS == 0 {
                limbs.push(limb);
                limb = 0;
            }
        }
        if limb_idx % BITS != 0 {
            limbs.push(limb);
        }
        Ok(Gf2Poly::from_limb_storage(limbs))
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn deriv_linear(a: Gf2Poly, b: Gf2Poly) {
            prop_assert_poly_eq!((&a + &b).derivative(), a.derivative() + b.derivative());
        }

        #[test]
        fn deriv_mul(a: Gf2Poly, b: Gf2Poly) {
            prop_assert_eq!((&a * &b).derivative(), &a.derivative() * &b + a * b.derivative());
        }

        #[test]
        fn deriv_twice(a: Gf2Poly) {
            prop_assert_poly_eq!(a.derivative().derivative(), Gf2Poly::zero());
        }

        #[test]
        fn string_roundtrip(a: Gf2Poly) {
            let s = a.to_string();
            let b: Gf2Poly = s.parse().unwrap();
            prop_assert_poly_eq!(a, b);
        }

        #[test]
        fn hamming_metric(a: Gf2Poly, b: Gf2Poly) {
            prop_assert!(a.hamming_weight() + b.hamming_weight() >= (a + b).hamming_weight());
        }

        #[test]
        fn shift_euclid(mut a: Gf2Poly, split in 0u64..128) {
            let hi = a.clone() >> split << split;
            let lo = a.truncated(split);
            prop_assert_poly_eq!(a, hi + lo);
        }

        #[test]
        fn reverse_invariant(a: Gf2Poly) {
            let roundtrip = a.reverse().reverse() << a.trailing_zeros().unwrap_or(0);
            prop_assert_poly_eq!(roundtrip, a);
        }

        #[test]
        fn reverse_homo(a: Gf2Poly, b: Gf2Poly) {
            prop_assert_poly_eq!(a.reverse() * b.reverse(), (a * b).reverse());
        }

        #[test]
        fn square_is_square(a: Gf2Poly) {
            prop_assert_poly_eq!(&a * &a, a.square())
        }

        #[test]
        fn eval_mul(a: Gf2Poly, b: Gf2Poly) {
            for x in [false, true] {
                prop_assert_eq!(a.eval(x) & b.eval(x), (&a * &b).eval(x))
            }
        }

        #[test]
        fn eval_add(a: Gf2Poly, b: Gf2Poly) {
            for x in [false, true] {
                prop_assert_eq!(a.eval(x) ^ b.eval(x), (&a + &b).eval(x))
            }
        }

        #[test]
        fn get_set(a: Gf2Poly) {
            let mut b = Gf2Poly::zero();
            for i in 0..=a.deg() {
                if a.get(i) {
                    b.set(i);
                }
            }
            prop_assert_poly_eq!(a, b)
        }
    }

    #[test]
    fn derivative() {
        let a: Gf2Poly = "dbf".parse().unwrap();
        let b = a.derivative();
        assert_eq!(b.to_string(), "455");
        assert_eq!(Gf2Poly::zero().derivative(), Gf2Poly::zero());
    }

    #[test]
    fn truncate() {
        let mut a: Gf2Poly = "dbf".parse().unwrap();
        a.truncate_mut(7);
        assert_eq!(a.to_string(), "3f");
        a.truncate_mut(0);
        assert_eq!(a.to_string(), "0");
    }

    #[test]
    fn reverse() {
        let a: Gf2Poly = "dbf".parse().unwrap();
        assert_eq!(a.reverse().to_string(), "fdb");
    }

    #[test]
    fn reverse_gap() {
        let a = (Gf2Poly::one() << 198u64) + (Gf2Poly::one() << 135u64);
        let trail = a.trailing_zeros().unwrap_or(0);
        assert_eq!(a.reverse().reverse() << trail, a);
    }

    #[test]
    fn degree() {
        let a: Gf2Poly = "3".parse().unwrap();
        assert_eq!(a.deg(), 1);
    }
}
