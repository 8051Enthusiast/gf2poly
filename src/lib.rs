#![debugger_visualizer(gdb_script_file = "gf2poly_dbg.py")]
#![cfg_attr(not(test), no_std)]
extern crate alloc;

pub use rand;

mod add;
mod div;
mod factor;
mod gcd;
mod matrix;
mod modulo;
mod mul;
mod shift;
#[cfg(test)]
mod test_util;

pub use modulo::Gf2PolyMod;

use alloc::vec::Vec;
use core::{fmt::Display, str::FromStr};
use rand::Rng;

pub type Limb = core::ffi::c_ulong;

// a plain vec is easier for debugging
#[cfg(debug_assertions)]
type LimbStorage = Vec<Limb>;
#[cfg(not(debug_assertions))]
type LimbStorage = smallvec::SmallVec<[Limb; 2]>;

#[cfg(debug_assertions)]
use alloc::vec as limbs;
#[cfg(not(debug_assertions))]
use smallvec::smallvec as limbs;

const BITS: u64 = Limb::BITS as u64;

const HALF_MASK: Limb = alternating_mask(BITS.ilog2() - 1);

fn normalized_len(limbs: &[Limb]) -> usize {
    limbs.len() - limbs.iter().rev().take_while(|&&x| x == 0).count()
}

fn deg(limbs: &[Limb]) -> u64 {
    let Some(last) = limbs.last() else {
        return 0;
    };
    limbs.len() as u64 * BITS - 1 - last.leading_zeros() as u64
}

fn limbs_for_deg(n: u64) -> usize {
    usize::try_from(n / BITS).unwrap() + 1
}

fn limb_degree(a: Limb) -> u8 {
    (BITS as u8).saturating_sub(1 + a.leading_zeros() as u8)
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

    fn normalize(&mut self) {
        self.limbs.truncate(normalized_len(&self.limbs));
        self.deg = deg(&self.limbs);
    }

    /// Degree of the polynomial. Returns degree 0 for the zero polynomial.
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let a: Gf2Poly = "7".parse().unwrap();
    /// assert_eq!(a.deg(), 2);
    /// ```
    pub fn deg(&self) -> u64 {
        self.deg
    }

    /// Whether the polynomial is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.limbs().is_empty()
    }

    /// Whether the polynomial is either the zero or one polynomial.
    pub fn is_constant(&self) -> bool {
        self.deg() == 0
    }

    /// Whether the polynomial is the 1 polynomial.
    pub fn is_one(&self) -> bool {
        self.is_constant() && !self.is_zero()
    }

    /// Gets the nth coefficient of the polynomial,
    /// with the 0th coefficient being the constant part.
    pub fn get(&self, n: u64) -> bool {
        let Ok(idx) = usize::try_from(n / BITS) else {
            return false;
        };
        match self.limbs().get(idx) {
            Some(limb) => {
                let res = n % BITS;
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
        let idx = n / BITS;
        let res = n % BITS;
        self.limbs[usize::try_from(idx).unwrap()] |= 1 << res;
    }

    /// Sets the nth coefficient of the polynomial to 0,
    /// with the 0th coefficient being the constant part.
    pub fn clear(&mut self, n: u64) {
        if n > self.deg {
            return;
        }
        let idx = n / BITS;
        let res = n % BITS;
        self.limbs[usize::try_from(idx).unwrap()] &= !(1 << res);
        if n == self.deg {
            self.normalize();
        }
    }

    /// Sets the coefficient with degree `n`` to the value in `coefficient`.
    pub fn insert(&mut self, n: u64, coefficient: bool) {
        if coefficient {
            self.set(n);
        } else {
            self.clear(n);
        }
    }

    /// Converts a vector of limbs into a polynomial, with least significant limb first.
    fn from_limb_storage(mut limbs: LimbStorage) -> Self {
        limbs.truncate(normalized_len(&limbs));
        let deg = deg(&limbs);
        Gf2Poly { limbs, deg }
    }

    /// Converts a slice of limbs into a polynomial, with least significant limb first.
    pub fn from_limbs(limbs: &[Limb]) -> Self {
        let len = normalized_len(limbs);
        let limb_vec = LimbStorage::from(&limbs[..len]);
        let deg = deg(&limbs[..len]);
        Gf2Poly {
            limbs: limb_vec,
            deg,
        }
    }

    /// Converts a single limb into a polynomial.
    pub fn from_limb(limb: Limb) -> Self {
        if limb == 0 {
            return Gf2Poly::zero();
        }
        let mut limbs = LimbStorage::with_capacity(1);
        limbs.push(limb);
        Gf2Poly {
            deg: limb_degree(limb) as u64,
            limbs,
        }
    }

    pub fn limbs(&self) -> &[Limb] {
        &self.limbs
    }

    /// Converts a slice of bytes into a polynomial, with least significant byte first.
    pub fn from_bytes(bytes: &[u8]) -> Gf2Poly {
        let limb_size = core::mem::size_of::<Limb>();
        let mut limbs = LimbStorage::with_capacity(bytes.len().div_ceil(limb_size));
        for chunk in bytes.chunks(core::mem::size_of::<Limb>()) {
            let mut limb = 0;
            for (i, &byte) in chunk.iter().enumerate() {
                limb |= (byte as Limb) << (i << 3);
            }
            limbs.push(limb);
        }
        Gf2Poly::from_limb_storage(limbs)
    }

    /// Converts a polynomial into a vector of bytes, with least significant byte first.
    /// The vector is always as short as possible, which means that this method returns
    /// an empty vector for 0.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity((self.deg() / 8 + 1) as usize);
        let Some((last, limbs)) = self.limbs().split_last() else {
            return Vec::new();
        };
        for limb in limbs {
            bytes.extend_from_slice(&limb.to_le_bytes());
        }
        let last_bytes = (self.deg() % BITS) / 8 + 1;
        bytes.extend_from_slice(&last.to_le_bytes()[..last_bytes as usize]);
        bytes
    }

    /// 0 (Zero).
    /// The polynomial with all coefficients 0.
    /// The unique polynomial that is divisible by every other polynomial.
    /// The unique fixpoint of the derivative operator on GF(2)\[x\].
    pub fn zero() -> Self {
        Default::default()
    }

    /// 1
    pub fn one() -> Self {
        Gf2Poly {
            deg: 0,
            limbs: core::iter::once(1).collect(),
        }
    }

    /// Returns the x of GF(2)\[x\].
    pub fn x() -> Self {
        Gf2Poly {
            deg: 1,
            limbs: core::iter::once(2).collect(),
        }
    }

    /// Returns x^n
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let a = Gf2Poly::x_to_the_power_of(3);
    /// assert_eq!(a.to_string(), "8");
    /// ```
    pub fn x_to_the_power_of(n: u64) -> Self {
        let mut x = Gf2Poly::zero();
        x.set(n);
        x
    }

    /// Calculates the derivative of self.
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let a: Gf2Poly = "dbf".parse().unwrap();
    /// let b = a.derivative();
    /// assert_eq!(b.to_string(), "455");
    /// assert_eq!(Gf2Poly::zero().derivative(), Gf2Poly::zero());
    /// ```
    pub fn derivative(&self) -> Self {
        const MASK: Limb = Limb::MAX / 3;
        let limbs = self.limbs().iter().map(|x| (x >> 1) & MASK).collect();
        Self::from_limb_storage(limbs)
    }

    /// Calculates the number of coefficients with value 1.
    pub fn hamming_weight(&self) -> u64 {
        self.limbs().iter().map(|x| x.count_ones() as u64).sum()
    }

    /// Calculates the polynomial with all coefficients reversed in place.
    /// If f is of degree n, this is x^n * f(1/x).
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let a: Gf2Poly = "dbf".parse().unwrap();
    /// assert_eq!(a.reverse().to_string(), "fdb");
    /// ```
    pub fn reverse(&self) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        let mut limbs;
        let bit_offset = (self.deg + 1) % BITS;
        if bit_offset == 0 {
            limbs = self
                .limbs()
                .iter()
                .rev()
                .copied()
                .map(Limb::reverse_bits)
                .collect();
        } else {
            limbs = LimbStorage::with_capacity(self.limbs().len());
            let mut prev = 0;
            let mut first = true;
            for x in self.limbs().iter().rev() {
                let rev = x.reverse_bits();
                if !first {
                    limbs.push((rev << bit_offset) | prev);
                } else {
                    first = false;
                }
                prev = rev >> (BITS - bit_offset);
            }
            if prev != 0 {
                limbs.push(prev);
            }
        }
        Self::from_limb_storage(limbs)
    }

    /// Calculates n for f such that x^n || f, i.e. the number of trailing zeros.
    /// Returns None if f is zero.
    pub fn trailing_zeros(&self) -> Option<u64> {
        if self.is_zero() {
            return None;
        }
        let mut zeros = 0;
        for limb in self.limbs().iter() {
            if *limb == 0 {
                zeros += BITS;
            } else {
                zeros += limb.trailing_zeros() as u64;
                break;
            }
        }
        Some(zeros)
    }

    /// Calculates f mod x^n, i.e. truncates f to the n least significant bits.
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let mut a: Gf2Poly = "dbf".parse().unwrap();
    /// a.truncate_mut(7);
    /// assert_eq!(a.to_string(), "3f");
    /// a.truncate_mut(0);
    /// assert_eq!(a.to_string(), "0");
    /// ```
    pub fn truncate_mut(&mut self, n: u64) {
        if n == 0 {
            *self = Self::zero();
            return;
        }
        if self.deg < n {
            return;
        }
        let n_bits = ((n - 1) % BITS) as usize;
        let Ok(n_limbs) = usize::try_from((n - 1) / BITS) else {
            return;
        };
        self.limbs.truncate(n_limbs + 1);
        self.limbs[n_limbs] &= (1 << n_bits) | ((1 << n_bits) - 1);
        self.normalize();
    }

    /// Calculates f mod x^n, i.e. truncates f to the n least significant bits.
    pub fn truncated(&self, n: u64) -> Self {
        if n == 0 {
            return Self::zero();
        }
        self.subrange(..n)
    }

    /// Calculates self^n, the nth power of self.
    pub fn power(&self, mut n: u64) -> Self {
        if n == 0 {
            return Gf2Poly::one();
        }

        if self.is_zero() {
            return Gf2Poly::zero();
        }

        let mut result = Gf2Poly::one();
        let mut base = self;
        let mut base_val;
        while n > 0 {
            if n % 2 == 1 {
                result *= base;
            }
            base_val = base.square();
            base = &base_val;
            n >>= 1;
        }
        result
    }

    /// Given `self`, calculates `self` * `self`.
    pub fn square(&self) -> Self {
        if self.is_zero() {
            return Self::zero();
        }

        let deg = self.deg * 2;
        let mut limbs = LimbStorage::with_capacity(limbs_for_deg(deg));
        let (last_limb, start_limbs) = self.limbs().split_last().unwrap();
        let lo_half = |x: Limb| x & HALF_MASK;
        let hi_half = |x: Limb| (x >> (BITS / 2)) & HALF_MASK;
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

    /// Given f, calculates its squareroot.
    /// Returns None if there is none.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero());
        }

        if self.deg % 2 != 0 {
            return None;
        }

        let deg = self.deg / 2;
        let mut limbs = LimbStorage::with_capacity(limbs_for_deg(deg));

        let is_square_limb = |x: Limb| (x & !alternating_mask(0)) == 0;

        for elements in self.limbs().chunks_exact(2) {
            let [lo, hi] = elements.try_into().unwrap();
            if !is_square_limb(lo) || !is_square_limb(hi) {
                return None;
            }
            limbs.push(unspacen_bits(lo) | (unspacen_bits(hi) << (BITS / 2)));
        }

        if self.limbs().len() % 2 != 0 {
            let last = *self.limbs().last().unwrap();
            if !is_square_limb(last) {
                return None;
            }
            limbs.push(unspacen_bits(last))
        }

        Some(Gf2Poly { deg, limbs })
    }

    /// Evaluates the polynomial at the point x.
    pub fn eval(&self, x: bool) -> bool {
        if !x {
            1 & *self.limbs().first().unwrap_or(&0) != 0
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
        let hi = self >> n;
        let lo = self.truncated(n);
        (hi, lo)
    }

    /// Returns a uniformly random polynomial of degree `deg` using the given RNG.
    pub fn random<R: Rng>(deg: u64, rng: &mut R) -> Self {
        if deg == 0 {
            return Self::one();
        }
        let len = limbs_for_deg(deg);
        let mut limbs = limbs![0; len];
        rng.fill(limbs.as_mut_slice());
        let last = &mut limbs[len - 1];
        let deg_part = deg % BITS;
        *last &= (1 << deg_part) - 1;
        *last |= 1 << deg_part;
        Gf2Poly { deg, limbs }
    }

    #[cfg(test)]
    fn is_normalized(&self) -> bool {
        self.limbs().is_empty() && self.deg == 0 || self.limbs().len() == limbs_for_deg(self.deg)
    }
}

// creates an integer where, for each bit with index i,
// the bit is set if the nth bit of the bit index is set.
// E.g. 0 => 0x5555..., 1 => 0x33333..., 2 => 0x0f0f0f...
const fn alternating_mask(n: u32) -> Limb {
    !0 / ((1 << (1 << n)) + 1)
}

fn spacen_bits(mut x: Limb) -> Limb {
    for i in (0..(Limb::BITS.ilog2() - 1)).rev() {
        x = (x | (x << (1 << i))) & alternating_mask(i);
    }
    x
}

fn unspacen_bits(mut x: Limb) -> Limb {
    for i in 0..(Limb::BITS.ilog2() - 1) {
        x = (x | (x >> (1 << i))) & alternating_mask(i + 1);
    }
    x
}

impl Display for Gf2Poly {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let Some((last, limbs)) = self.limbs().split_last() else {
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
        let mut limbs = LimbStorage::with_capacity(s.len() / BITS as usize + 1);
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
    use rand::SeedableRng;

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
        fn sqrt_is_inverse_of_square(a: Gf2Poly) {
            prop_assert_poly_eq!(a.square().sqrt().unwrap(), a);
        }

        #[test]
        fn non_square_has_no_root(a: Gf2Poly) {
            prop_assume!(!a.is_zero());
            let non_square = a.square() * Gf2Poly::x();
            prop_assert!(non_square.sqrt().is_none());
        }

        #[test]
        fn exponent_homo(a: Gf2Poly, n in 0..64u64, m in 0..64u64) {
            prop_assert_poly_eq!(a.clone().power(n) * a.clone().power(m), a.power(n + m));
        }

        #[test]
        fn power_homo(a: Gf2Poly, b: Gf2Poly, n in 0..64u64) {
            prop_assert_poly_eq!((&a * &b).power(n), a.power(n) * b.power(n));
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

        #[test]
        fn random_polynomial(n in 0..1024u64, seed: u64) {
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            let random = Gf2Poly::random(n, &mut rng);
            prop_assert!(random.is_normalized());
            prop_assert_eq!(random.deg(), n);
        }

        #[test]
        fn bytes_roundtrip(a: Gf2Poly) {
            let bytes = a.to_bytes();
            let b = Gf2Poly::from_bytes(&bytes);
            prop_assert!(bytes.last() != Some(&0));
            prop_assert_eq!(b.to_bytes(), bytes);
            prop_assert_poly_eq!(a, b);
        }
    }

    #[test]
    fn reverse_gap() {
        let a = (Gf2Poly::one() << 198u64) + (Gf2Poly::one() << 135u64);
        let trail = a.trailing_zeros().unwrap_or(0);
        assert_eq!(a.reverse().reverse() << trail, a);
    }
}
