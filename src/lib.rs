#![debugger_visualizer(gdb_script_file = "gf2poly_dbg.py")]
#![cfg_attr(not(test), no_std)]

extern crate alloc;

pub use rand;

mod add;
mod div;
mod factor;
mod gcd;
mod matrix;
mod mul;
mod shift;
#[cfg(test)]
mod test_util;
use core::{fmt::Display, str::FromStr};
use rand::Rng;

pub type Limb = core::ffi::c_ulong;
// a plain vec is easier for debugging
#[cfg(debug_assertions)]
type LimbStorage = alloc::vec::Vec<Limb>;
#[cfg(not(debug_assertions))]
type LimbStorage = smallvec::SmallVec<[Limb; 2]>;

const BITS: usize = Limb::BITS as usize;

const HALF_MASK: Limb = alternating_mask(BITS.ilog2() - 1);

fn normalized_len(limbs: &[Limb]) -> usize {
    limbs.len() - limbs.iter().rev().take_while(|&&x| x == 0).count()
}

fn deg(limbs: &[Limb]) -> u64 {
    let Some(last) = limbs.last() else {
        return 0;
    };
    limbs.len() as u64 * BITS as u64 - 1 - last.leading_zeros() as u64
}

fn limbs_for_deg(n: u64) -> usize {
    usize::try_from(n as u64 / BITS as u64).unwrap() + 1
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

    /// Gets the nth coefficient of the polynomial,
    /// with the 0th coefficient being the constant part.
    pub fn get(&self, n: u64) -> bool {
        let Ok(idx) = usize::try_from(n / BITS as u64) else {
            return false;
        };
        match self.limbs().get(idx) {
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

    /// Sets the nth coefficient of the polynomial to 0,
    /// with the 0th coefficient being the constant part.
    pub fn clear(&mut self, n: u64) {
        if n > self.deg {
            return;
        }
        let idx = n / BITS as u64;
        let res = n % BITS as u64;
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
        let len = normalized_len(&limbs);
        let limb_vec = LimbStorage::from(&limbs[..len]);
        let deg = deg(&limbs);
        Gf2Poly {
            limbs: limb_vec,
            deg,
        }
    }

    pub fn limbs(&self) -> &[Limb] {
        &self.limbs
    }

    /// Converts a slice of bytes into a polynomial, with least significant byte first.
    pub fn from_bytes(bytes: &[u8]) -> Gf2Poly {
        let limb_size = core::mem::size_of::<Limb>();
        let mut limbs = LimbStorage::with_capacity(bytes.len() / limb_size + 1);
        for chunk in bytes.chunks(core::mem::size_of::<Limb>()) {
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

    /// Calculates the derivative of self.
    pub fn derivative(&self) -> Self {
        const MASK: Limb = Limb::MAX / 3;
        let limbs = self.limbs().iter().map(|x| x >> 1 & MASK).collect();
        Self::from_limb_storage(limbs)
    }

    /// Calculates the number of coefficients with value 1.
    pub fn hamming_weight(&self) -> u64 {
        self.limbs().iter().map(|x| x.count_ones() as u64).sum()
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
        self.normalize();
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

    /// Calculates self^n, the nth power of self.
    pub fn power(self, mut n: u64) -> Self {
        if n == 0 {
            return Gf2Poly::one();
        }

        if self.is_zero() {
            return Gf2Poly::zero();
        }

        let mut result = Gf2Poly::one();
        let mut base = self;
        while n > 0 {
            if n % 2 == 1 {
                result *= &base;
            }
            base = base.square();
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
            limbs.push(unspacen_bits(lo) | unspacen_bits(hi) << (BITS / 2));
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
        let hi = self.clone() >> n;
        let lo = self.truncated(n);
        (hi, lo)
    }

    /// Returns a uniformly random polynomial of degree `deg` using the given RNG.
    pub fn random<R: Rng>(deg: u64, rng: &mut R) -> Self {
        if deg == 0 {
            return Self::one();
        }
        let len = limbs_for_deg(deg);
        let mut limbs = LimbStorage::with_capacity(len);
        limbs.resize(len, 0);
        rng.fill(limbs.as_mut_slice());
        let last = &mut limbs[len - 1];
        let deg_part = deg % BITS as u64;
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
        fn exponent_homo(a: Gf2Poly, n in 0..128u64, m in 0..128u64) {
            prop_assert_poly_eq!(a.clone().power(n) * a.clone().power(m), a.power(n + m));
        }

        #[test]
        fn power_homo(a: Gf2Poly, b: Gf2Poly, n in 0..128u64) {
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
        fn random_polynomial(n in 0..1024u64) {
            let mut rng = rand::rng();
            let random = Gf2Poly::random(n, &mut rng);
            prop_assert!(random.is_normalized());
            prop_assert_eq!(random.deg(), n);
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
