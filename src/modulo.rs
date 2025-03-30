use core::borrow::Borrow;

use crate::Gf2Poly;

impl Gf2Poly {
    /// Calculates `lhs` * `rhs` modulo `self`.
    pub fn mod_mul(&self, lhs: &Self, rhs: &Self) -> Self {
        let mut res = lhs * rhs;
        if res.deg() >= self.deg() {
            res %= self;
        }
        res
    }

    /// Calculates the inverse of `elem` modulo `self`.
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let modulus: Gf2Poly = "211".parse().unwrap();
    /// let a: Gf2Poly = "b5".parse().unwrap();
    /// let inv = modulus.mod_inv(&a).unwrap();
    /// assert_eq!(inv.to_string(), "181");
    /// ```
    pub fn mod_inv(&self, elem: &Self) -> Option<Self> {
        if self.is_one() {
            return Some(Gf2Poly::zero());
        }
        let (gcd, [_, inv]) = self.clone().xgcd(elem.clone());
        if !gcd.is_one() { None } else { Some(inv) }
    }

    /// Calculates `lhs` * `rhs`^-1 modulo `self`, or None
    /// if `rhs` does not have an inverse modulo `self`.
    pub fn mod_div(&self, lhs: &Self, rhs: &Self) -> Option<Self> {
        let inv = self.mod_inv(rhs)?;
        Some(self.mod_mul(lhs, &inv))
    }

    /// Calculates `elem` * `elem` modulo `self`.
    pub fn mod_square(&self, elem: &Self) -> Self {
        let mut square = elem.square();
        if square.deg() >= self.deg() {
            square %= self;
        }
        square
    }

    /// Calculates `base`^`n` modulo `self`.
    pub fn mod_power(&self, base: &Self, mut n: u64) -> Self {
        if self.is_one() {
            return Gf2Poly::zero();
        }

        if n == 0 {
            return Gf2Poly::one();
        }

        if base.is_zero() {
            return Gf2Poly::zero();
        }

        let mut base = base;
        let mut base_val;
        let mut result = Gf2Poly::one();
        while n > 0 {
            if n % 2 == 1 {
                result = self.mod_mul(&result, base);
            }
            base_val = self.mod_square(base);
            base = &base_val;
            n >>= 1;
        }
        result
    }
}

// If we have a fixed modulus, we can generally speed up modular reductions using Barrett reduction.
// See for example Intel 323102:
// https://github.com/tpn/pdfs/blob/a8256e545a4f4fd31342d3750a3a85b7c58ee44f/Fast%20CRC%20Computation%20for%20Generic%20Polynomials%20Using%20PCLMULQDQ%20Instruction%20-%20Intel%20(December%2C%202009).pdf
// Here's why Barrett reduction works with polynomials:
// Let g be a polynomial of degree n and f a polynomial of degree < n.
// Also let q, r be such that f * x^n = g*q + r with deg(r) < n.
// (note that adding something below degree n on the left side does not change q).
// Now let q', r' be such that x^2n = g*q' + r' with deg(r') < n.
// We have f * (g*q' + r') = x^n * (g*q + r*x^n).
// Both sides have degree 3n, so if we take the most n+1 significant bits, we get
// (f * g * q') >> 2*n = (x^n * g * q) >> 2*n
// since both deg(f * r') < 2n and deg(x^n * r) < 2n.
// We also have deg(g) = n, and since the low 2n bits from the result don't matter,
// we can change the left expression to the equivalent ((f * q' >> n) * g) >> n.
// The right side is also just (q * g) >> n.
// Note that h |-> (h * g) >> n is a linear map and forms a triangular matrix,
// which means it is invertible. Hence, from
// ((f * q' >> n) * g) >> n = (q * g) >> n
// we get q = (f * q') >> n.

/// This type precalculates the factor for Barrett reduction, which can be used to
/// speed up calculations using the same modulus.
/// This takes a Borrow<Gf2Poly>, which can for example be a Gf2Poly, a &Gf2Poly or
/// an Arc<Gf2Poly>
pub struct Gf2PolyMod<T: Borrow<Gf2Poly>> {
    modulus: T,
    barrett_reducer: Gf2Poly,
}

impl<T: Borrow<Gf2Poly>> Gf2PolyMod<T> {
    /// Constructs a new [Gf2PolyMod] and precalculates data for faster remaindering.
    /// Panics if `modulus` is zero.
    pub fn new(modulus: T) -> Self {
        if modulus.borrow().is_zero() {
            panic!("Zero modulus is not allowed.");
        }

        let barrett_reducer =
            Gf2Poly::x_to_the_power_of(2 * modulus.borrow().deg()) / modulus.borrow();
        Gf2PolyMod {
            modulus,
            barrett_reducer,
        }
    }

    /// Returns a reference to the original modulus this [Gf2PolyMod] was constructed with.
    pub fn modulus(&self) -> &Gf2Poly {
        &self.modulus.borrow()
    }

    /// Returns the modulus value this [Gf2PolyMod] was constructed with.
    pub fn modulus_value(self) -> T {
        self.modulus
    }

    /// Returns the degree of the modulus.
    pub fn deg(&self) -> u64 {
        self.modulus().deg()
    }

    fn barrett_step(&self, upper_half: Gf2Poly) -> Gf2Poly {
        // barrett_reducer is the q' above, and upper_half is f
        (upper_half * &self.barrett_reducer) >> self.deg()
    }

    fn barrett_remainder(&self, poly: &Gf2Poly) -> Gf2Poly {
        let quotient = self.barrett_step(poly >> self.deg());
        poly - quotient * self.modulus()
    }

    /// Calculates the remainder of `elem` when divided by `self`.
    ///
    /// ## Example
    /// ```rust
    /// # use gf2poly::{Gf2Poly, Gf2PolyMod};
    /// let m: Gf2Poly = "b1".parse().unwrap();
    /// let modulus = Gf2PolyMod::new(m);
    /// let a: Gf2Poly = "2f7b7".parse().unwrap();
    /// let remainder = modulus.remainder(&a);
    /// assert_eq!(remainder.to_string(), "3");
    /// ```
    pub fn remainder(&self, elem: &Gf2Poly) -> Gf2Poly {
        if elem.deg() < self.deg() {
            return elem.clone();
        }

        if self.modulus().is_one() {
            return Gf2Poly::zero();
        }

        let step = self.deg();
        if elem.deg() < 2 * step {
            return self.barrett_remainder(elem);
        }

        let last_segment = elem.deg() / step;
        let range = |segment: u64| segment * step..(segment + 1) * step;
        let mut remainder = elem.subrange(range(last_segment));

        for segment in (0..last_segment).rev() {
            remainder <<= step;
            remainder += elem.subrange(range(segment));
            remainder = self.barrett_remainder(&remainder);
        }

        remainder
    }

    /// Performs multiplication of `lhs` with `rhs` modulo `self`.
    pub fn mul(&self, lhs: &Gf2Poly, rhs: &Gf2Poly) -> Gf2Poly {
        let product = lhs * rhs;
        if product.deg() < self.deg() {
            return product;
        }
        self.remainder(&product)
    }

    /// Performs squaring of `elem` modulo `self`.
    pub fn square(&self, elem: &Gf2Poly) -> Gf2Poly {
        let square = elem.square();
        if square.deg() < self.deg() {
            return square;
        }
        self.remainder(&square)
    }

    /// Convenience method for [Gf2Poly::mod_inv];
    /// Is not more efficient than using self.modulus() directly.
    pub fn inverse(&self, elem: &Gf2Poly) -> Option<Gf2Poly> {
        self.modulus().mod_inv(elem)
    }

    /// Convenience method for [Gf2Poly::mod_div].
    /// Is not more efficient than using self.modulus() directly.
    pub fn div(&self, lhs: &Gf2Poly, rhs: &Gf2Poly) -> Option<Gf2Poly> {
        self.modulus().mod_div(lhs, rhs)
    }
}

#[cfg(test)]
mod tests {
    use crate::prop_assert_poly_eq;

    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn modmul_is_mul_remainder(modulus: Gf2Poly, a: Gf2Poly, b: Gf2Poly) {
            let res = modulus.mod_mul(&a, &b);
            let rem = &a * &b % &modulus;
            prop_assert_poly_eq!(res, rem);
        }

        #[test]
        fn modulus_mul(modulo: Gf2Poly, a: Gf2Poly, b: Gf2Poly) {
            prop_assume!(!modulo.is_zero());
            let res1 = modulo.mod_mul(&a, &b);
            let modulus = Gf2PolyMod::new(&modulo);
            let res2 = modulus.mul(&a, &b);
            prop_assert_poly_eq!(res1, res2);
        }

        #[test]
        fn modular_inv(modulus: Gf2Poly, elem: Gf2Poly) {
            let inv = modulus.mod_inv(&elem);
            let elem = &elem % &modulus;
            if let Some(inv) = inv {
                prop_assert_eq!(modulus.mod_mul(&elem, &inv), Gf2Poly::one() % modulus);
            }
        }

        #[test]
        fn modulus_remainder(modulus: Gf2Poly, elem: Gf2Poly) {
            prop_assume!(!modulus.is_zero());
            let modulus = Gf2PolyMod::new(modulus);
            let res1 = modulus.remainder(&elem);
            let res2 = elem % &modulus.modulus;
            prop_assert_poly_eq!(res1, res2);
        }

        #[test]
        fn modular_div(modulo: Gf2Poly, lhs: Gf2Poly, rhs: Gf2Poly) {
            prop_assume!(!rhs.is_zero());
            let div = modulo.mod_div(&lhs, &rhs);
            if let Some(div) = div {
                prop_assert_eq!(modulo.mod_mul(&div, &rhs), lhs % modulo);
            }
        }

        #[test]
        fn exponent_homo(modulo: Gf2Poly, a: Gf2Poly, n in 0..128u64, m in 0..128u64) {
            prop_assert_poly_eq!(
                modulo.mod_mul(&modulo.mod_power(&a, n), &modulo.mod_power(&a, m)),
                modulo.mod_power(&a, n + m));
        }

        #[test]
        fn power_homo(modulo: Gf2Poly, a: Gf2Poly, b: Gf2Poly, n in 0..128u64) {
            prop_assert_poly_eq!(
                modulo.mod_power(&modulo.mod_mul(&a, &b), n),
                modulo.mod_mul(&modulo.mod_power(&a, n), &modulo.mod_power(&b, n))
            );
        }
    }
}
