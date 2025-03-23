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
        if !gcd.is_one() {
            return None;
        } else {
            return Some(inv);
        }
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

#[cfg(test)]
mod tests {
    use crate::prop_assert_poly_eq;

    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn modmul_associativity(modulo: Gf2Poly, a: Gf2Poly, b: Gf2Poly, c: Gf2Poly) {
            prop_assert_poly_eq!(
                modulo.mod_mul(&modulo.mod_mul(&a, &b), &c),
                modulo.mod_mul(&a, &modulo.mod_mul(&b, &c))
            );

        }

        #[test]
        fn modmul_commutativity(modulo: Gf2Poly, a: Gf2Poly, b: Gf2Poly) {
            prop_assert_poly_eq!(
                modulo.mod_mul(&a, &b),
                modulo.mod_mul(&b, &a)
            );
        }

        #[test]
        fn modmul_identity(modulo: Gf2Poly, a: Gf2Poly) {
            prop_assert_poly_eq!(
                modulo.mod_mul(&a, &Gf2Poly::one()),
                a % modulo
            );
        }

        #[test]
        fn modmul_distributivity(modulo: Gf2Poly, a: Gf2Poly, b: Gf2Poly, c: Gf2Poly) {
            prop_assert_poly_eq!(
                modulo.mod_mul(&a, &(&b + &c)),
                modulo.mod_mul(&a, &b) + modulo.mod_mul(&a, &c)
            );
        }

        #[test]
        fn modular_inv(modulo: Gf2Poly, elem: Gf2Poly) {
            prop_assume!(!elem.is_zero());
            let inv = modulo.mod_inv(&elem);
            let elem = &elem % &modulo;
            if let Some(inv) = inv {
                prop_assert_eq!(modulo.mod_mul(&elem, &inv), Gf2Poly::one() % modulo);
            }
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
