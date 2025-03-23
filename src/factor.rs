use rand::SeedableRng;

use crate::Gf2Poly;
use alloc::vec::Vec;
use alloc::vec;

// Wikipedia has an excellent page on factorization of polynomials over finite fields:
// https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields

impl Gf2Poly {
    fn square_free_factorization_impl(&self, multiplicity_multiplier: u64) -> Vec<(Gf2Poly, u64)> {
        let mut factors = Vec::new();

        // the gcd of a polynomial with its derivative decreases the multiplicities
        // of all roots (except for squares, see later where we recurse) by one
        let mut duplicate_part = self.clone().gcd(self.derivative());
        // if we divide the our original by the polynomial with multiplicities reduced,
        // the result will have multiplicities of most 1, ie it is squarefree
        let mut squarefree = self / &duplicate_part;

        let mut multiplicity = 0;
        while squarefree.deg() > 0 {
            // here, for each factor with multiplicity e_i, the result will have a multiplicity
            // * at most e_i - 1 (because of the reduced multiplicities in duplicate_part)
            // * at most 1 (because squarefree only has multiplicity at most 1)
            // the result is that we have multiplicity 1 for each (non-square) e_i > 1 and 0 for others
            let other_uniq_factors = duplicate_part.clone().gcd(squarefree.clone());
            // and if we remove those from squarefree, we get the complement, ie all those with e_i = 1
            let factor = &squarefree / &other_uniq_factors;

            // prepare for next iteration
            duplicate_part /= &other_uniq_factors;
            squarefree = other_uniq_factors;
            multiplicity += multiplicity_multiplier;
            if factor.deg() > 0 {
                factors.push((factor, multiplicity));
            }
        }

        // if we have a square in characteristics 2, its derivative will be zero:
        // due to the freshman's dream, the square will only have mononomials of
        // even degree, and the derivative of them will always be zero because
        // they will be multiplied by 2.
        // this can be worked around by just dividing all exponents by 2 and running
        // the same algorithm again, with multiplicities doubled as everything
        // is really a square
        if duplicate_part.deg() > 0 {
            let mut t = duplicate_part
                .sqrt()
                .expect("duplicate_part should be a square at this point")
                .square_free_factorization_impl(multiplicity_multiplier * 2);
            factors.append(&mut t);
        }

        factors
    }

    /// Computes the square-free factorization of the polynomial.
    /// Given a polynomial f = prod_{i=0}^n {f_i}^i, gives back
    /// pairs (f_i, i) with f != 1 such that f is squarefree.
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let a: Gf2Poly = "78f314da3a4".parse().unwrap();
    /// let result = a.square_free_factorization();
    /// let [(a1, 1), (a2, 2), (a3, 3)] = &result[..] else {
    ///     panic!()
    /// };
    /// assert_eq!(a1.to_string(), "ed");
    /// assert_eq!(a2.to_string(), "da");
    /// assert_eq!(a3.to_string(), "b5");
    /// ```
    pub fn square_free_factorization(&self) -> Vec<(Gf2Poly, u64)> {
        if self.is_zero() {
            panic!("Cannot factorize zero");
        }

        let mut factors = self.square_free_factorization_impl(1);
        factors.sort_unstable_by_key(|(_, mult)| *mult);
        factors
    }

    /// Computes the distinct degree factorization of the polynomial.
    /// The input polynomial must be squarefree.
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let a: Gf2Poly = "3813c0be6".parse().unwrap();
    /// let result = a.distinct_degree_factorization();
    /// let [(a1, 1), (a2, 2), (a5, 5), (a12, 12)] = &result[..] else {
    ///     panic!();
    /// };
    /// assert_eq!(a1.to_string(), "6");
    /// assert_eq!(a2.to_string(), "7");
    /// assert_eq!(a5.to_string(), "37");
    /// assert_eq!(a12.to_string(), "1764ac5");
    /// ```
    pub fn distinct_degree_factorization(&self) -> Vec<(Gf2Poly, u64)> {
        if self.is_zero() {
            panic!("Cannot factorize zero");
        }

        let mut current;
        let mut current_ref = self;
        // hypersquare is here x^(2^k) mod current_ref for the kth iteration of the loop (starting from 1)
        let mut hypersquare = Gf2Poly::x();
        let mut degree = 0;
        let mut factors = Vec::new();

        while current_ref.deg() > 0 {
            degree += 1;
            if current_ref.deg() == degree {
                factors.push((current_ref.clone(), degree));
                break;
            }
            hypersquare = current_ref.mod_square(&hypersquare);

            // the `product` variable `x^(2^k) + x` is the product of all polynomials whose
            // degree divides `k`
            let product = &hypersquare + Gf2Poly::x();
            // note that `product` is actually modulo `current_ref`, but since we calculate
            // the gcd with `current_ref`, it doesn't make a difference for the gcd calculation.
            // as a result, we get all polynomials dividing `current_ref` whose degree divides `k`/`degree`
            let common = product.gcd(current_ref.clone());
            if common.deg() == 0 {
                continue;
            }

            // make sure to remove the current polynomials because otherwise they will be found again
            // in iteration 2*k
            current = current_ref / &common;
            current_ref = &current;
            factors.push((common, degree));
        }

        factors.sort_unstable_by_key(|(_, mult)| *mult);
        factors
    }

    // the trace Tr_L/K(a) of a field extension L/K where L = K/f for some irreducible polynomial
    // f can be calculate as `sum(a^(2^i) for i in range(deg(f)))` and the result will be
    // an element of the subfield K.
    // here we do that calculation modulo `self` where `self` is not necessarily irreducible
    // and with a given degree, not the one of the polynomial, hence _pseudo_ trace
    fn pseudo_trace(&self, deg: u64, a: &Gf2Poly) -> Gf2Poly {
        let mut square_ref = a;
        let mut square;
        if a.deg() >= self.deg() {
            square = a % self;
            square_ref = &square;
        }
        let mut res = Gf2Poly::zero();
        for _ in 0..deg - 1 {
            res += square_ref;
            square = self.mod_square(square_ref);
            square_ref = &square;
        }
        res += square_ref;
        res
    }

    // typically, for odd characteristics we typically raise a power that can be either
    // 1 or -1 (or 0 in rare cases) modulo an irreducible power, however this doesn't work out in char 2
    // as 1 == -1. instead we use the trace, which should be 0 or 1 with roughly same probability
    fn same_degree_factorization_split(&self, deg: u64, a: &Gf2Poly) -> Option<(Gf2Poly, Gf2Poly)> {
        // by the chinese remainder theorem, this is the same as taking the trace in each individual
        // irreducible factor, where it should be either 0 or 1 modulo that factor
        let trace = self.pseudo_trace(deg, a);
        // if we take the gcd with self, we get all those factors where the result was 0 and not 1,
        // which should be roughly half of them
        let common = self.clone().gcd(trace);
        if common.deg() == 0 || common.deg() == self.deg() {
            return None;
        }
        let complement = self / &common;
        Some((common, complement))
    }

    /// Takes a degree `deg` and a product of distinct irreducible
    /// polynomials of degree `deg` and returns the individual factors.
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let a: Gf2Poly = "1764ac5".parse().unwrap();
    /// let result = a.same_degree_factorization(12);
    /// let [a1, a2] = &result[..] else {
    ///     panic!();
    /// };
    /// assert_eq!(a1.to_string(), "1823");
    /// assert_eq!(a2.to_string(), "1a63");
    /// ```
    pub fn same_degree_factorization(&self, deg: u64) -> Vec<Gf2Poly> {
        if self.is_zero() {
            panic!("Cannot factorize zero");
        }

        if deg == 0 {
            panic!("Cannot factorize into zero degree polynomials");
        }

        if self.deg() == deg {
            return vec![self.clone()];
        }

        if deg == 1 && self.deg() == 2 {
            return vec![Gf2Poly::x(), Gf2Poly::x() + Gf2Poly::one()];
        }

        let mut factors = Vec::new();
        let mut todo = vec![self.clone()];
        // chosen by fair dice roll
        let mut rng = rand::rngs::StdRng::seed_from_u64(0x19b88918bffa85d);
        while let Some(current) = todo.pop() {
            if current.deg() == deg {
                factors.push(current);
                continue;
            }

            loop {
                let mut rand = Self::random(current.deg(), &mut rng);
                rand.clear(current.deg());

                if let Some((common, complement)) =
                    current.same_degree_factorization_split(deg, &rand)
                {
                    todo.push(common);
                    todo.push(complement);
                    break;
                }
            }
        }
        factors.sort_unstable();
        factors
    }

    /// Computes the factorization of the polynomial.
    /// The result is a vector of irreducible polynomials and their multiplicities.
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// // 2^2 * 3^2 * 7^2 * d * 29 * 2f * bc9
    /// let a: Gf2Poly = "1d9247f2c".parse().unwrap();
    /// let result = a.factor();
    /// let [(g_d, 1), (g_29, 1), (g_2f, 1), (g_bc9, 1), (g_2, 2), (g_3, 2), (g_7, 2)] = &result[..] else  {
    ///     panic!()
    /// };
    /// assert_eq!(g_d.to_string(), "d");
    /// assert_eq!(g_29.to_string(), "29");
    /// assert_eq!(g_2f.to_string(), "2f");
    /// assert_eq!(g_bc9.to_string(), "bc9");
    /// assert_eq!(g_2.to_string(), "2");
    /// assert_eq!(g_3.to_string(), "3");
    /// assert_eq!(g_7.to_string(), "7");
    /// ```
    pub fn factor(&self) -> Vec<(Gf2Poly, u64)> {
        if self.is_zero() {
            panic!("Cannot factorize zero");
        }

        let square_free_factors = self.square_free_factorization();
        let mut factors = Vec::new();
        for (square_free_factor, mult) in square_free_factors {
            for (distinct_degree_factor, deg) in square_free_factor.distinct_degree_factorization()
            {
                let same_degree_factors = distinct_degree_factor.same_degree_factorization(deg);
                factors.extend(same_degree_factors.into_iter().map(|x| (x, mult)));
            }
        }
        factors
    }
}

#[cfg(test)]
mod tests {
    use crate::{Gf2Poly, prop_assert_poly_eq};
    use proptest::prelude::*;


    proptest! {
        #[test]
        fn square_free_factorization(a: Gf2Poly) {
            prop_assume!(!a.is_zero());
            let results = a.square_free_factorization();
            for (i, (first, _)) in results.iter().enumerate() {
                for (second, _) in results.iter().skip(i + 1) {
                    prop_assert_poly_eq!(first.clone().gcd(second.clone()), Gf2Poly::one());
                }
            }
            let mut product = Gf2Poly::one();
            for (res, mul) in results {
                product *= &res.power(mul);
            }
            prop_assert_poly_eq!(product, a);
        }

        #[test]
        fn distinct_degree_factorization(mut a: Gf2Poly) {
            prop_assume!(!a.is_zero());
            // make sure `a` is squarefree
            a /= a.derivative().gcd(a.clone());
            let results = a.distinct_degree_factorization();
            for (i, (first, deg)) in results.iter().enumerate() {
                prop_assert_eq!(first.deg() % deg, 0);
                for (second, _) in results.iter().skip(i + 1) {
                    prop_assert_poly_eq!(first.clone().gcd(second.clone()), Gf2Poly::one());
                }
            }
            let mut product = Gf2Poly::one();
            for (res, _) in results {
                product *= &res;
            }
            prop_assert_poly_eq!(product, a);
        }

        #[test]
        fn factorization(a: Gf2Poly, trace_base: Gf2Poly) {
            prop_assume!(!a.is_zero());
            let factors = a.factor();
            for (i, (first, _)) in factors.iter().enumerate() {
                for (second, _) in factors.iter().skip(i + 1) {
                    prop_assert_poly_eq!(first.clone().gcd(second.clone()), Gf2Poly::one());
                }
            }

            for (factor, _) in factors.iter() {
                let trace = factor.pseudo_trace(factor.deg(), &trace_base);
                prop_assert!(trace.is_constant());
            }

            let mut product = Gf2Poly::one();
            for (res, mul) in factors {
                product *= &res.power(mul);
            }
            prop_assert_poly_eq!(product, a);
        }
    }
}
