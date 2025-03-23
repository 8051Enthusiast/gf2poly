use crate::{BITS, Gf2Poly, Limb};
use proptest::prelude::*;
#[macro_export]
macro_rules! prop_assert_poly_eq {
    ($lhs:expr, $rhs:expr) => {{
        let lhs: crate::Gf2Poly = $lhs;
        let rhs: crate::Gf2Poly = $rhs;
        ::proptest::prelude::prop_assert_eq!(&lhs, &rhs);
        ::proptest::prelude::prop_assert!(lhs.is_normalized());
        ::proptest::prelude::prop_assert!(rhs.is_normalized());
    }};
}

fn exp_dist(expected_val: f64) -> impl Strategy<Value = u64> {
    (0.0f64..1.0f64)
        .prop_filter("u must be nonzero", |&u| u > 0.0)
        .prop_map(move |u| (-u.ln() * expected_val).floor().max(0.0) as u64)
}

pub fn gf2poly_deg(expected_deg: f64) -> BoxedStrategy<Gf2Poly> {
    exp_dist(expected_deg)
        .prop_flat_map(|deg| {
            let non_leading_limbs = (deg / BITS) as usize;
            let last_bit = (deg % BITS) as usize;
            let lower_mask = (last_bit > 0)
                .then_some(((1 as Limb) << last_bit) - 1)
                .unwrap_or(1);

            (
                prop::collection::vec(any::<Limb>(), non_leading_limbs),
                any::<Limb>(),
            )
                .prop_map(move |(mut limbs, candidate_last)| {
                    let lower_bits = candidate_last & lower_mask;
                    // make sure we also generate 0
                    let last_limb = lower_bits | (((deg != 0) as Limb) << last_bit);
                    if last_limb != 0 {
                        limbs.push(last_limb);
                    }
                    Gf2Poly {
                        deg,
                        limbs: limbs.into(),
                    }
                })
        })
        .boxed()
}

impl Arbitrary for Gf2Poly {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with((): Self::Parameters) -> Self::Strategy {
        gf2poly_deg(64.0)
    }
}
