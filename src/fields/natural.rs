
use rand;
use std::borrow::Borrow;

use ::fields::Field;
use ::fields::PrimeField;
use ::fields::{Encode, Decode};
use ::numtheory::{mod_pow, mod_inverse};

#[derive(Clone,Debug,PartialEq)]
pub struct NaturalPrimeField<T>(pub T);

impl Field for NaturalPrimeField<i64> 
{
    type E = i64;
    
    fn zero(&self) -> Self::E {
        0_i64
    }
    
    fn one(&self) -> Self::E {
        1_i64
    }
    
    fn add<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E {
        (a.borrow() + b.borrow()) % self.0
    }
    
    fn sub<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E {
        (a.borrow() - b.borrow()) % self.0
    }
    
    fn mul<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E {
        (a.borrow() * b.borrow()) % self.0
    }
    
    fn pow<A: Borrow<Self::E>>(&self, a: A, e: u32) -> Self::E {
        mod_pow(*a.borrow(), e, self.0)
    }
    
    fn inv<A: Borrow<Self::E>>(&self, a: A) -> Self::E {
        mod_inverse(*a.borrow(), self.0)
    }
    
    fn eq<L: Borrow<Self::E>, R: Borrow<Self::E>>(&self, lhs: L, rhs: R) -> bool {
        (lhs.borrow() % self.0) == (rhs.borrow() % self.0)
    }
    
    fn sample_with_replacement<R: rand::Rng>(&self, count: usize, rng: &mut R) -> Vec<Self::E> {
        use rand::distributions::Sample;
        let mut range = rand::distributions::range::Range::new(0, self.0 - 1); // TODO why -1?
        (0..count).map(|_| range.sample(rng)).collect()
    }
}

impl PrimeField for NaturalPrimeField<i64> {
    type P = u32;
    
    fn new(prime: Self::P) -> Self {
        NaturalPrimeField(prime as i64)
    }
}

impl Encode<u32> for NaturalPrimeField<i64> {
    fn encode(&self, x: u32) -> Self::E {
        x as i64
    }
}

impl Decode<u32> for NaturalPrimeField<i64> {
    fn decode<E: Borrow<Self::E>>(&self, x: E) -> u32 {
        let x = x.borrow() % self.0; // TODO get rid of this -- should be normalised during computations
        if x >= 0 {
            x as u32
        } else {
            (x + self.0) as u32
        }
    }
}

#[cfg(test)]
all_fields_test!(NaturalPrimeField<i64>);
