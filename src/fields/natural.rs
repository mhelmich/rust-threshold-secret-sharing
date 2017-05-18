
use rand;

use ::fields::Field;
use ::fields::PrimeField;
use ::fields::{Encode, Decode};
use ::numtheory::{mod_pow, mod_inverse};

#[derive(Debug, PartialEq)]
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
    
    fn add(&self, a: Self::E, b: Self::E) -> Self::E {
        (a + b) % self.0
    }
    
    fn sub(&self, a: Self::E, b: Self::E) -> Self::E {
        (a - b) % self.0
    }
    
    fn mul(&self, a: Self::E, b: Self::E) -> Self::E {
        (a * b) % self.0
    }
    
    fn pow(&self, a: Self::E, e: u32) -> Self::E {
        mod_pow(a, e, self.0)
    }
    
    fn inv(&self, a: Self::E) -> Self::E {
        mod_inverse(a, self.0)
    }
    
    fn eq(&self, lhs: Self::E, rhs: Self::E) -> bool {
        (lhs % self.0) == (rhs % self.0)
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
    fn decode(&self, x: Self::E) -> u32 {
        let x = x % self.0; // TODO get rid of this -- should be normalised during computations
        if x >= 0 {
            x as u32
        } else {
            (x + self.0) as u32
        }
    }
}

#[cfg(test)]
all_fields_test!(NaturalPrimeField<i64>);
