// Copyright (c) 2017 rust-threshold-secret-sharing developers

/// Euclidean GCD implementation (recursive). The first member of the returned
/// triplet is the GCD of `a` and `b`.
pub fn gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let n = a / b;
        let c = a % b;
        let r = gcd(b, c);
        (r.0, r.2, r.1 - r.2 * n)
    }
}

// TODO see exercise 4.10 in Shoup
pub fn binary_egcd(mut a: i64, mut b: i64) -> (i64, i64, i64) {
    // simple cases
    if a == 0 {
        return (b, 0, 1);
    }
    if b == 0 {
        return (a, 1, 0);
    }

    let mut u = 1;
    let mut v = 0;
    let mut s = 0;
    let mut t = 1;

    // find greatest power r of 2 dividing both a and b
    let mut r = 0;
    while (a | b) & 1 == 0 {
        a >>= 1;
        b >>= 1;
        r += 1;
    }

    let alpha = a;
    let beta = b;

    while a & 1 == 0 {
        a >>= 1;

        if (u | v) & 1 == 0 {
            u >>= 1;
            v >>= 1;
        } else {
            u = (u + beta) >> 1;
            v = (v - alpha) >> 1;
        }
    }

    while a != b {
        if b & 1 == 0 {
            b >>= 1;
            if (s | t) & 1 == 0 {
                s >>= 1;
                t >>= 1;
            } else {
                s = (s + beta) >> 1;
                t = (t - alpha) >> 1;
            }
        } else if b < a {
            a = b;
            b = a;
            u = s;
            v = t;
            s = u;
            t = v;
        } else {
            b -= a;
            s -= u;
            t -= v;
        }
    }

    (a << r, s, t)
}

/// Inverse of `k` in the *Zp* field defined by `prime`.
pub fn mod_inverse(k: i64, prime: i64) -> i64 {
    let k2 = k % prime;
    let r = if k2 < 0 {
        -gcd(prime, -k2).2
    } else {
        gcd(prime, k2).2
    };
    (prime + r) % prime
}
// pub fn mod_inverse(k: i64, prime: i64) -> i64 {
//     let k2 = k % prime;
//     let r = if k2 < 0 {
//         -binary_egcd(prime, -k2).2
//     } else {
//         binary_egcd(prime, k2).2
//     };
//     (prime + r) % prime
// }

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(12, 16), (4, -1, 1));
    }

    #[test]
    pub fn test_binary_egcd() {
        assert_eq!(binary_egcd(10, 4), (2, 1, -2));
    }

    #[test]
    fn test_mod_inverse() {
        assert_eq!(mod_inverse(3, 7), 5);
    }
}
