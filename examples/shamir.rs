// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.
extern crate threshold_secret_sharing as tss;

use tss::*;

fn main() {

    type MyField = NaturalPrimeField<i64>;
    let field = MyField::new(41); // any large enough prime will do

    let ref tss = tss::ShamirSecretSharing {
        threshold: 9,
        share_count: 20,
        field: field,
    };

    let secret = 5;
    let all_shares = tss.share(secret);

    let reconstruct_share_count = 10;
    assert!(reconstruct_share_count >= tss.reconstruct_limit());

    let indices: Vec<usize> = (0..reconstruct_share_count).collect();
    let shares: &[i64] = &all_shares[0..reconstruct_share_count];
    let recovered_secret = tss.reconstruct(&indices, shares);

    println!("The recovered secret is {}", recovered_secret);
    assert_eq!(recovered_secret, secret);

}
