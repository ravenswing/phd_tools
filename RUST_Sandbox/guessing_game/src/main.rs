
use rand::Rng; // Rng is a trait!
use std::cmp::Ordering; // Ordering is an enum!
use std::io; // import the io module from standard library

// Main function
fn main() {
    println!("-- GUESS THE NUMBER --");
    
    // Create a random number between 1 - 100 INCLUSIVE
    let secret_number = rand::thread_rng().gen_range(1..=100); 

    // loop creates an infinite loop, until break
    loop {
    
        println!("Please input your guess...");
        // Create mutable guess variable
        let mut guess = String::new();

        // Access stdin function from io module
        io::stdin()
            // readline method on Stdin handle created by stdin function
            .read_line(&mut guess) // &mut indicates a mutable reference 
            // handling of Result (enum) arising from read_line
            .expect("Failed to read line!");
        
        /* Have to change the type of the guess variable!
        * u32 = unsigned interger (32-bit)
        * reassigning variables in this way == SHADOWING
        * .trim()  - String method that eliminates whitespace
        * .parse() - String method that converts it
        *          - also returns Result, hence .expect() is needed
        * guess:   - colon indicates that type will be changed? */
        // let guess: u32 = guess.trim().parse().expect("Please type a number!");

        /* ALTERNATIVELY: incorporate error handlingi
         * Add match to catch the two output of the Result enum
         * _ = catch all i.e. all error codes*/
        let guess: u32 = match guess.trim().parse() {
            Ok(num) => num,
            Err(_) => continue,
        };

        println!("Your guess: {guess}"); // print with placeholders
        
        /* cmp method returns VARIANTS of Ordering enum
        * match is used to compare the output, along different ARMS
        * the arms match PATTERNS (in this case the variants of Ordering. */ 
        match guess.cmp(&secret_number) {
            Ordering::Less => println!("Too small!"),
            Ordering::Greater => println!("Too big!"),
            Ordering::Equal => {
                println!("You win!"); // print message
                break; // exit loop ^^
            }
        }
    }
}
