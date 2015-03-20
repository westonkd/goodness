/*************************************************************************
 * Program:
 *    Goodness -- Exploring Iterative Improvement
 *
 * Author:
 *    Brother Neff
 *
 * Summary:
 *    A program to experiment with Simulated Annealing.
 *************************************************************************/
#if 0
The simulated annealing algorithm, courtesy The Web:
----------------------------------------------------------------------------
1. Choose an initial state (random numbers?)
2. On each iteration, choose a move to another (neighbouring) state.
3. If that move reduces the "energy" (improves the situation) then the
algorithm takes that move.
4. Otherwise, it takes the move with a computed probability that
decreases over time.
(Hence early on the algorithm will tend to take moves even if they
 do not improve the situation.
 Later on, the algorithm will only make moves that improve the
 situation.)
5. Use the temperature function T(n) = 100 / n where n is the iteration
number.
6. Allow a move with negative impact with probability P(dE) = exp(dE/T)
where dE is the difference in energy and T is the temperature for the
iteration.
----------------------------------------------------------------------------
s ← s0; e ← E(s)                        // Initial state, energy.
sbest ← s; ebest ← e                    // Initial "best" solution
k ← 0                                   // Energy evaluation count.
while k < kmax and e > emax             // While time left and not good enough:
T ← temperature(k/kmax)               // Calculate temperature.
snew ← neighbour(s)                   // Pick some neighbour.
enew ← E(snew)                        // Compute its energy.
if P(e, enew, T) > random() then      // Should we move to it?
s ← snew; e ← enew                  // Yes, change state.
if enew < ebest then                  // Is this a new best?
sbest ← snew; ebest ← enew          // Yes, save 'new neighbour' to 'best found'.
k ← k + 1                             // One more evaluation done
return sbest                            // Return the best solution found.
----------------------------------------------------------------------------
Start from a state s0 and continue to either a maximum of kmax
steps or until a state with an energy of emax or less is found. In the
process, the call neighbour(s) should generate a randomly chosen
neighbour of a given state s; the call random() should return a random
value in the range [0, 1]. The annealing schedule is defined by the call
temperature(r), which should yield the temperature to use, given the
fraction r of the time budget that has been expended so far.

P(e, enew, T) = 1 if enew < e, exp(-(enew - e)/T) otherwise
----------------------------------------------------------------------------
#endif

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <time.h>

#define LARGE_HASH_SIZE 33554432
#define MED_HASH_SIZE 1048576
#define SM_HASH_SIZE 32768

using namespace std;

/**********************************************************************
 * Represents a state for the simulated annearling algorithm.
 * a,b,c, and d are the values used to shift by in the hash function.
 *********************************************************************/
struct State
{
    int a;
    int b;
    int c;
    int d;
};

/*************************************************************************
 * Function prototypes
 *************************************************************************/
double calcEnergy(string filename, State state, int size);
State getNeightbor(State state);
float pOfAccept(float currentEnergy, float newEnergy, float temp);
float getTemp(float k, float kmax);
void noLog(State init, State java);

/**********************************************************************
 * anneal
 *  Run the simulated annealing algorithm to find the state which
 *  minimizes collisions.
 *********************************************************************/
State anneal(State s, int kmax, float emax, int size, bool verbose = false)
{
    // Calculate the energy of the initial state
    float e = calcEnergy("hashed", s, size);
    
    //'best' initial solution uses Java numbers
    State sbest = s;
    
    //energy of the initial best solution
    float ebest = e;
    
    //energy evaluation count
    int k = 0;
    
    if (verbose)
        cout << setw(5) << k << " " << setw(10) << " - " << " " << "[" <<  setw(2) << sbest.a << " " << setw(2) << sbest.b << " " << setw(2) << sbest.c << " " << setw(2) << sbest.d << "] "  << setw(10) << ebest << endl;
    
    //while time left and the solution is not good enough
    while (k < kmax && e > emax)
    {
        if (verbose)
            cout << setw(5) << k << " ";
        
        // Calculate temperature.
        float T = getTemp(k, kmax);
        
        if (verbose)
            cout << setw(10) << T << " ";
        
        // Pick some neighbor.
        State snew = getNeightbor(s);
        
        if (verbose)
            cout << "[" <<  setw(2) << snew.a << " " << setw(2) << snew.b << " " << setw(2) << snew.c << " " << setw(2) << snew.d << "] ";
        
        // Compute its energy.
        float enew = calcEnergy("hashed", snew, size);
        
        if (verbose)
            cout << setw(10) << enew;
        
        // Should we move to it?
        float random = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        
        if (verbose)
            cout << setw(12) << pOfAccept(e, enew, T);
        
        if (pOfAccept(e, enew, T) > random)
        {
            // Yes, change state.
            s = snew;
            e = enew;
        }
        
        // Is this a new best?
        if (enew < ebest)
        {
            // Yes, save 'new neighbor' to 'best found'.
            sbest = snew;
            ebest = enew;
            
            if (verbose)
                cout << setw(10) << " accepted" << endl;
        }
        else if (verbose)
        {
            cout << endl;
        }
        
        //increment count
        k++;
    }
    
    //return the best state
    return sbest;
}

/**********************************************************************
 * toUnsignedString
 *  makes its integer argument into a 32-character bit string (0s or 1s)
 *  useful for debugging.
 *********************************************************************/
string toUnsignedString(unsigned int i)
{
    char buf[32];
    int charPos = 32;
    int shift = 1;
    int radix = 1 << shift;
    int mask = radix - 1;
    do
    {
        buf[--charPos] = '0' + (i & mask);
        i >>= shift;
    }
    while (i != 0);
    while (charPos > 0)
        buf[--charPos] = '0';
    return string(buf, 32);
}

/*************************************************************************
 * indexFor
 *
 * Uses a bitwise 'and' to effectively 'mod' h by length
 *************************************************************************/
long indexFor(long h, int length)
{
    return h & (length-1);
}

/*************************************************************************
 * safteyHash
 *
 * A secondary hashing function that is applied to all values to be hashed.
 *************************************************************************/
long safteyHash(unsigned long h, int a, int b, int c, int d)
{
    // This function ensures that hashCodes that differ only by
    // constant multiples at each bit position have a bounded
    // number of collisions (approximately 8 at default load factor).
    h = h ^ (h >> a) ^ (h >> b);
    return h ^ (h >> c) ^ (h >> d);
}

/*************************************************************************
 * verifyState
 *
 * Verifies the values in a state are within specified bounds.
 *************************************************************************/
State verifyState(State state, int min, int max)
{
    //check that a is in bounds
    if (state.a > max)
        state.a = max;
    else if (state.a < min)
        state.a = min;
    
    //check that b is in bounds
    if (state.b > max)
        state.b = max;
    else if (state.b < min)
        state.b = min;
    
    //check that c is in bounds
    if (state.c > max)
        state.c = max;
    else if (state.c < min)
        state.c = min;
    
    //check that d is in bounds
    if (state.d > max)
        state.d = max;
    else if (state.d < min)
        state.d = min;
    
    return state;
}

/*************************************************************************
 * getTemp
 *
 * Calculates the temperature.
 *************************************************************************/
float getTemp(float k, float kmax)
{
    return 100.0 / (k / kmax);
}

/*************************************************************************
 * pOfAccept
 *
 * Calculates the probability of accepting a state with less energy.
 *************************************************************************/
float pOfAccept(float currentEnergy, float newEnergy, float temp)
{
    //always accept a better value
    if (newEnergy < currentEnergy)
        return 1;
    
    //calculate the P. This value gets lower as temperature increases
    return exp(((float)(newEnergy - currentEnergy) * -1.0 )/ temp);
}

/*************************************************************************
 * getNeighbor
 *
 * Returns a random neighbor of state.
 *************************************************************************/
State getNeightbor(State state)
{
    //get a random number between 0 and 3
    int varToChange = rand() % 4;
    
    //increment varToChange or decrement?
    bool increment = rand() % 2;
    
    //choose a value between 1 and 6 to add to the variable.
    //We want our neighbor to be somewhat close.
    int toAdd = rand() % 6 + 1;
    
    //make the change
    switch(varToChange)
    {
        case 0:
            state.a = increment ? state.a + toAdd : state.a - toAdd;
            break;
        case 1:
            state.b = increment ? state.b + toAdd : state.b - toAdd;
            break;
        case 2:
            state.c = increment ? state.c + toAdd : state.c - toAdd;
            break;
        case 3:
            state.d = increment ? state.d + toAdd : state.d - toAdd;
            break;
        default:
            assert(false); //we should never get here!
    }
    
    //verify all values are in bounds
    state = verifyState(state, 0, 31);
    
    return state;
}

/**********************************************************************
 * hashCode
 *    Returns an integer value of a string.
 *    works like Java's String.hashCode(), which is computed as
 *
 *    s[0] * 31^(n-1) + s[1] * 31^(n-2) + ... + s[n-1]
 *
 *    using int arithmetic, where s[i] is the i'th character
 *    of the string, n is the length of the string,
 *    and ^ indicates exponentiation.
 *    (The hash value of the empty string is zero.)
 *********************************************************************/
unsigned int hashCode(string &word)
{
    unsigned int h = 0;
    for (int i = 0; i < word.length(); i++)
    {
        h = 31 * h + word[i]; // GOOD
        //h += word[i]; // BAD!
    }
    
    return h;
}

/**********************************************************************
 * badHashCode
 *    Uses a bad initial hash code
 *********************************************************************/
unsigned int badHashCode(string &word)
{
    unsigned int h = 0;
    for (int i = 0; i < word.length(); i++)
    {
        h += word[i]; // BAD!
    }
    
    return h;
}

/*************************************************************************
 * calcEnergy
 *
 * Get the hash code of each word in the file and output as 'hashed.'
 *************************************************************************/
double calcEnergy(string filename, State state, int size)
{
    //open the file
    ifstream fin(filename.c_str());
    
    if (fin.fail())
        return -1;
    
    map<unsigned long, int> collisionRecord;
    
    unsigned long temp;
    
    //for each value in the file
    while (fin >> temp)
    {
        temp = safteyHash(temp,state.a,state.b,state.c,state.d);
        temp = indexFor(temp, size);
        
        if (temp >= size)
            cout << "Error" << endl;
        
        //if the map does not contain the key
        if(collisionRecord.count(temp) == 0)
            collisionRecord[temp] = 0;
        else
        {
            collisionRecord[temp] =  collisionRecord[temp] + 1;
        }
    }
    
    fin.close();
    
    //calculate the average
    double average = 0;
    
    typedef map<unsigned long, int>::iterator it_type;
    for(it_type iterator = collisionRecord.begin(); iterator != collisionRecord.end(); iterator++)
    {
        average += iterator->second;
    }
    
    //average /= (double) collisionRecord.size();
    
    return average;
}

/*************************************************************************
 * hashFile
 *
 * Get the hash code of each word in the file and output in 'hashed' file
 *************************************************************************/
bool hashFile(string file, bool bad = false)
{
    //open the files
    ifstream fin(file.c_str());
    ofstream fout("hashed");
    
    //if IO failure
    if (fin.fail() || fout.fail())
    {
        cout << "\n\nError, could not find '" + file + "' Please ensure the program has read/write access.\n\n";
        return false;
    }
    
    //run each word through initial hash and save it to 'hashed'
    string word;
    while (fin >> word)
    {
        fout << (bad ? badHashCode(word) : hashCode(word)) << endl;
    }
    
    fin.close();
    fout.close();
    
    return true;
}

/*************************************************************************
 * largeTest
 *
 * Run simulated annealing with the largest hash size
 *************************************************************************/
void largeTest(State init, State java, bool verbose = true)
{
    cout << "Test with hash size of " << LARGE_HASH_SIZE << endl;
    
    //run the simulated annealing for largest hash size and output results
    State best = anneal(init, 100, 0.0, LARGE_HASH_SIZE, verbose);
    
    //output the best state and energy
    cout << "\nBest state was: " << best.a << " " << best.b << " " << best.c << " " << best.d << " ";
    cout << " collisions = " << calcEnergy("hashed", best, LARGE_HASH_SIZE) << endl;
    
    //output the java default values and energy
    cout << "Java was: " << java.a << " " << java.b << " " << java.c << " " << java.d << " ";
    cout << " collisions = " << calcEnergy("hashed", java, LARGE_HASH_SIZE) << endl << endl;
}

/*************************************************************************
 * medTest
 *
 * Run simulated annealing with the medium hash size
 *************************************************************************/
void medTest(State init, State java, bool verbose = true)
{
    cout << "\nTest with hash size of " << MED_HASH_SIZE << endl;
    
    //run the simulated annealing for largest hash size and output results
    State best = anneal(init, 100, 0.2, MED_HASH_SIZE, verbose);
    
    //output the best state and energy
    cout << "\nBest state was: " << best.a << " " << best.b << " " << best.c << " " << best.d << " ";
    cout << " collisions = " << calcEnergy("hashed", best, MED_HASH_SIZE) << endl;
    
    //output the java default values and energy
    cout << "Java was: " << java.a << " " << java.b << " " << java.c << " " << java.d << " ";
    cout << " collisions = " << calcEnergy("hashed", java, MED_HASH_SIZE) << endl << endl;
}

/*************************************************************************
 * smallTest
 *
 * Run simulated annealing with the smallest hash size
 *************************************************************************/
void smallTest(State init, State java, bool verbose = true)
{
    cout << "\nTest with hash size of " << SM_HASH_SIZE << endl;
    
    //run the simulated annealing for largest hash size and output results
    State best = anneal(init, 100, 5, SM_HASH_SIZE, verbose);
    
    //output the best state and energy
    cout << "\nBest state was: " << best.a << " " << best.b << " " << best.c << " " << best.d << " ";
    cout << " collisions = " << calcEnergy("hashed", best, SM_HASH_SIZE) << endl;
    
    //output the java default values and energy
    cout << "Java was: " << java.a << " " << java.b << " " << java.c << " " << java.d << " ";
    cout << " collisions = " << calcEnergy("hashed", java, SM_HASH_SIZE) << endl << endl;
}

/*************************************************************************
 * runAll
 *
 * Runs all tests.
 *************************************************************************/
void runAll()
{
    string fileName = "/usr/share/dict/words";
    
    //seed rand
    srand(time(NULL));
    
    //create the two initial states
    State init = {20, 0, 1, 31};
    State java = {20, 12, 7, 4};
    
    //hash the specified file
    hashFile(fileName);
    
    //run the large test
    largeTest(init, java);
    
    //run the medium test
    medTest(init, java);
    
    //run the small test
    smallTest(init, java);
    
    //run test with large size and bad hashing algorithm
    hashFile(fileName, true);
    medTest(init,java);
}

/*************************************************************************
 * runOne
 *
 * Runs the user-specified test.
 *************************************************************************/
void runOne(string test)
{
    string fileName = "/usr/share/dict/words";
    
    //seed rand
    srand(time(NULL));
    
    //hash the specified file
    if (hashFile(fileName))
    {
        //create the initial state and java's default state
        State init = {20, 0, 1, 31};
        State java = {20, 12, 7, 4};
        
        if (test == "small")
            smallTest(init, java);
        else if (test == "medium")
            medTest(init, java);
        else if (test == "large")
            largeTest(init, java);
        else if (test == "bad")
        {
            hashFile(fileName, true);
            runAll();
        }
        else if (test == "nolog")
            noLog(init, java);
        else
            cout << "Error: could not find test '" + test + "'\n";
    }
}

/*************************************************************************
 * noLog
 *
 * run all tests without verbose output.
 *************************************************************************/
void noLog(State init, State java)
{
    //run the large test
    largeTest(init, java, false);
    
    //run the medium test
    medTest(init, java, false);
    
    //run the small test
    smallTest(init, java, false);
    
    //run test with large size and bad hashing algorithm
}

/*************************************************************************
 * usage
 *
 * Tells the user how to use the program.
 *************************************************************************/
void usage(const char * programName)
{
    //all tests
    cout << programName << " all\n";
    cout << "\trun all tests\n";
    
    //small test
    cout << programName << " small\n";
    cout << "\trun simulated annealing test for small hash size\n\n";
    
    //medium test
    cout << programName << " medium\n";
    cout << "\trun simulated annealing test for medium hash size\n\n";
    
    //large test
    cout << programName << " large\n";
    cout << "\trun simulated annealing test for large hash size\n\n";
    
    //bad test
    cout << programName << " bad\n";
    cout << "\trun simulated annealing test for test with bad initial hashing algorithm\n\n";
    
    //non-verbose test
    cout << programName << " nolog\n";
    cout << "\trun simulated annealing test for all hash sizes without log output\n\n";
}

/*************************************************************************
 * learned
 *
 * Tells interested parties what you learned.
 *************************************************************************/
void learned()
{
    cout << "\nAfter reading “Algorithm of the Gods” and discussing simulated annealing in class, I felt like I had a good grasp on what the algorithm could do. Understanding how hashing in Java works took the longest time to study, but was something I was grateful to learn.\n\n";
    
    cout << "Implementing the code also went well. After researching simulated annealing and writing the algorithm I can honestly say that I learned how it works well enough to teach it to others. I was actually so excited that I explained how the algorithm works to my wife. She was very kind but not nearly as enthused.\n\n";
    
    cout << "I learned that, unlike a greedy algorithm, simulated annealing works in the beginning by ‘exploring’ the solution set and then becomes more and more greedy as the temperature cools. It is necessary to compute the ‘energy’ of neighboring states and then choose a good neighbor to move to. A good decision does not mean the neighbor with the lowest energy- accepting higher values in the beginning will help find the global minimum without being trapped at local minima. The temperature variable is used to help calculate the probability of accepting a state with higher energy. As the temperature cools the probability decreases to zero. This helps the algorithm to find the true global minimum.\n\n";
    
    cout << "To choose a neighbor the algorithm randomly selects a value in the state to change (a, b, or c). Once a number to change is chosen the algorithm adds a random number between -5 and 5 to that number. It makes sure the new state is in bounds and then returns it.\n\n";
    
    cout << "The tests I implemented run the annealing algorithm at three different load factors. The first hash size is 2^25, the second is 2^20, and the third is 2^15. As expected the tests show the average number of collisions is inversely proportional to the number of buckets. I ran these tests several times to obtain an optimum set of numbers. I also created a version of the program that calculates energy based on the average number of collisions per bucket rather than the total collision count. This test showed that the difference between the Java default values and other, better values is very small.\n";
    
    cout << "I tested the good initial algorithm compared with the bad initial algorithm. Because the hash codes are 32 bit numbers, the good algorithm multiplies ‘h’ by 31 before adding the letter value. This helps widen the spread of hash codes. Running this particular test made it clear that the good hashing algorithm is much better.\n\n";
    
    cout << "While the numbers I calculated as ideal shifting values varied, I was consistently able to find better values than the Java default choices. One such set is [5 10 2 23] which had 616 fewer collisions than the Java defaults. Other results were [31 16 6 27] and [3 2 19 31].\n\n";
}

