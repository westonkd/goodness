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

#define LARGE_HASH_SIZE 1048576
#define MED_HASH_SIZE 32768
#define SM_HASH_SIZE 1024

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

/**********************************************************************
 * anneal
 *  Run the simulated annealing algorithm to find the state wich 
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
        
        // Pick some neighbour.
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
            // Yes, save 'new neighbour' to 'best found'.
            sbest = snew;
            ebest = enew;
            cout << setw(10) << " accepted" << endl;
        }
        else
        {
            cout << endl;
        }
        
        k++;
    }

    //return the best state
    return sbest;
}

/**********************************************************************
 * toUnsignedString
 *  makes its integer argument into a 32-character bitstring (0s or 1s)
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
 * TODO
 *
 *************************************************************************/
long indexFor(long h, int length)
{
    return h & (length-1);
}

/*************************************************************************
 * safteyHash
 * A secondary hashing function that is applied to all values to be hashed.
 *
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
 * verifies the values in a state are within specified bounds. If 
 * a value is over or under max or min, the number wraps around.
 *************************************************************************/
State verifyState(State state, int min, int max)
{
//    //check that a is in bounds
//    if (state.a > max)
//        state.a = state.a % max;
//    else if (state.a < min)
//        state.a = state.a + max;
//    
//    //check that b is in bounds
//    if (state.b > max)
//        state.b = state.b % max;
//    else if (state.b < min)
//        state.b = state.b + max;
//    
//    //check that c is in bounds
//    if (state.c > max)
//        state.c = state.c % max;
//    else if (state.c < min)
//        state.c = state.c + max;
//    
//    //check that d is in bounds
//    if (state.d > max)
//        state.d = state.d % max;
//    else if (state.d < min)
//        state.d = state.d + max;

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
//
    return state;
}

/*************************************************************************
 * getTemp
 *
 * calculates the temperature.
 *************************************************************************/
float getTemp(float k, float kmax)
{
    return 100.0 / (k / kmax);
}

/*************************************************************************
 * pOfAccept
 *
 * calculates the probability of accepting a state with less energy.
 *************************************************************************/
float pOfAccept(float currentEnergy, float newEnergy, float temp)
{
    //always accept a better value
    if (newEnergy < currentEnergy)
        return 1;
    
    //calculate the P. This value gets lower as temperature increases
    return exp((newEnergy - currentEnergy) * -1 / temp);
}

/*************************************************************************
 * getNeighbor
 *
 * returns a random neightbor of state
 *************************************************************************/
State getNeightbor(State state)
{
    //get a random number between 0 and 3
    int varToChange = rand() % 4;
    
    //increment varToChange or decrement?
    bool increment = rand() % 2;
    
    //choose a value between 1 and 6 to add to the variable.
    //We want our niehgbor to be somewhat close.
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
 *    returns an integer value of a string.
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

/*************************************************************************
 * calcEnergy
 *
 * Get the hash code of each word in the file and output as 'hashed'
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
    
    //cout << "using " << average << " / " << collisionRecord.size() << " = ";
    average /= (double) collisionRecord.size();

    //return the average collisions
    return average;
}

/*************************************************************************
 * hashFile
 *
 * Get the hash code of each word in the file and output in 'hashed' file
 *************************************************************************/
bool hashFile(string file)
{
    //open the files
    ifstream fin(file.c_str());
    ofstream fout("hashed");
    
    //if IO failure
    if (fin.fail() || fout.fail())
    {
        cout << "Error, could not find '" + file + "'\n";
        return false;
    }
    
    //run each word through initial hash and save it to 'hashed'
    string word;
    while (fin >> word)
    {
        fout << hashCode(word) << endl;
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
void largeTest(State init, State java)
{
    cout << "Test with hash size of " << LARGE_HASH_SIZE << endl;
    
    //run the simulated annealing for largest hash size and output results
    State best = anneal(init, 100, 0.1, LARGE_HASH_SIZE, true);
    
    //output the best state and energy
    cout << "\nBest state was: " << best.a << " " << best.b << " " << best.c << " " << best.d << " ";
    cout << calcEnergy("hashed", best, LARGE_HASH_SIZE) << endl;
    
    //output the java default values and energy
    cout << "Java was: " << java.a << " " << java.b << " " << java.c << " " << java.d << " ";
    cout << calcEnergy("hashed", java, LARGE_HASH_SIZE);
}

/*************************************************************************
 * medTest
 *
 * Run simulated annealing with the medium hash size
 *************************************************************************/
void medTest(State init, State java)
{
    cout << "\nTest with hash size of " << MED_HASH_SIZE << endl;
    
    //run the simulated annealing for largest hash size and output results
    State best = anneal(init, 100, 0.1, MED_HASH_SIZE, true);
    
    //output the best state and energy
    cout << "\nBest state was: " << best.a << " " << best.b << " " << best.c << " " << best.d << " ";
    cout << calcEnergy("hashed", best, MED_HASH_SIZE) << endl;
    
    //output the java default values and energy
    cout << "Java was: " << java.a << " " << java.b << " " << java.c << " " << java.d << " ";
    cout << calcEnergy("hashed", java, MED_HASH_SIZE);
}

/*************************************************************************
 * smallTest
 *
 * Run simulated annealing with the smallest hash size
 *************************************************************************/
void smallTest(State init, State java)
{
    cout << "\nTest with hash size of " << SM_HASH_SIZE << endl;
    
    //run the simulated annealing for largest hash size and output results
    State best = anneal(init, 100, 0.1, SM_HASH_SIZE, true);
    
    //output the best state and energy
    cout << "\nBest state was: " << best.a << " " << best.b << " " << best.c << " " << best.d << " ";
    cout << calcEnergy("hashed", best, SM_HASH_SIZE) << endl;
    
    //output the java default values and energy
    cout << "Java was: " << java.a << " " << java.b << " " << java.c << " " << java.d << " ";
    cout << calcEnergy("hashed", java, SM_HASH_SIZE);
}

/*************************************************************************
 * runOne
 *
 * Runs one using the largest hash size
 *************************************************************************/
void runOne(string test, string fileName = "words")
{
    //seed rand
    srand(time(NULL));
    
    //hash the specified file
    hashFile(fileName);
    
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
        cout << "bad";
    else
        cout << "Error: could not find test '" + test + "'\n";
}

/*************************************************************************
 * runAll
 *
 * Runs all tests.
 *************************************************************************/
void runAll(string fileName = "words")
{
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
    cout << "\trun simulated annealing test for test with bad initial hashing algorithm\n";
}

/*************************************************************************
 * learned
 *
 * Tells interested parties what you learned.
 *************************************************************************/
void learned()
{
}
