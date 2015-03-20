/*************************************************************************
 * Program:
 *    Goodness -- Exploring Iterative Improvement
 *
 * Author:
 *    Gage Peterson (with a lot of help from Weston and Chad)
 *
 * Summary:
 *    A program to experiment with Simulated Annealing.
 *************************************************************************/

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
 * Represents a state for the simulated annealing algorithm.
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
        h = 31 * h + word[i];
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
    
    //output the Java default values and energy
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
    
    //output the Java default values and energy
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
    cout << "I got the double dose of this algorithm both in AI and in this class! It was quite interesting to see how well something so random can work. I spent a long time ";
}

