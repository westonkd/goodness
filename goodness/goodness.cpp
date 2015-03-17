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

#define HASH_SIZE 1048576

using namespace std;

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
 * hash
 *
 *
 *************************************************************************/
int safteyHash(unsigned int h)
{
    // This function ensures that hashCodes that differ only by
    // constant multiples at each bit position have a bounded
    // number of collisions (approximately 8 at default load factor).
    h = h ^ (h >> 20) ^ (h >> 12);
    return h ^ (h >> 7) ^ (h >> 4);
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
   }
    
   return h % HASH_SIZE;
}

/*************************************************************************
 * calcEnergy
 *
 * Get the hash code of each word in the file and output as 'hashed'
 *************************************************************************/
double calcEnergy(string filename)
{
    //open the file
    ifstream fin(filename.c_str());
    
    if (fin.fail())
        return -1;

    map<int,int> collisionRecord;
    
    int temp;
    
    //for each value in the file
    while (fin >> temp)
    {
        temp = safteyHash(temp);
        
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
    
    typedef map<int, int>::iterator it_type;
    for(it_type iterator = collisionRecord.begin(); iterator != collisionRecord.end(); iterator++)
    {
        average += iterator->second;
    }
    
    average /= (double) collisionRecord.size();

    //return the average collisions
    return average;
}

/*************************************************************************
 * hashFile
 *
 * Get the hash code of each word in the file and output as 'hashed'
 *************************************************************************/
void hashFile(string file)
{
    ifstream fin(file.c_str());
    ofstream fout("hashed");
    
    if (fin.fail() || fout.fail())
    {
        cerr << "Error reading file";
        return;
    }
    
    string word;
    while (fin >> word)
    {
        fout << hashCode(word) << endl;
    }
    
    fin.close();
    fout.close();
}

/*************************************************************************
 * runOne
 *
 * Runs one test.
 *************************************************************************/
void runOne(string test)
{
}

/*************************************************************************
 * runAll
 *
 * Runs all tests.
 *************************************************************************/
void runAll()
{
    hashFile("words");
    cout << "Average number of collisions: " << calcEnergy("hashed") << endl;
}

/*************************************************************************
 * usage
 *
 * Tells the user how to use the program.
 *************************************************************************/
void usage(const char * programName)
{
}

/*************************************************************************
 * learned
 *
 * Tells interested parties what you learned.
 *************************************************************************/
void learned()
{
}
