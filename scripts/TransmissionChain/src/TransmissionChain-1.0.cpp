//============================================================================
// Name        : TransmissionChain-1.0.cpp
// Version     : 1.0
// Description : 1.0 Implement a neutral transmission chain.
//               Each infection is described by an exponential birth-death process.
//               Transmission occurs randomly once the population has reached
//               a certain size.
//============================================================================

#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <time.h>
#include <iomanip>
#include <algorithm>
#include <map>
#include <stdio.h>
#include <vector>
using namespace std;

//=========================================//
// Model parameters                        //
//=========================================//

double MU=1; // genome-wide mutation rate

//=========================================//
// Global variables, structs, functions    //
//=========================================//

// Data structures

// Mutation type.
// Contains an identifier number
// and a fitness effect.
struct mutation_t{
	long long mutID;
	double mutEffect;
};
// Variant type.
// Contains an identifier number
// and a fitness effect.
struct variant_t{
	long long mutID;
	double mutEffect=0;
	double mutFreq=0;
	long long mutNumIndividuals=0;
};
// Clone type.
// By default, begins with one individual
// and no mutations.
struct clone_t{
	vector <mutation_t> mutations;
	double fitness=1;
	long long numIndividuals=1;
};
// Population type.
// By default, begins with no clones
// and no individuals.
struct population_t{
	vector <clone_t> clones;
	long long numClones=0;
	long long numIndividuals=0;
	double totFitness=0;
	long long iMutations=0;
};

// Random number generators
auto const seed = random_device()();
mt19937 generator(seed);
uniform_real_distribution<double> uniform(0.,1.);


// Helper functions
void initiatePopulation(population_t*);
void exportPopulation(population_t*, string);
void exportVariants(population_t*, string);
void birth(population_t*);
void death(population_t*);
void childClone(population_t*, clone_t*);
double mutationFitness(void);

int main() {

	// Allocate memory for a new population.
	population_t *population = new population_t;

	// Print population size as a sanity check.
	printf ("Population size: %lld\n", population->numIndividuals);

	// Initiate the population by adding a clone.
	initiatePopulation(population);

	// Print population size again after adding a clone.
	printf ("Population size: %lld\n", population->numIndividuals);

	// Print the fitness of the first clone.
	printf ("Fitness of first clone: %f\n", population->clones[0].fitness);

	// Export the population and the variants it contains.
	birth(population);
	birth(population);
	birth(population);
	birth(population);
	birth(population);
	exportPopulation(population, "test.pop");
	exportVariants(population, "test.var");

	// Clear allocated memory.
	delete population;

	printf ("Hello world!");
	return 0;
}

// initiatePopulation
// Given a new instance of the population class,
// populate it with one clone with a single individual.
// This clone carries no mutations.
void initiatePopulation(population_t* population){
	// Create a new clone and add it to the list.
	clone_t *newClone = new clone_t;
	population->clones.push_back(*newClone);
	population->numClones++;
	population->numIndividuals += newClone->numIndividuals;
	population->totFitness = 1;
	population->iMutations = 0;
	delete(newClone);
	return;
};

// exportPopulation
// Given a population and a file name,
// write the clones in a population out to the specified file.
void exportPopulation(population_t* population, string file){
	FILE * out;
	out = fopen(file.c_str(),"w");
	fprintf(out, "POPULATION\n");
	fprintf(out, "NUMCLONES\n%lld\n", population->numClones);
	fprintf(out, "NUMINDIVIDUALS\n%lld\n", population->numIndividuals);
	fprintf(out, "TOTFITNESS\n%f\n", population->totFitness);
	fprintf(out, "IMUTATIONS\n%lld\n", population->iMutations);
	// Iterate through each clone and print clone properties.
	for(vector <clone_t>::iterator iClone=population->clones.begin();
			iClone!=population->clones.end(); ++iClone){
		fprintf(out, "CLONE\n");
		fprintf(out, "NUMINDIVIDUALS\n%lld\n", iClone->numIndividuals);
		fprintf(out, "FITNESS\n%f\n", iClone->fitness);
		fprintf(out, "MUTATIONS\n");
		// Iterate through and print the mutations
		// associated with each clone.
		for(vector <mutation_t>::iterator iMut=iClone->mutations.begin();
				iMut!=iClone->mutations.end(); ++iMut){
			fprintf(out, "%lld\t%f\n", iMut->mutID, iMut->mutEffect);
		}
		fprintf(out, "END CLONE\n");
	}
	fclose(out);
	return;
};

// exportVariants
// Given a population and a file name,
// tally the total frequency of each mutation
// and export each mutation, its frequency, and fitness
// to the specified file.
void exportVariants(population_t* population, string file){
	FILE * out;
	out = fopen(file.c_str(),"w");
	vector <variant_t> variants;
	// Populate the vector of variants by going through the population.
	// Iterate through each clone in the population.
	for(vector <clone_t>::iterator iClone=population->clones.begin();
				iClone!=population->clones.end(); ++iClone){
		// Iterate through each mutation in the clone.
		for(vector <mutation_t>::iterator iMut=iClone->mutations.begin();
				iMut!=iClone->mutations.end(); ++iMut){
			// Check the list of variants.
			// If the mutation is already present,
			// then add the number of individuals in the clone.
			bool mutFound=false;
			for(vector <variant_t>::iterator iVar=variants.begin();
					iVar!=variants.end(); ++iVar){
				if(iVar->mutID==iMut->mutID){
					iVar->mutNumIndividuals += iClone->numIndividuals;
					mutFound=true;
				}
			}
			if(!mutFound){
				// Otherwise, add a new mutation to the list of variants.
				variant_t *newVariant = new variant_t;
				newVariant->mutID = iMut->mutID;
				newVariant->mutEffect = iMut->mutEffect;
				newVariant->mutNumIndividuals = iClone->numIndividuals;
				variants.push_back(*newVariant);
				delete(newVariant);
			}
		}
	}
	// Print out the vector of variants.
	for(vector <variant_t>::iterator iVar=variants.begin();
			iVar!=variants.end(); ++iVar){
		iVar->mutFreq = (double) iVar->mutNumIndividuals / (double) population->numIndividuals;
		fprintf(out, "%lld\t%f\t%f\n", iVar->mutID, iVar->mutFreq, iVar->mutEffect);
	}
	fclose(out);
	return;
};

// birth
// Given a population, choose an individual to give birth,
// taking its fitness into account.
// Determine whether a mutation occurs to form a new clone.
// If so, then generate a new mutation and add it to a new clone.
// Add that clone to the population and update the population.
// If not, then increment the clone that gives birth
// and update the population.
void birth(population_t* population){
	// Randomly select an individual based on its fitness.
	// double fitness = uniform(generator) * population->totFitness;
	double fitness = population->totFitness;
	// Search for the chosen clone by iterating through clones
	// and tallying the cumulative fitness.
	double iFitness = 0;
	for(vector <clone_t>::iterator iClone=population->clones.begin();
				iClone!=population->clones.end(); ++iClone){
		iFitness += (double) iClone->numIndividuals * iClone->fitness;
		if(iFitness >= fitness){
			// Determine whether a mutation occurs.
			bool mutate=false;
			if(MU>uniform(generator)){
				mutate=true;
			}
			// If no mutations occur, then increment the clone
			// and update the population tracking parameters.
			if(!mutate){
				iClone->numIndividuals += 1;
				population->numIndividuals += 1;
				population->totFitness += iClone->fitness;
			}
			// If a mutation occurs, then create a new clone
			// based on the parent clone.
			// Add it to the population and update tracking parameters.
			else{
				childClone(population, &(*iClone));
			}
			break;
		}
	}
	return;
};

// death
// Given a population, choose a random individual to die.
// Identify the clone of which it is a part,
// and decrement the size of that clone.
// Remove the clone if it is no longer populated.
void death(population_t* population){
	// Randomly select an individual.
	long long individual =
			(long long) (uniform(generator) * (double) population->numIndividuals);
	// Search for the chosen individual among the populated clones.
	long long iNumIndividuals = 0;
	for(vector <clone_t>::iterator iClone=population->clones.begin();
			iClone!=population->clones.end(); ++iClone){
		iNumIndividuals += iClone->numIndividuals;
		// Decrement the number of individuals in the appropriate clone.
		// Also adjust the population tracking parameters.
		if(iNumIndividuals > individual){
			iClone->numIndividuals -= 1;
			population->numIndividuals-=1;
			population->totFitness-=iClone->fitness;
			// If there are no more individuals left in the clone,
			// then remove the clone.
			if(iClone->numIndividuals==0){
				population->clones.erase(iClone);
			}
		}
		break;
	}

	return;
};

// childClone
// Given a population and a parent clone that produces a new clone through mutation,
// create a copy of the parent clone,
// generate the mutation that differentiates the clones,
// add the mutation to the child clone, and return the child clone.
void childClone(population_t* population, clone_t* parent){
	// Generate the new clone by copying parameters from the parent clone.
	// Copy the parent fitness to the child.
	// The number of individuals in the child clone is 1 by default.
	clone_t *child = new clone_t;
	child->fitness = parent->fitness;
	// Iterate through the mutations in the parent and copy them to the child.
	for(vector <mutation_t>::iterator iMut=parent->mutations.begin();
			iMut!=parent->mutations.end(); ++iMut){
		child->mutations.push_back(*iMut);
	}
	// Generate a new mutation and add it to the child.
	mutation_t *mutation = new mutation_t;
	mutation->mutID = population->iMutations;
	mutation->mutEffect = mutationFitness();
	child->mutations.push_back(*mutation);
	population->iMutations++;
	child->fitness += mutation->mutEffect;
	// Add the child clone to the population.
	population->clones.push_back(*child);
	// Update population tracking parameters.
	population->totFitness += child->fitness;
	population->numIndividuals += 1;
	population->numClones += 1;
	delete(mutation);
	delete(child);
	return;
};

// mutationFitness
// Returns a mutation fitness.
// Currently, all mutations have selective coefficient 0.
double mutationFitness(void){
	return (double) 0;
}
