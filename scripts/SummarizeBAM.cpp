//============================================================================
// Name        : SummarizeBAM.cpp
// Version     : 1.21
// Description : 1.21 Modify types to account for very high coverage.
//               Ignore secondary read alignments.
//           1.2 Fix base position counter to account for read orientation.
//               Modify file format to remove header row and include genome position.
//           1.1 Add read trimming functionality.
//           1.0 Given a BAM file and reference sequence,
//               tally the base counts at each position
//               as well as the average base qualities and read positions.
//               Output a summary of each base at each position
//               as well as the consensus FASTA for the alignment file.
//============================================================================

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <map>
#include <cstring>

using namespace std;

// RUN PARAMETERS
string SAM="";
string REFFASTA="";
string OUTFILE="";
string OUTFASTA="";
int BASEQTHRESHOLD=20;
int MAPQTHRESHOLD=20;
int LEFTTRIM=0;
int RIGHTTRIM=0;

struct PositionBase_t{
	long long Count;
	long long TotalQuality;
	long long TotalReadPosition; // 1-indexed read position
};

// FUNCTIONS
int ArgsParse(int argc, char *argv[]);
void PrintUsage();
void PrintParameters();
void SetDebug();
int ReadMultiFasta(string filename,
		vector<string> *sequencenames,
		vector<string> *sequences);
vector<string> StringSplit(string s, char c);
string ExpandCIGAR(string cigar);

int main(int argc, char *argv[]) {

	//==================================================
	// Parse command-line arguments.
	//==================================================

	if(ArgsParse(argc, argv) != 0){
		PrintUsage();
		return 1;
	}

	PrintParameters();


	//==================================================
	// Read in reference sequence.
	//==================================================

	printf("Reading reference.\n");
	vector<string> RefNames;
	vector<string> RefSequences;
	if(ReadMultiFasta(REFFASTA,&RefNames, &RefSequences) != 0){
		printf("Error: reference sequence does not exist.\n");
		return 1;
	}

	//==============================================================
	// Create a data structure to store the BAM summary information.
	//==============================================================

	printf("Initializing data structure.\n");

	// Structure is a nested series of maps, called as follows:
	// BAMSummary[string chr][int pos][char base].
	//   (Count/TotalQuality/TotalReadPosition)

	map<string, map<int, map<char, PositionBase_t> > > BAMSummary;
	const char BASES[4]={'A','C','G','T'};

	// Also initiate a data structure to store total coverage at each position.
	map<string, map<int, int> > BAMCoverage;

	// Initialize the data structures.
	for(unsigned int i=0; i<RefSequences.size();i++){
		for(unsigned int j=0; j<RefSequences[i].size();j++){
			BAMCoverage[RefNames[i]][j]=0;
			for(unsigned int k=0; k<sizeof(BASES)/sizeof(char);k++){
				PositionBase_t PositionBase;
				BAMSummary[RefNames[i]][j][BASES[k]]=PositionBase;
				BAMSummary[RefNames[i]][j][BASES[k]].Count=0;
				BAMSummary[RefNames[i]][j][BASES[k]].TotalQuality=0;
				BAMSummary[RefNames[i]][j][BASES[k]].TotalReadPosition=0;
			}
		}
	}

	//==================================================
	// Read in BAM file and tally reads.
	//==================================================

	printf("Reading SAM file.\n");

	// Open the file.
	ifstream fin(SAM.c_str(), ios::in);

	// Tally the number of reads that contain indels and are not counted.
	int NumIndels=0;

	if(fin){

		// Read in the file line by line.
		string line;

		while(getline(fin, line)){

			// Split the tab-delimited line.
			vector<string> fields=StringSplit(line,'\t');

			// Store the information in the appropriate formats.
			// Fields are hard-coded based on BAM file format.
			// Convert sequence positions from one-indexed to zero-indexed.
			int Flag=atoi(fields[1].c_str())-1;
			string Chr=fields[2];
			int StartPos=atoi(fields[3].c_str())-1;
			int MapQ=atoi(fields[4].c_str())-1;
			string Cigar=fields[5];
			int TLen=atoi(fields[8].c_str());
			string Read=fields[9];
			string Quality=fields[10];

			// Exclude reads that did not map based on the CIGAR string.
			if(Cigar=="*"){
				continue;
			}

			// Exclude secondary read alignments,
			// i.e. those with FLAG above 256.
			if(Flag>256){
				continue;
			}

			// Expand the CIGAR string to account for soft clippings.
			string ExpandedCigar=ExpandCIGAR(Cigar);
			// Verify that the expanded CIGAR string matches the read length.
			if(ExpandedCigar.size() != Read.size()){
				printf("CIGAR parsing error.\n");
			}
			// Determine from the CIGAR string whether indels are present.
			// Tally the number of reads that contain indels and are not counted.
			bool indel=false;
			string ForbiddenOperations="IDP";
			for(unsigned int i=0; i<Cigar.size(); i++){
				if(strchr(ForbiddenOperations.c_str(), Cigar[i])!=NULL){
					indel=true;
					NumIndels++;
				}
			}

			// Consider only reads that map to the reference,
			// do not contain indels,
			// and have mapping quality above a certain threshold.
			if(find(RefNames.begin(), RefNames.end(), Chr)!=RefNames.end() &&
					MapQ>MAPQTHRESHOLD && indel==false){

				// Count the number of bases that are soft-clipped
				// from each end of the read.
				int LeftClip=0;
				int RightClip=0;
				for(unsigned int i=0; i<ExpandedCigar.size();i++){
					if(ExpandedCigar[i]=='S'){
						LeftClip++;
					}
					else{
						break;
					}
				}
				for(unsigned int i=0; i<ExpandedCigar.size();i++){
					if(ExpandedCigar[ExpandedCigar.size()-1-i]=='S'){
						RightClip++;
					}
					else{
						break;
					}
				}

				// Iterate along length of read and quality scores.
				// Tally bases that are correctly aligned based on the CIGAR string.
				// Adjust base numbering based on CIGAR string information.
				// Take into account read orientation based on TLen field.
				int NumAligned=0;
				for(unsigned int i=0; i<Read.size();i++){
					if(ExpandedCigar[i]!='I' && ExpandedCigar[i]!='S'){
						int RefPos=StartPos + NumAligned;
						// Tally only bases that exceed the quality threshold
						// and that have not been trimmed.
						if(((int) Quality[i])-33 > BASEQTHRESHOLD &&
								i >= LeftClip + LEFTTRIM &&
								i < Read.size() - RightClip - RIGHTTRIM){
							BAMCoverage[Chr][RefPos]++;
							BAMSummary[Chr][RefPos][Read[i]].Count++;
							BAMSummary[Chr][RefPos][Read[i]].TotalQuality+=
								((int) Quality[i])-33;
							// Tally base position in read,
							// accounting for read orientation.
							if(TLen>=0){
								BAMSummary[Chr][RefPos][Read[i]].TotalReadPosition+=i+1;
							}
							else{
								BAMSummary[Chr][RefPos][Read[i]].TotalReadPosition+=
										Read.size()-(i+1);
							}
						}
						NumAligned++;
					}
				}
			}

		}
	}
	else{
		printf("Error: SAM file does not exist.\n");
		return 1;
	}

	// Close the file.
	fin.close();

	//==============================================================
	// Output summary of base frequencies at each position.
	//==============================================================

	printf("Writing base frequencies.\n");

	ofstream fout(OUTFILE.c_str(), ios::out);
	//fout << "Chr\tPos\tBase\tRefBase\tGenomePos\tCount\tAvgQ\tAvgReadPos" << endl;

	// Iterate through the summary data structure
	// and output the desired values using 1-indexed read positions.
	long long GenomicPosition=0;
	for(unsigned int i=0; i<RefSequences.size();i++){
		for(unsigned int j=0; j<RefSequences[i].size();j++){
			GenomicPosition++;
			for(unsigned int k=0; k<sizeof(BASES)/sizeof(char);k++){
				if(BAMSummary[RefNames[i]][j][BASES[k]].Count > 0){
					fout << RefNames[i] << "\t" <<
							j+1 << "\t" <<
							BASES[k] << "\t" <<
							RefSequences[i][j] << "\t" <<
							GenomicPosition << "\t" <<
							BAMSummary[RefNames[i]][j][BASES[k]].Count <<
							"\t" <<
							(float) BAMSummary[RefNames[i]][j][BASES[k]].TotalQuality/
							BAMSummary[RefNames[i]][j][BASES[k]].Count <<
							"\t" <<
							(float) BAMSummary[RefNames[i]][j][BASES[k]].TotalReadPosition/
							BAMSummary[RefNames[i]][j][BASES[k]].Count<<
							"\t" << endl;
				}
				// For positions with 0 counts, replace the "nan" with 0.
				else{
					fout << RefNames[i] << "\t" <<
						j+1 << "\t" <<
						BASES[k] << "\t" <<
						RefSequences[i][j] << "\t" <<
						GenomicPosition << "\t" <<
						BAMSummary[RefNames[i]][j][BASES[k]].Count <<
						"\t" <<
						"0" <<
						"\t" <<
						"0"<<
						"\t" << endl;
				}
			}
		}
	}

	fout.close();

	//==============================================================
	// Output consensus FASTA file for the alignment.
	//==============================================================

	// Iterate through the summary data structure
	// and determine the consensus base at each position.

	if(OUTFASTA != ""){

		printf("Writing consensus reference.\n");

		ofstream foutf(OUTFASTA.c_str(), ios::out);

		for(unsigned int i=0; i<RefSequences.size();i++){
			foutf << ">" << RefNames[i] << endl;
			for(unsigned int j=0; j<RefSequences[i].size();j++){
				char maxbase='N';
				for(unsigned int k=0; k<sizeof(BASES)/sizeof(char);k++){
					if(BAMSummary[RefNames[i]][j][BASES[k]].Count >
						BAMSummary[RefNames[i]][j][maxbase].Count){
						maxbase=BASES[k];
					}
				}
				foutf << maxbase;
				// Insert a line break in the sequence every 70 bases.
				if((j+1)%70 == 0){
					foutf << endl;
				}
			}
			foutf << endl;
		}

		foutf.close();
	}

	printf("Number of reads containing indels: %d\n", NumIndels);
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	return 0;
}

// ArgsParse
// Parses command-line arguments.
// Returns 1 if any argument conditions are violated.
int ArgsParse(int argc, char *argv[]){

	// If the only argument is debug,
	// set all parameters to the debug state.
	if(argc==2 && strcmp(argv[1],"debug")==0){
		SetDebug();
		return 0;
	}

	// Ensure that there are an even number of arguments,
	// leaving aside the program name.
	// The only exception is if the only argument is debug,
	// which sets all parameters to the debug state.
	if((argc - 1) % 2 != 0){
		printf("Invalid number of arguments.\n");
		return 1;
	}
	// Check the structure of arguments.
	for(int i=1; i<argc; i++){
		// Verify that every other argument is a flag.
		if(i%2 != 0){
			if(argv[i][0] != '-' || strlen(argv[i])!=2){
				printf("Invalid use of argument flags.\n");
				return 1;
			}
		}
	}

	// Parse each pair of arguments.
	for(int i=0; i<(argc-1)/2; i++){

		string flag=argv[2*i+1];
		string arg=argv[2*i+2];

		// Parse the flag string.
		switch(flag[1]){
		// -i input SAM file
		case 'i':
			SAM = arg;
			break;
		// -r reference FASTA
		case 'f':
			REFFASTA = arg;
			break;
		// -o output summary file
		case 'o':
			OUTFILE = arg;
			break;
		// -f FASTA file summary of SAM consensus
		case 's':
			OUTFASTA = arg;
			break;
		// -Q base quality threshold
		case 'Q':
			BASEQTHRESHOLD = atoi(arg.c_str());
			if(BASEQTHRESHOLD > 40 || BASEQTHRESHOLD < 0){
				printf("Invalid -Q base quality threshold.\n");
				return 1;
			}
			break;
		// -q mapping quality threshold
		case 'q':
			MAPQTHRESHOLD = atoi(arg.c_str());
			break;
		// -l left read trim
		case 'l':
			LEFTTRIM = atoi(arg.c_str());
			break;
		// -r right read trim
		case 'r':
			RIGHTTRIM = atoi(arg.c_str());
			break;
		}
	}

	// Check that the required arguments exist.
	if(SAM==""){
		printf("Invalid arguments. Specify SAM file.\n");
		return 1;
	}
	if(REFFASTA==""){
		printf("Invalid arguments. Specify reference sequence.\n");
		return 1;
	}
	if(OUTFILE==""){
		printf("Invalid arguments. Specify output file.\n");
		return 1;
	}
	return 0;
}

// PrintParameters
// When called, prints the parameters for the run.
void PrintParameters(){
	cout << "RUN PARAMETERS" << endl;
	cout << "input file: " << SAM << endl;
	cout << "reference: " << REFFASTA << endl;
	cout << "output file: " << OUTFILE << endl;
	if(OUTFASTA != ""){
		cout << "output FASTA: " << OUTFASTA << endl;
	}
	cout << "base quality threshold: " << BASEQTHRESHOLD << endl;
	cout << "mapping quality threshold: " << MAPQTHRESHOLD << endl;
	cout << "left read trimming: " << LEFTTRIM << endl;
	cout << "right read trimming: " << RIGHTTRIM << endl;
	cout << endl;
}

// PrintUsage
// When called, prints the usage statement for this program.
void PrintUsage(){
	printf("\n\n");
	printf("Usage: SummarizeBAM -i input.sam -f ref.fasta -o out.summary\n");
	printf("\n");
	printf("Input options (defaults in parentheses):\n");
	printf("  -s FILE\twrite consensus sequence to FILE\n");
	printf("  -Q INT\tminimum base quality for a base to be tallied [20]\n");
	printf("  -q INT\tminimum mapping quality for a read to be tallied [20]\n");
	printf("  -l INT\tnum bases to trim from 5' (left) end of each read, after soft clipping [0]\n");
	printf("  -r INT\tnum bases to trim from 3' (right) end of each read, after soft clipping [0]\n");
	printf("\n\n");
}

// SetDebug
// Sets all parameters to their debug state.
void SetDebug(){
	SAM="sam.test";
	REFFASTA="ref.test";
	OUTFILE="out.test";
	OUTFASTA="seq.test";
}

// ReadMultiFasta
// Given a file name for a FASTA file containing multiple sequences,
// as well as a location to store multiple sequences,
// reads in the sequences as multiple strings to the given location
// and returns the FASTA header names as a vector of strings.
// Returns 1 if the file does not exist.
int ReadMultiFasta(string filename,
		vector<string> *sequencenames,
		vector<string> *sequences){

	// Open the file.
	ifstream f_in(filename.c_str(), ios::in);

	string line;
	string header;
	string sequence;
	int numsequences=0;

	if(f_in){
		// Read in the file line by line,
		// storing lines that begin with '>' as the sequence name
		while(getline(f_in, line)){

			if(line[0] == '>') {
				if(numsequences != 0){
					(*sequences).push_back(sequence);
				}
				numsequences += 1;
				sequence = "";
				// Store the first word of the sequence name,
				// with the > character removed.
				(*sequencenames).push_back(
						StringSplit(StringSplit(line,' ')[0],'>')[0]);
				continue;
			}
			sequence += line;
		}
		(*sequences).push_back(sequence);
	}
	else{
		return 1;
	}

	// Close the file.
	f_in.close();

	return 0;

}

//
// StringSplit
// Takes in a string and a character delimiter
// and returns a vector of strings split at that character.
vector<string> StringSplit(string s, char c){
	vector<string> splits;
	string s0;
	unsigned int i=0;

	while(i < s.length()){
		// Skip through delimiter characters at the beginnings of lines.
		while(s[i] == c && i < s.length() - 1){
			i++;
		}
		// Iterate through actual characters until you encounter c.
		while(i < s.length() && s[i] != c){
			s0 += s[i];
			i++;
		}
		// Once c is encountered, stop and save the string, then reset it.
		if(s0.size() > 0){
			splits.push_back(s0);
			s0 = "";
		}
		i++;
	}

	return splits;
}

//
// ExpandCIGAR
// Takes in a CIGAR string as specified by the SAM/BAM file standard
// and expands it into a string the length of the read describing its features.
// For instance, 1S5M4S becomes SMMMMMSSSS
string ExpandCIGAR(string cigar){
	string expanded="";
	string NumBases="";
	// All allowed operations characters.
	string CigarOperations="MIDNSHP=X";
	// Operations characters that sum to the length of the sequence.
	string CigarOperationsSubset="MIS=X";

	// Iterate through the CIGAR string and check for operation characters.
	for(unsigned int i=0; i<cigar.length(); i++){
		// Store non-operation characters that indicate the number of bases.
		if(strchr(CigarOperations.c_str(), cigar[i])==NULL){
			NumBases+=cigar[i];
		}
		// After reaching one of the useful operation characters,
		// add the appropriate number of characters to the expanded string.
		else if(strchr(CigarOperationsSubset.c_str(), cigar[i])!=NULL){
			unsigned int Num=atoi(NumBases.c_str());
			for(unsigned int j=0; j<Num; j++){
				expanded+=cigar[i];
			}
			NumBases="";
		}
		// If the operation character should not be added to the string,
		// reset the number of bases.
		else{
			NumBases="";
		}
	}
	return expanded;
}
