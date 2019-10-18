//============================================================================
// Name        : AnnotateVariants.cpp
// Version     : 1.13
// Description : 1.13 Output position of a variant in a codon.
//          1.12 Output reference and alternate codons.
//          1.11 Fix problems arising from translating base N.
//           1.1 Allow multiple annotations for a single variant.
//               For instance, one might be in M1 and M2.
//           1.0 Given a BED format file, a reference, and list of sites,
//               annotate those sites as synonymous or nonsynonymous.
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
string VARFILE="";
string REFFASTA="";
string REFBED="";
string OUTFILE="";

bool DEBUG=false;

struct Annotation_t{
	string Chr="";
	int ChrStart=0;
	int ChrEnd=0;
	string Name="";
	int NumExons=0;
	vector <int> ExonSizes;
	vector <int> ExonStarts;
};


// FUNCTIONS
int ArgsParse(int argc, char *argv[]);
void PrintUsage();
void PrintParameters();
void SetDebug();
int ReadMultiFasta(string filename,
		vector<string> *sequencenames,
		vector<string> *sequences);
int ReadBED(string filename, vector<Annotation_t> *annotations);
vector<string> StringSplit(string s, char c);
int InitializeCodonTable(map<char, map<char, map<char,char> > > *codontable);
char TranslateCodon(string codonseq, map<char, map<char, map<char,char> > > *codontable);
char GenePosToBase(int genepos, string *refseqp, Annotation_t annotation);

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

	//==================================================
	// Read in BED annotation file.
	// Designed to follow specification here:
	// https://genome.ucsc.edu/FAQ/FAQformat.html
	//==================================================

	printf("Reading BED file.\n");
	vector<Annotation_t> Annotations;
	if(ReadBED(REFBED,&Annotations) != 0){
		printf("Error: BED annotation does not exist.\n");
		return 1;
	}

	//==================================================
	// Set up codon table to use for translations.
	//==================================================

	// Initialize a data structure to store the codon table.
	map<char, map<char, map<char,char> > > CodonTable;
	if(InitializeCodonTable(&CodonTable)!=0){
		printf("Error in codon table initialization.\n");
		return 1;
	}

	//==================================================
	// Read in variant file and annotate variants.
	//==================================================

	printf("Processing variant file.\n");

	// Open the variant file.
	ifstream fin(VARFILE.c_str(), ios::in);

	// Open the output file.
	ofstream fout(OUTFILE.c_str(), ios::out);

	// List the DNA bases.
	const char BASES[4]={'T','C','A','G'};


	if(fin){

		// Read in the file line by line.
		string line;

		while(getline(fin, line)){

			// Split the tab-delimited line.
			vector<string> fields=StringSplit(line,'\t');

			// Check that the appropriate fields are present.
			if(fields.size() < 4){
				printf("Variant file does not contain sufficient fields.\n");
				return 1;
			}

			// Store the information in the appropriate formats.
			// Fields in the variant file are hard-coded to follow a one-indexed
			// Chr Pos Base RefBase format.
			// Change to a zero-indexed format to accord with the BED specs.
			string Chr=fields[0];
			int Pos=atoi(fields[1].c_str())-1;
			string AltBase=fields[2];
			string RefBase=fields[3];

			// Iterate through the vector of annotations.
			// Look for matches to the specified variant.
			string Gene="none";
			int Codon=-1;
			int CodonIntPos=-1;
			char RefAA='Z';
			char AltAA='Z';
			int Syn=-1;
			int FourfoldSyn=-1;
			string RefCodonSeq="NNN";
			string AltCodonSeq="NNN";
			for(unsigned int i=0; i<Annotations.size(); i++){
				if(Chr!=Annotations[i].Chr){
					continue;
				}
				if(Pos>=Annotations[i].ChrStart &&
						Pos<Annotations[i].ChrEnd){

					// Iterate through the exons and determine whether
					// the variant falls in any of them.
					for(int j=0; j<Annotations[i].NumExons; j++){
						if(Pos>=Annotations[i].ExonStarts[j]+
								Annotations[i].ChrStart &&
								Pos<Annotations[i].ChrStart+
								Annotations[i].ExonStarts[j]+
								Annotations[i].ExonSizes[j]){

							Gene=Annotations[i].Name;

							// Determine the position of the base in the gene
							// and calculate its codon number.
							int GenePos=Pos-Annotations[i].ExonStarts[j]-
									Annotations[i].ChrStart;
							for(int k=0; k<j; k++){
								GenePos+=Annotations[i].ExonSizes[k];
							}
							Codon = (int) GenePos/3;
							CodonIntPos = GenePos % 3;

							// Determine the sequence of the reference codon.
							for(unsigned int k=0; k<RefNames.size(); k++){
								if(Chr==RefNames[k]){
									// For each position in the codon,
									// determine the base in the reference.
									int GeneCodonStart=GenePos-CodonIntPos;
									for(int m=0; m<3; m++){
										RefCodonSeq[m]=GenePosToBase(GeneCodonStart+m,
												&RefSequences[k], Annotations[i]);
									}

									// Implement a check to ensure that
									// the annotation is being read correctly.
									if(RefSequences[k][Pos]!=
											GenePosToBase(GenePos, &RefSequences[k],
													Annotations[i])){
										printf("Error in reading annotation.\n");
										return 1;
									}

									// Verify that the base in the reference
									// corresponds to the reference base
									// given in the variant file.
									if(RefSequences[k][Pos] != RefBase[0]){
										printf("Invalid reference base at position %d.\n", Pos);
										return 1;
									}
								}
							}

							// Translate the reference codon.
							RefAA=TranslateCodon(RefCodonSeq, &CodonTable);

							// Translate the alternate codon.
							AltCodonSeq=RefCodonSeq;
							AltCodonSeq[CodonIntPos]=AltBase[0];
							AltAA=TranslateCodon(AltCodonSeq, &CodonTable);

							// Determine whether the change was synonymous.
							Syn=(RefAA==AltAA) ? 1 : 0;

							// Determine whether the site is fourfold synonymous.
							string PossibleCodonSeq=RefCodonSeq;
							FourfoldSyn=1;
							for(int k=0; k<4; k++){
								PossibleCodonSeq[CodonIntPos]=BASES[k];
								if(TranslateCodon(PossibleCodonSeq, &CodonTable)!=RefAA){
									FourfoldSyn=0;
								}
							}
						}
					}
				}

				if(DEBUG){
					cout << Chr << "\t" << Pos+1 << "\t" <<
							AltBase << "\t" << RefBase << "\t" <<
							Gene << "\t" << Codon+1 << "\t" <<
							RefAA << "\t" << AltAA << "\t" <<
							Syn << "\t" << FourfoldSyn << "\t" <<
							RefCodonSeq << "\t" << AltCodonSeq << "\t" <<
							CodonIntPos << endl;
				}

				// Write the annotation information to a new file.
				// The codon number is 1-indexed.
				fout << line << "\t" << Gene << "\t" << Codon+1 << "\t" <<
						RefAA << "\t" << AltAA << "\t" <<
						Syn << "\t" << FourfoldSyn << "\t" <<
						RefCodonSeq << "\t" << AltCodonSeq << "\t" <<
						CodonIntPos << endl;

				// Iterate through the vector of annotations.
				// Reset gene annotations.
				Gene="none";
				Codon=-1;
				CodonIntPos=-1;
				RefAA='Z';
				AltAA='Z';
				Syn=-1;
				FourfoldSyn=-1;
				RefCodonSeq="NNN";
				AltCodonSeq="NNN";

			}


		}
	}
	else{
		printf("Error: Variant file does not exist.\n");
		return 1;
	}

	// Close the files.
	fin.close();
	fout.close();


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
		// -i input variant file
		case 'i':
			VARFILE = arg;
			break;
		// -r reference FASTA
		case 'f':
			REFFASTA = arg;
			break;
		// -b reference BED file
		case 'b':
			REFBED = arg;
			break;
		// -o output annotated variant file
		case 'o':
			OUTFILE = arg;
			break;
		}
	}

	// Check that the required arguments exist.
	if(VARFILE==""){
		printf("Invalid arguments. Specify variant file.\n");
		return 1;
	}
	if(REFFASTA==""){
		printf("Invalid arguments. Specify reference sequence.\n");
		return 1;
	}
	if(REFBED==""){
		printf("Invalid arguments. Specify reference BED annotation file.\n");
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
	cout << "input file: " << VARFILE << endl;
	cout << "reference: " << REFFASTA << endl;
	cout << "reference annotation: " << REFBED << endl;
	cout << "output file: " << OUTFILE << endl;
	cout << endl;
}

// PrintUsage
// When called, prints the usage statement for this program.
void PrintUsage(){
	printf("\n\n");
	printf("Usage: AnnotateVariants -i input.txt -f ref.fasta -b ref.bed -o out.txt\n");
	printf("\n");
	printf("Required inputs:\n");
	printf("  -i FILE\ttab-delimited variant file in form chr, pos, base, refbase\n");
	printf("  -f FILE\tFASTA format reference sequence\n");
	printf("  -b FILE\tBED format annotation file for reference sequence\n");
	printf("  -o FILE\toutput; tab-delimited variant file\n");
	printf("\n\n");
}

// SetDebug
// Sets all parameters to their debug state.
void SetDebug(){
	VARFILE="var.test";
	REFFASTA="ref.test";
	REFBED="bed.test";
	OUTFILE="out.test";
	DEBUG=true;
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

// ReadBED
// Given a file name for a BED format file containing sequence annotations,
// as well as a location to store the annotations,
// reads in the annotations and stores the information in appropriate form.
// Returns 1 if the file does not exist.
int ReadBED(string filename, vector<Annotation_t> *annotations){
	// Open the file.
	ifstream f_in(filename.c_str(), ios::in);

	string line;
	string header;
	string sequence;

	if(f_in){
		// Read in the file line by line,
		// storing lines that begin with '>' as the sequence name
		while(getline(f_in, line)){

			vector<string> fields=StringSplit(line,'\t');

			// Check that annotation contains the necessary fields.
			if(fields.size()>=12){

				// Verify that genes are on the positive strand of the vRNA.
				if(fields[5]!="+"){
					printf("This script does not accept negative-sense genes.\n");
					return 1;
				}

				// Store relevant information in Annotation_t format.
				Annotation_t annotation;

				annotation.Chr=fields[0];
				annotation.ChrStart=atoi(fields[1].c_str());
				annotation.ChrEnd=atoi(fields[2].c_str());
				annotation.Name=fields[3];
				annotation.NumExons=atoi(fields[9].c_str());

				vector<string> exonsizes=StringSplit(fields[10],',');
				vector<string> exonstarts=StringSplit(fields[11],',');
				for(int i=0; i<annotation.NumExons; i++){
					annotation.ExonSizes.push_back(atoi(exonsizes[i].c_str()));
					annotation.ExonStarts.push_back(atoi(exonstarts[i].c_str()));
				}

				(*annotations).push_back(annotation);
			}
			else{
				printf("BED file is not in accepted format.\n");
				return 1;
			}
		}
	}
	else{
		return 1;
	}

	// Close the file.
	f_in.close();

	return 0;
}

// InitializeCodonTable
// Given a series of nested maps to store a codon table,
// returns an initialized codon table.
int InitializeCodonTable(map<char, map<char, map<char,char> > > *codontable){

	const char BASES[4]={'T','C','A','G'};
	string CodonsToAAs="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRR"
			"IIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

	int iCodon=0;
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			for(int k=0; k<4; k++){
				(*codontable)[BASES[i]][BASES[j]][BASES[k]]=CodonsToAAs[iCodon];
				iCodon++;
			}
		}
	}
	return 0;
}


// TranslateCodon
// Given a three-base codon and a codon table,
// return the one-letter amino-acid code for the codon.
// Return X if the codon is invalid.
char TranslateCodon(string codonseq, map<char, map<char, map<char,char> > > *codontable){
	// Verify that the codon has three bases and contains no N's..
	if(codonseq.size()==3 && codonseq.find('N')==std::string::npos){
		return (*codontable)[codonseq[0]][codonseq[1]][codonseq[2]];
	}
	else{
		return 'X';
	}
}

// GenePosToBase
// Given a position in a gene (i.e. the three positions corresponding to a codon),
// as well as a reference sequence and the BED-format annotation for that sequence,
// return the base corresponding to that position in the gene.
char GenePosToBase(int genepos, string *refseqp, Annotation_t annotation){
	int ChrPos=0;
	int GenePosCounter=0;
	// Iterate through each exon and determine whether the desired gene position
	// is present.
	for(int i=0; i<annotation.NumExons; i++){

		// Offset the position in the chromosome based on the exon start site.
		ChrPos=annotation.ChrStart + annotation.ExonStarts[i];
		if(genepos>=GenePosCounter &&
				genepos<GenePosCounter+annotation.ExonSizes[i]){
			ChrPos+=genepos-GenePosCounter;
			break;
		}
		// If not, then increment the counter for positions in the gene.
		else{
			GenePosCounter+=annotation.ExonSizes[i];
		}
	}
	return (*refseqp)[ChrPos];
}
