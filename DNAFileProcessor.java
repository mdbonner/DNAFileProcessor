/********************************************************************************
 * Written by: Mark Bonner														*
 * 																				*
 * This program will perform various computations on genome sequence 			*
 * files. The genome input file must be of FASTA format, with the accessionID	*
 * at the top line and the entire genome sequence following immediately			*
 * on the next line. Output files will be created for complementation,			*
 * transcription, and translation computations.									*
 * 																				*
 * Computations performed:														*
 * 1. Count individual DNA nucleotides 											*
 * 2. Compute and return the GC content of the genome							*
 * 3. Complement and reverse a genome											*
 * 4. Transcribe DNA sequence to RNA											*
 * 5. Translate RNA sequence to proteins										*
 * 6. Calculate the protein mass of the protein sequence						*
 * 7. Count Hamming Distance (point mutations) between two DNA sequences		*
 * 																				*
 * This program is continually being updated with new computations.				*
 * 																				*
 *******************************************************************************/
 
import java.util.*;
import java.io.*;

public class DNAFileProcessor {

	public static void main(String[] args) throws FileNotFoundException {
		Scanner console = new Scanner(System.in);
		prompt();
		
		// Get file FASTA file contents
		HashMap<String, String> info = new HashMap<>(getFile(console));
		
		String fileName = info.get("fileName");
		String accessionID = info.get("ID");
		String dna = info.get("genome");
		
		System.out.println("\nComputations performed on " + accessionID);
		
		// Counting DNA Nucleotides
		System.out.println("\nNucleotide count in genome: ");
		countNucleotides(accessionID, dna);
		
		// Computing GC Content
		System.out.println("\nGuanine/Cytosine content of the DNA: ");
		computingGCContent(dna);
		
		// Complementing a strand of DNA
		System.out.println("\nReverse Complemented DNA: ");
		reverseComplement(fileName, accessionID, dna);	
		
		// Transcribe DNA into RNA
		System.out.println("\nDNA Transcribed to RNA: ");
		transcription(fileName, accessionID, dna);
		
		// Translating RNA into Protein
		System.out.println("\nRNA Sequence Translated into Proteins: ");
		String rna = transcription(dna);
		String protein = translateRNAIntoProtein(fileName, accessionID, rna);
		
		// Calculating Protein Mass
		System.out.println("\nTotal monoisotopic mass of protein string: ");
		proteinMass(protein);
		
		// Counting Point Mutations
		System.out.println("\n\nEnter y/n if you would like to count point mutations "
				+ "between two DNA sequences of equal length: ");
		if (console.nextLine().equalsIgnoreCase("y")) {
			countingPointMutations(dna, accessionID, console);
		}
		
	}
	
	// Get FASTA file input
	public static HashMap<String, String> getFile(Scanner console) throws FileNotFoundException {
		HashMap<String, String> info = new HashMap<>();
		
		// Ask the user for the name of the file they wish to process
		String inputName = console.nextLine();
		// Open the file by appending .fasta to the postfix
		Scanner fileScanner = new Scanner(new File(inputName + ".fasta"));
		// Copy accessionID from genome sequence
		String accessionID = fileScanner.nextLine();	

		// Copy entire DNA sequence into a variable
		String dna = "";
		while (fileScanner.hasNextLine()) {
			dna += fileScanner.nextLine();
		}

		info.put("fileName", inputName);
		info.put("ID", accessionID);
		info.put("genome", dna);

		return info;
	}
	
	/* Given: A protein string P.
	 * Return: The total weight of P.
	 */
	public static void proteinMass(String p) throws FileNotFoundException {
		Scanner protein = new Scanner(p);
		char[] proteinString = protein.next().toCharArray();
		
		// constructing the data structure for the Monoisotopic Mass Table
		HashMap<Character, Double> massTable = buildMassTable();
		
		double weight = 0.0;
		
		/* iterate through each amino acid in the protein string and get the
		 * monoisotopic mass for that amino acid from the data structure
		 */
		for (int i = 0; i < proteinString.length; i++) {
			weight += massTable.getOrDefault(proteinString[i], 0.0);
		}
		System.out.printf("%.3f", weight);
		System.out.print(" Da");
		protein.close();
	}
	
	/* method to create a data structure for each amino acid
	 * and it's corresponding monoisotopic mass
	 */
	public static HashMap<Character, Double> buildMassTable() throws FileNotFoundException {
		Scanner file = new Scanner(new File("monoisotopic_mass_table.txt"));
		HashMap<Character, Double> massTable = new HashMap<>();
		
		while (file.hasNext()) {
			Character aminoAcid = file.next().charAt(0);
			Double mass = file.nextDouble();
			massTable.put(aminoAcid, mass);
		}
		return massTable;
	}	
	
	/* Given: Two DNA strings s and t of equal length.
	 * Return: The Hamming distance dH(s,t).
	 */
	public static void countingPointMutations(String dna1, String ID1, Scanner console) throws FileNotFoundException {
		System.out.println("Input file name of second DNA sequence, following the guide above: ");

		// Get file FASTA file contents
		HashMap<String, String> info = new HashMap<>(getFile(console));

		//String fileName = info.get("fileName");
		String ID2 = info.get("ID");
		String dna2 = info.get("genome");
		int hammingDistance = 0;
		
		char[] array1 = new char[dna1.length()];
		char[] array2 = new char[dna2.length()];
		
		// copy each nucleobase from the DNA string into an array
		for (int i = 0; i < dna1.length(); i++) {
			array1[i] = dna1.charAt(i);
			array2[i] = dna2.charAt(i);
		}
		
		/* interate through each element of the arrays and compare the nucleobases,
		 * if the nucleobases match, then increment the hammingDistance variable
		 */
		for (int n = 0; n < array1.length; n++) {
			if (array1[n] != array2[n]) {
				hammingDistance++;
			}
		}		
		System.out.println("\nPoint Mutations between given DNA sequences: ");
		System.out.println(ID1 + "\n" + ID2);
		System.out.println(hammingDistance);		
	}
	
	/* Given: A DNA string dna.
	 * Return: The GC-content of the given dna sequence.
	 */
	public static void computingGCContent(String dna) {
		
		/* count each guanine and cytosine nucleobase in the DNA string and 
		 * print out the GC-content percentage
		 */
		int count = 0;
		int length = dna.length();

		for (int i = 0; i < length - 1; i++) {
			if (dna.charAt(i) == 'G' || dna.charAt(i) == 'C') {
				count++;
			}
		}
		double GCPercentage = (double) count / length * 100;
		System.out.printf("%.6f", GCPercentage);
		System.out.println("%");
	}
	
	/* Given: A DNA string dna.
	 * Return: The reverse complement dna^c of dna.
	 */
	public static void reverseComplement(String fileName, String id, String dna) throws FileNotFoundException {
		PrintStream outputFile = new PrintStream(new File(fileName + "_OUTPUT_COMPLEMENT.txt"));

		// building a data structure for the complement of each DNA symbol
		HashMap<Character, Character> complement = new HashMap<>();
		complement.put('A', 'T');
		complement.put('T', 'A');
		complement.put('C', 'G');
		complement.put('G', 'C');

		String reverseComplement = "";

		/* begin the for-loop from the end of the given DNA string 
		 * and complement each symbol into the reverse complement DNA string.		 * 
		 */
		for (int i = dna.length() - 1; i >= 0; i--) {
			reverseComplement += complement.get(dna.charAt(i));
		}
		
		outputFile.print(id + "\n" + reverseComplement);
		
		if (reverseComplement.length() < 101) {
			System.out.println("***output file created***");
			System.out.println(reverseComplement);
		} else {
			System.out.println("***output file created***");
		}
	}
	
	/* Given: A DNA string dna.
	 * Return: The transcribed RNA string of given dna.
	 */
	public static void transcription(String fileName, String id, String dna) throws FileNotFoundException {
		PrintStream outputFile = new PrintStream(new File(fileName + "_OUTPUT_TRANSCRIPTION.txt"));
		
		String rna = dna.replace('T', 'U');
		outputFile.print(id + "\n" + rna);
		
		if (rna.length() < 101) {
			System.out.println("***output file created***");
			System.out.println(rna);
		} else {
			System.out.println("***output file created***");
		}
	}
	
	public static String transcription(String dna) {
		String rna = dna.replace('T', 'U');
		return rna;
	}
	
	/* Given: An RNA string rna corresponding to a strand of mRNA.
	 * Return: The protein string encoded by rna.
	 */
	public static String translateRNAIntoProtein(String fileName, String id, String rna) throws FileNotFoundException {
		PrintStream outputFile = new PrintStream(new File(fileName + "_OUTPUT_RNA_TO_PROTEIN.txt"));
		HashMap<String, String> rnaCodonTable = buildRNACodonStructure();

		String[] arr = new String[rna.length() / 3];

		// copy all codons in the RNA string into an array
		for (int i = 0, start = 0, end = 3; i < arr.length; i++, start += 3, end += 3) {
			arr[i] = rna.substring(start, end);
		}

		String proteinString = "";

		// use the RNA Codon Table data structure to encode each codon into it's corresponding amino acid
		for (int i = 0; i < arr.length; i++) {
			proteinString += rnaCodonTable.get(arr[i]);
		}
		
		outputFile.print(id + "\n" + proteinString);
		outputFile.close();
		
		if (proteinString.length() < 101) {
			System.out.println("***output file created***");
			System.out.print(proteinString + "\n");
		} else {
			System.out.println("***output file created***");
		}
		return proteinString;
	}
	
	/* this method uses an input file containing the amino acid alphabet
	 * and it's corresponding codons to build a data structure for encoding
	 */
	public static HashMap<String, String> buildRNACodonStructure() throws FileNotFoundException {
		Scanner table = new Scanner(new File("RNA_codon_table.txt"));

		HashMap<String, String> rnaCodonTable = new HashMap<String, String>();
		
		while (table.hasNext()) {
			String rnaString = table.next();
			String aminoAcid = table.next();
			rnaCodonTable.put(rnaString, aminoAcid);
		}
		return rnaCodonTable;
	}
	
	/* Given: A DNA string dna.
	 * Return: Four integer counts for the respective number of times 
	 * 		   that the symbols 'A', 'C', 'G', and 'T' occur in dna.
	 */
	public static void countNucleotides(String id, String dna) {
		int[] array = new int[4];

		/* loop through each char in the DNA string and 
		 * increment the corresponding element of the array
		 */
		for (int i = 0; i < dna.length(); i++) {
			if (dna.charAt(i) == 'A') {
				array[0]++;
			} else if (dna.charAt(i) == 'C') {
				array[1]++;
			} else if (dna.charAt(i) == 'G') {
				array[2]++;
			} else {
				array[3]++;
			}
		}
		System.out.println(array[0] + " Adenine");
		System.out.println(array[1] + " Cytosine");
		System.out.println(array[2] + " Guanine");
		System.out.println(array[3] + " Thymine");
	}
	
	// Initial prompt to the user with instructions
	public static void prompt() {
		System.out.println("This program will perform various computations on a genome sequence of FASTA format.\n");
		System.out.println("Please follow these steps before continuing:\n");
		System.out.println("1. Move the FASTA file you wish to process into the same folder as this program.\n");
		System.out.println("2. Copy only the name of the FASTA file--be sure to exclude the .fasta file extension.\n");
		System.out.println("3. Paste the name of the file you wish to process "
				+ "then press enter (use CTRL + V to paste if right-click does not work): ");
	}
}
