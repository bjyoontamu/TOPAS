------------------------------------------------------------------------
		TOPAS: Network-based RNA structural alignment algorithm
------------------------------------------------------------------------


--------------------------Matlab implementation--------------------------

 Usage:
		[output_aseq1, output_aseq2]= TOPAS(input_sequences, input_prob1, input_prob2, input_prob3, input_parm1, input_parm2)

 
 Input:	
        input_sequences		the input sequences in FASTA format
      	input_prob1		the input probability of base pairing for sequence1 with format (base1 base2 probability)
      	input_prob2		the input probability of base pairing for sequence2 with format (base1 base2 probability)
      	input_prob3		the input probability of base alignment with format (base1 base2 probability)
      	input_parm1		the input parameter alpha for structure similarity
      	input_parm2		the input parameter beta for connected similarity

 Output:
        output_aseq1		the output aligned sequence 1
        output_aseq2		the output aligned sequence 2

 
 Example:
        [align1,align2]= TOPAS(‘seqs.fa’,’bp.1’,’bp.2’,’ap.12’, 0.5, 0.48)


 Input file format:
	The file input_sequences is to describe sequences in the FASTA format.
	The input_prob1 is to describe the base pairing probability in a format as follows:
               node1	node2	base-pairing probability

	The input_prob2 is to describe the base pairing probability with the same format as input_prob1.
	The input_prob3 is to describe the alignment probability in a format as follows:
               node1	node2	alignment probability

 Input parameters:
	The input_parm1 is weight parameter alpha to scale structure similarity. The range of alpha is 0≤ alpha≤ 1
	The input_parm2 is weight parameter alpha to scale connected similarity. The range of beta is 0≤ beta≤ 1



 For more information on the algorithm, please see:

 Chun-Chi Chen, Hyundoo Jeong, Xiaoning Qian and Byung-Jun Yoon, “TOPAS: network-based structural alignment of RNA sequences”

 Contact: aky3100@tamu.edu, bjyoon@ece.tamu.edu

 All Rights Reserved.
