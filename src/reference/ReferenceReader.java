package reference;

import java.io.File;
import java.io.FileNotFoundException;

import dataStructure.ChromosomeOffset;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class ReferenceReader {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			IndexedFastaSequenceFile ifsf = new IndexedFastaSequenceFile(
					new File(
							"/home/dboehm/NewReference_2012_06_24/share/reference/genomes/NCBI_GRCh37/NCBI_GRCh37.fa"));
			System.out.println(ifsf.isIndexed());
			ReferenceSequence rs = ifsf.getSequence("chr1");
			SAMSequenceDictionary ssd = ifsf.getSequenceDictionary();
			SAMSequenceRecord ssr = ssd.getSequence(0);
			System.out.println(ssr.getSequenceLength());
			System.out.println(rs.getContigIndex() + " " + rs.length());
			System.out.println(rs.getName());
			// byte[] sequence = rs.getBases();
			int count = 0;
			// for (byte c : sequence) {
			// System.out.print((char) c);
			// if (count++ % 100 == 0) System.out.println();
			// }
			ReferenceSequence rf = ifsf
					.getSubsequenceAt("chr1", 100000, 100999);
			byte[] sequence = rf.getBases();
			for (byte c : sequence) {
				System.out.print((char) c);
				count++;
				if (count % 10 == 0)
					System.out.print(" ");
				if (count % 100 == 0)
					System.out.println();
			}

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
